import logging
import os
import re
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from io import StringIO

import requests
import yaml
from Bio import SeqIO
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import Select
from Bio.PDB.mmcifio import MMCIFIO
from requests.adapters import HTTPAdapter, Retry

logger = logging.getLogger(__name__)


def _log_interval(total: int) -> int:
    """Return how often to log progress: every 10 items or every 20%, whichever is smaller."""
    return max(1, min(10, total // 5))


def _build_session(pool_maxsize: int) -> requests.Session:
    """Builds a single reusable Session with connection pooling and retry/backoff."""
    session = requests.Session()
    retry_kwargs = dict(
        total=3,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
    )
    try:
        retries = Retry(**retry_kwargs, allowed_methods=frozenset(["GET", "POST"]))
    except TypeError:
        retries = Retry(**retry_kwargs, method_whitelist=frozenset(["GET", "POST"]))
    adapter = HTTPAdapter(pool_connections=1, pool_maxsize=pool_maxsize, max_retries=retries)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def _build_uniprot_url(uid: str) -> str:
    if uid.startswith("UPI"):
        return f"https://rest.uniprot.org/uniparc/{uid}.fasta"
    return f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"


def _fetch_fasta_text(uid: str, session: requests.Session) -> str:
    url = _build_uniprot_url(uid)
    with session.get(url, timeout=(5, 30)) as response:
        if not (response.status_code == 200 and response.text.startswith(">")):
            raise ValueError(f"Sequence not found for UID: {uid}")
        return response.text


def download_fasta_record(uid: str, session: requests.Session):
    text = _fetch_fasta_text(uid, session)
    return next(SeqIO.parse(StringIO(text), "fasta"))


def download_fasta_to_cache(uid: str, cache_dir: str, session: requests.Session):
    save_path = os.path.join(cache_dir, f"{uid}.fasta")
    tmp_path = save_path + ".tmp"
    os.makedirs(cache_dir, exist_ok=True)

    text = _fetch_fasta_text(uid, session)
    with open(tmp_path, 'w') as f:
        f.write(text)
    os.rename(tmp_path, save_path)

    return next(SeqIO.parse(StringIO(text), "fasta"))


def _collect_cached_records(accessions, cache_dir: str) -> tuple[list, list[str]]:
    records = []
    uncached: list[str] = []
    for acc in accessions:
        cache_path = os.path.join(cache_dir, f"{acc}.fasta")
        if os.path.exists(cache_path):
            cached = get_records_from_fasta(cache_path)
            if cached:
                records.append(cached[0])
                continue
        uncached.append(acc)
    return records, uncached


def _download_records_parallel(accessions, max_workers: int, cache_dir: str | None, initial_records: list):
    if not accessions:
        return []

    with _build_session(pool_maxsize=max_workers) as session:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            if cache_dir:
                future_to_acc = {
                    executor.submit(download_fasta_to_cache, acc, cache_dir, session): acc
                    for acc in accessions
                }
            else:
                future_to_acc = {
                    executor.submit(download_fasta_record, acc, session): acc
                    for acc in accessions
                }
            records: list = []
            target_total = len(initial_records) + len(accessions)
            success = len(initial_records)
            for future in as_completed(future_to_acc):
                acc = future_to_acc[future]
                try:
                    record = future.result()
                    records.append(record)
                    success += 1
                    if success % 100 == 0 or success == target_total:
                        logger.info("Fetched %d/%d records ...", success, target_total)
                except Exception as e:
                    logger.error("Error fetching %s: %s", acc, e)
    return records


def get_records_from_accessions(accessions, max_workers: int = 8, cache_dir: str | None = None):
    logger.info("Downloading %d records from Uniprot/Uniparc...", len(accessions))

    if cache_dir:
        cached_records, uncached_accessions = _collect_cached_records(accessions, cache_dir)
        if uncached_accessions:
            logger.info("%d cached, downloading %d remaining.", len(cached_records), len(uncached_accessions))
    else:
        cached_records = []
        uncached_accessions = list(accessions)

    downloaded_records = _download_records_parallel(uncached_accessions, max_workers, cache_dir, cached_records)
    return cached_records + downloaded_records


def get_records_from_fasta(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    return records


def write_fasta_from_records(records, out_path):
    SeqIO.write([record for record in records if record is not None], out_path, "fasta")
    return out_path


def run_download_tasks(
    uids: list[str],
    output_dir: str,
    max_workers: int = 16,
) -> list[str]:
    total = len(uids)
    interval = _log_interval(total)
    successful_uids = []

    with _build_session(pool_maxsize=max_workers) as session:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_task = {
                executor.submit(download_alphafold_structure, uid, output_dir, session): uid
                for uid in uids
            }
            completed = 0
            for future in as_completed(future_to_task):
                uid = future_to_task[future]
                completed += 1
                if future.result():
                    successful_uids.append(uid)
                else:
                    logger.warning("Skipping entry '%s': could not download AlphaFold structure.", uid)
                if completed % interval == 0 or completed == total:
                    logger.info("Processing %d/%d structures...", completed, total)

    failed = total - len(successful_uids)
    logger.info("Downloaded %d/%d structure files (%d failed).", len(successful_uids), total, failed)
    return successful_uids


_UNIPROT_ACCESSION_RE = re.compile(
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}"
)


def extract_uid_from_recordID(record_id):
    match = _UNIPROT_ACCESSION_RE.search(record_id)
    if match:
        return match.group()
    else:
        logger.warning("Could not extract UniProt accession from record ID: %s", record_id)
        return None


def extract_region_from_recordID(record_id):
    region_pattern = re.compile(r"/(\d+)-(\d+)")
    match = region_pattern.search(record_id)
    if match:
        return int(match.group(1)), int(match.group(2))
    else:
        logger.warning("Could not extract region from record ID: %s", record_id)
        return None, None


def download_alphafold_structure(
    accession_id: str,
    output_dir: str,
    session: requests.Session,
) -> bool:
    save_path = os.path.join(output_dir, f"{accession_id}.cif")
    if os.path.exists(save_path):
        logger.info("File already exists for '%s', skipping download.", accession_id)
        return True
    cif_url = _fetch_alphafold_cif_url(accession_id, session)
    if not cif_url:
        return False
    tmp_path = save_path + ".tmp"
    if not download_file(cif_url, tmp_path, session):
        return False
    os.rename(tmp_path, save_path)
    return True


def _fetch_alphafold_cif_url(accession_id: str, session: requests.Session):
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession_id}"
    try:
        with session.get(api_url, timeout=30) as response:
            response.raise_for_status()
            predictions = response.json()
    except requests.RequestException as e:
        logger.error("Error fetching AlphaFold CIF URL for %s: %s", accession_id, e)
        return False

    if not predictions:
        return False
    url = predictions[0].get("cifUrl")
    if not url:
        return False
    return url


def read_yaml_config(config_path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config


def download_file(url, output_path, session: requests.Session):
    try:
        with session.get(url, stream=True, timeout=10) as response:
            response.raise_for_status()
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            with open(output_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
        return True
    except Exception as e:
        logger.error("Failed to download %s: %s", url, e)
        return False


def copy_file(src: str, dst: str) -> bool:
    try:
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy(src, dst)
        return True
    except Exception as e:
        logger.error("Failed to copy %s to %s: %s", src, dst, e)
        return False


def filter_fasta_file(fasta_path: str, ids: list, out_path: str) -> str:
    records = get_records_from_fasta(fasta_path)
    filtered_records = []
    found_records = 0

    for id in ids:
        id_found = False
        for record in records:
            if extract_uid_from_recordID(record.id) == id:
                filtered_records.append(record)
                found_records += 1
                id_found = True
                break
        if not id_found:
            logger.warning("ID '%s' not found in the FASTA file: %s", id, fasta_path)
    logger.info("Found %d/%d ids in the FASTA file: %s", found_records, len(ids), fasta_path)
    write_fasta_from_records(filtered_records, out_path)
    return out_path


def extract_region_of_protein(path_to_protein: str, region: tuple[int, int], output_path: str | None = None) -> str:
    logger.info("Extracting region %s from protein file: %s", region, path_to_protein)

    class SelectRegion(Select):
        def accept_residue(self, residue):
            return region[0] <= residue.id[1] <= region[1]

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("protein", path_to_protein)
    if output_path is None:
        root, ext = os.path.splitext(path_to_protein)
        output_path = f"{root}_roi{ext}"
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(output_path, select=SelectRegion())
    return output_path