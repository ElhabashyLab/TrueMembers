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
        total=2,
        backoff_factor=0.5,
        status_forcelist=[429, 500, 502, 503, 504],
    )
    try:
        retries = Retry(**retry_kwargs, allowed_methods=frozenset(["GET", "POST"]))
    except TypeError:
        retries = Retry(**retry_kwargs, method_whitelist=frozenset(["GET", "POST"]))
    adapter = HTTPAdapter(pool_connections=pool_maxsize, pool_maxsize=pool_maxsize, max_retries=retries)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update({"User-Agent": "True-Members-Pipeline/1.0 (contact: aron.wichtner@tuebingen.mpg.de)"})
    return session


def _build_uniprot_url(uid: str) -> str:
    if uid.startswith("UPI"):
        return f"https://rest.uniprot.org/uniparc/{uid}.fasta"
    return f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"


def _fetch_fasta_text(uid: str, session: requests.Session) -> str:
    url = _build_uniprot_url(uid)
    with session.get(url, timeout=(5, 15)) as response:
        if not (response.status_code == 200 and response.text.startswith(">")):
            raise ValueError(f"Sequence not found for UID: {uid}")
        return response.text


def download_fasta_record(uid: str, session: requests.Session):
    text = _fetch_fasta_text(uid, session)
    return next(SeqIO.parse(StringIO(text), "fasta"))


def download_fasta_to_cache(uid: str, cache_dir: str, session: requests.Session):
    save_path = os.path.join(cache_dir, f"{uid}.fasta")
    tmp_path = save_path + f".tmp.{os.getpid()}"
    os.makedirs(cache_dir, exist_ok=True)

    text = _fetch_fasta_text(uid, session)
    with open(tmp_path, 'w') as f:
        f.write(text)
    os.rename(tmp_path, save_path)

    return next(SeqIO.parse(StringIO(text), "fasta"))


def _chunked(seq: list, size: int):
    for i in range(0, len(seq), size):
        yield seq[i:i + size]


def _build_batch_query(ids: list[str], is_uniparc: bool) -> str:
    field = "upi" if is_uniparc else "accession_id"
    return " OR ".join(f"{field}:{uid}" for uid in ids)


def _fetch_batch_fasta_text(ids: list[str], session: requests.Session, is_uniparc: bool) -> str:
    base_url = "https://rest.uniprot.org/uniparc/stream" if is_uniparc else "https://rest.uniprot.org/uniprotkb/stream"
    query = _build_batch_query(ids, is_uniparc)
    with session.get(base_url, params={"query": query, "format": "fasta"}, timeout=(5, 60)) as response:
        response.raise_for_status()
        return response.text


def download_fasta_batch_to_cache(ids: list[str], cache_dir: str | None, session: requests.Session, is_uniparc: bool):
    text = _fetch_batch_fasta_text(ids, session, is_uniparc)
    records = list(SeqIO.parse(StringIO(text), "fasta"))

    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
        for record in records:
            uid = extract_uid_from_recordID(record.id) or record.id
            save_path = os.path.join(cache_dir, f"{uid}.fasta")
            tmp_path = save_path + f".tmp.{os.getpid()}"
            SeqIO.write([record], tmp_path, "fasta")
            os.rename(tmp_path, save_path)

    return records


def _download_records_batched(accessions: list[str], max_workers: int, cache_dir: str | None,
                                initial_count: int, total_count: int, batch_size: int):
    if not accessions:
        return [], []

    uniparc_ids = [a for a in accessions if a and a.startswith("UPI")]
    uniprotkb_ids = [a for a in accessions if a and not a.startswith("UPI")]

    chunks = [(chunk, False) for chunk in _chunked(uniprotkb_ids, batch_size)]
    chunks += [(chunk, True) for chunk in _chunked(uniparc_ids, batch_size)]

    records: list = []
    found_ids: set = set()
    success = initial_count

    with _build_session(pool_maxsize=max_workers) as session:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_chunk = {
                executor.submit(download_fasta_batch_to_cache, chunk, cache_dir, session, is_uniparc): chunk
                for chunk, is_uniparc in chunks
            }
            for future in as_completed(future_to_chunk):
                chunk = future_to_chunk[future]
                try:
                    chunk_records = future.result()
                    records.extend(chunk_records)
                    for record in chunk_records:
                        found_ids.add(extract_uid_from_recordID(record.id) or record.id)
                    success += len(chunk_records)
                    logger.info("Fetched %d/%d records (batch of %d accessions)...", success, total_count, len(chunk))
                except Exception as e:
                    logger.error("Error fetching batch of %d accessions: %s", len(chunk), e)

    missing = [a for a in accessions if a not in found_ids]
    return records, missing


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


def _download_from_local_db(accessions: list[str], local_fasta_path: str, cache_dir: str | None, index_dir: str = "") -> tuple[list, list[str]]:
    if not accessions:
        return [], []

    from src.scripts.uniprot_local import UniprotLocalDB

    index_dir = index_dir or os.path.join(os.getcwd(), ".uniprot_index")
    db = UniprotLocalDB(local_fasta_path, index_dir=index_dir)
    db.ensure_index()

    records: list = []
    missing: list[str] = []
    for acc in accessions:
        record = db.get_sequence(acc)
        if record is None:
            missing.append(acc)
            continue
        records.append(record)
        if cache_dir:
            save_path = os.path.join(cache_dir, f"{acc}.fasta")
            tmp_path = save_path + f".tmp.{os.getpid()}"
            os.makedirs(cache_dir, exist_ok=True)
            SeqIO.write([record], tmp_path, "fasta")
            os.rename(tmp_path, save_path)

    logger.info("Local DB: %d found, %d missing from %s", len(records), len(missing), local_fasta_path)
    return records, missing


def get_records_from_accessions(accessions, max_workers: int = 8, cache_dir: str | None = None,
                                 batch_size: int = 100, use_batching: bool = True,
                                 local_fasta_path: str = "", local_fasta_index_dir: str = ""):
    logger.info("Downloading %d records from Uniprot/Uniparc...", len(accessions))

    if cache_dir:
        cached_records, uncached_accessions = _collect_cached_records(accessions, cache_dir)
        if uncached_accessions:
            logger.info("%d cached, downloading %d remaining.", len(cached_records), len(uncached_accessions))
    else:
        cached_records = []
        uncached_accessions = list(accessions)

    if not uncached_accessions:
        return cached_records

    local_records: list = []
    still_missing = uncached_accessions
    if local_fasta_path:
        local_records, still_missing = _download_from_local_db(
            uncached_accessions, local_fasta_path, cache_dir, local_fasta_index_dir
        )

    if not still_missing:
        return cached_records + local_records

    if use_batching:
        total = len(accessions)
        batched_records, missing = _download_records_batched(
            still_missing, max_workers, cache_dir, len(cached_records) + len(local_records), total, batch_size
        )
        if missing:
            logger.info("%d/%d accessions not returned by first batch pass; retrying once.",
                         len(missing), len(still_missing))
            retried, still_missing2 = _download_records_batched(
                missing, max_workers, cache_dir, len(cached_records) + len(local_records) + len(batched_records), total, batch_size
            )
            batched_records += retried
            if still_missing2:
                logger.info("%d accessions still not resolved after batch retry — skipped.", len(still_missing2))
        return cached_records + local_records + batched_records

    downloaded_records = _download_records_parallel(still_missing, max_workers, cache_dir, cached_records)
    return cached_records + local_records + downloaded_records


def get_records_from_fasta(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    return records


def write_fasta_from_records(records, out_path):
    SeqIO.write([record for record in records if record is not None], out_path, "fasta")
    return out_path


def run_download_tasks(
    uids: list[str],
    output_dir: str,
    max_workers: int = 8,
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
    tmp_path = save_path + f".tmp.{os.getpid()}"
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