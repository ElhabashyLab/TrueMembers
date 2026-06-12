from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import shutil
import requests
import yaml
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.PDBIO import Select
from io import StringIO
import re


def _log_interval(total: int) -> int:
    """Return how often to log progress: every 10 items or every 20%, whichever is smaller."""
    return max(1, min(10, total // 5))


def download_fasta_record(uid: str) -> SeqRecord:
    """
    Download a FASTA record by accession from Uniprot or Uniparc.
    """
    if uid.startswith("UPI"):
        url = f"https://rest.uniprot.org/uniparc/{uid}.fasta"
    else:
        url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
        
    response = requests.get(url, timeout=10)
    if not (response.status_code == 200 and response.text.startswith(">")):
        raise ValueError(f"Sequence not found for UID: {uid}")

    record = next(SeqIO.parse(StringIO(response.text), "fasta"))
    return record


def get_records_from_accessions(accessions):
    """
    Fetch FASTA records from Uniprot/Uniparc via its REST API concurrently.
    """
    print(f"Downloading {len(accessions)} records from Uniprot/Uniparc...")
    records = []
    
    with ThreadPoolExecutor(max_workers=16) as executor:
        future_to_acc = {executor.submit(download_fasta_record, acc): acc for acc in accessions}
        success = 0
        for future in as_completed(future_to_acc):
            acc = future_to_acc[future]
            try:
                record = future.result()
                records.append(record)
                success += 1
                if success % 100 == 0 or success == len(accessions):
                    print(f"Fetched {success}/{len(accessions)} records ...")
            except Exception as e:
                print(f"Error fetching {acc}: {e}")
                
    return records


def get_records_from_fasta(fasta_file):
    """
    Parse FASTA records from a local FASTA file.
    """
    records = list(SeqIO.parse(fasta_file, "fasta"))
    return records


def write_fasta_from_records(records, out_path):
    """Writes the sequences records of records to a fasta file."""
    SeqIO.write([record for record in records if record is not None], out_path, "fasta")
    return out_path


def run_download_tasks(
    uids: list[str],
    output_dir: str,
    max_workers: int = 16,
) -> list[str]:
    """Run download tasks concurrently."""
    total = len(uids)
    interval = _log_interval(total)
    successful_uids = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_task = {
            executor.submit(download_alphafold_structure, uid, output_dir): uid
            for uid in uids
        }
        completed = 0
        for future in as_completed(future_to_task):
            uid = future_to_task[future]
            completed += 1
            if future.result():
                successful_uids.append(uid)
            else:
                print(f"Skipping entry '{uid}': could not download AlphaFold structure.")
            if completed % interval == 0 or completed == total:
                print(f"Processing {completed}/{total} structures...")

    failed = total - len(successful_uids)
    print(f"Downloaded {len(successful_uids)}/{total} structure files ({failed} failed).")
    return successful_uids


# Matches the official UniProt accession format:
# https://www.uniprot.org/help/accession_numbers
_UNIPROT_ACCESSION_RE = re.compile(
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}"
)

def extract_uid_from_recordID(record_id):
    """
    Extracts the UniProt accession from a record ID.
    """
    match = _UNIPROT_ACCESSION_RE.search(record_id)
    if match:
        return match.group()
    else:
        print(f"Could not extract UniProt accession from record ID: {record_id}")
        return None
    
    
def extract_region_from_recordID(record_id):
    """
    Extracts the region (start, end) from a record ID in the format '{uid}/<start>-<end>'.
    """
    region_pattern = re.compile(r"/(\d+)-(\d+)")
    match = region_pattern.search(record_id)
    if match:
        return int(match.group(1)), int(match.group(2))
    else:
        print(f"Could not extract region from record ID: {record_id}")
        return None, None


def download_alphafold_structure(
    accession_id: str,
    output_dir: str
) -> bool:
    """
    Download an AlphaFold-predicted CIF structure.
    """
    cif_url = _fetch_alphafold_cif_url(accession_id)
    if not cif_url:
        return False
    save_path = os.path.join(output_dir, f"{accession_id}.cif")
    if not download_file(cif_url, save_path):
        return False
    return True


def _fetch_alphafold_cif_url(accession_id: str):
    """
    Queries the AlphaFold DB REST API and returns the CIF download URL.
    """
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession_id}"
    try:
        response = requests.get(api_url, timeout=30)
        response.raise_for_status()
        predictions = response.json()
    except requests.RequestException as e:
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


def download_file(url, output_path):
    """
    Downloads a file from a given URL and saves it to a specified local path.
    """
    try:
        response = requests.get(url, stream=True, timeout=10)
        response.raise_for_status()
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'wb') as f:
            # download file in 8KB chunks
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        return True
    except Exception as e:
        print(f"Failed to download {url}: {e}")
        return False
    

def copy_file(src: str, dst: str) -> bool:
    """
    Copies a file from src to dst.
    """
    try:
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy(src, dst)
        return True
    except Exception as e:
        print(f"Failed to copy {src} to {dst}: {e}")
        return False


def filter_fasta_file(fasta_path: str, ids: list, out_path: str) -> str:
    """Filters a FASTA file to include only records with IDs in the provided list.

    Args:
        fasta_path (str): Path to the input FASTA file.
        ids (list): List of IDs to filter by.
        out_path (str): Path to the output filtered FASTA file.

    Returns:
        str: Path to the output filtered FASTA file.
    """
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
            print(f"ID '{id}' not found in the FASTA file: {fasta_path}")
    print(f"Found {found_records}/{len(ids)} ids in the FASTA file: {fasta_path}")
    write_fasta_from_records(filtered_records, out_path)
    return out_path


def extract_region_of_protein(path_to_protein: str, region: tuple[int, int], output_path: str | None = None) -> str:
    """
    Extract a residue range from a CIF structure file.

    Parses a CIF structure, filters residues to the specified range, and writes
    the trimmed structure to an output file.

    Args:
        path_to_protein: Path to the input CIF structure file.
        region: Residue range (start, end) to extract, 1-based inclusive.
        output_path: Destination path. If None, appends '_roi' suffix to the input path.

    Returns:
        Path to the written output file.
    """
    print(f"Extracting region {region} from protein file: {path_to_protein}")
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
