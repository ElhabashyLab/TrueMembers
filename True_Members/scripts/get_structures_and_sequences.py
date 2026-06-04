from scripts.utils import *
import pandas as pd
import os
    
    
def get_all_structures_and_sequences(seed: list, bit_scores: list[float], number_tracker: dict, root_dir: str) -> dict:
    
    # 0. Get annotation files
    annotation_files = get_annotation_files(root_dir, seed, bit_scores)
    number_tracker["annotation_files"] = len(annotation_files)
    
    # 1. Finding all unique RepIDs across annotation files
    unique_rep_ids = get_unique_ids_from_annotation_files(annotation_files)
    number_tracker["unique_rep_ids"] = len(unique_rep_ids)
    
    # 2. Extract uniprot ids from RepIDs for AlphaFold downloads
    uids_set = set()
    for id in unique_rep_ids:
        uid = extract_uid_from_recordID(id)
        if uid:
            uids_set.add(uid)
    uids = list(uids_set)
    print(f"Total unique uniprot UIDs for AlphaFold: {len(uids)}/{len(unique_rep_ids)}\n")
    number_tracker["unique_uniprot_uids"] = len(uids)

    # 3. Download structures from AlphaFold
    structures_output_dir = f"{root_dir}/data/alphafold_structures"
    succesful_uids = run_download_tasks(uids, structures_output_dir)
    number_tracker["successful_structure_downloads"] = len(succesful_uids)
    
    # 4. Download all records which have a successfully downloaded structure
    records = get_records_from_accessions(succesful_uids)
    number_tracker["downloaded_records"] = len(records)
    
    # 5. Write records to fasta
    fasta_path = write_fasta_from_records(records, f"{root_dir}/data/all_sequences_full.fasta")
    
    return number_tracker
    
    
def get_annotation_files(root_dir: str, seed: list, bit_scores: list[float]) -> list[str]:
    protein_uid = seed[0]
    r_start = seed[1]
    r_end = seed[2]
    prefix = f"{protein_uid}_{r_start}-{r_end}"
    
    annotation_files = []
    for b in bit_scores:
        batch_key = f"_b{b}"
        batch_prefix = f"{prefix}{batch_key}"
        annotation_file = os.path.join(root_dir, f"data/evcouplings/output/{batch_prefix}/align/{batch_prefix}_annotation.csv")
        if os.path.exists(annotation_file):
            annotation_files.append(annotation_file)
        else:
            print(f"Warning: Annotation file {annotation_file} does not exist.")
    
    return annotation_files


def get_unique_ids_from_annotation_files(annotation_files: list[str]) -> list[str]:
    """Gets the unique RepIDs from all annotation files in the specified list.
    Args:
        annotation_files (list[str]): A list of paths to annotation CSV files.

    Returns:
        list[str]: A list of unique RepIDs.
    """
    all_rep_ids = set()
    for file_path in annotation_files:
        rep_ids = pd.read_csv(file_path)["RepID"].dropna().drop_duplicates().to_list()
        all_rep_ids.update(rep_ids)
    unique_rep_ids = list(all_rep_ids)
    return unique_rep_ids
