from scripts.utils import *
import pandas as pd
import os
    
    
def get_all_structures_and_sequences(seed: list, bit_scores: list[float], number_tracker: dict, root_dir: str) -> dict:
    
    # Get annotation files
    annotation_files = get_annotation_files(root_dir, seed, bit_scores)
    number_tracker["annotation_files_jackhmmer"] = len(annotation_files)
    
    # Finding all unique IDs across annotation files
    unique_ids = get_unique_ids_from_annotation_files(annotation_files)
    number_tracker["unique_ids_jackhmmer"] = len(unique_ids)
    
    # Get unique IDs from foldseek_webserver_output
    foldseek_webserver_output_path = f"{root_dir}/data/filtering/foldseek/foldseek_webserver_output.csv"
    foldseek_webserver_ids_df = pd.read_csv(foldseek_webserver_output_path)[["target"]].drop_duplicates()
    
    # Combining all unique IDs
    foldseek_webserver_ids = set(foldseek_webserver_ids_df["target"].tolist())
    combined_unique_ids = set(unique_ids).union(foldseek_webserver_ids)
    unique_ids = list(combined_unique_ids)
    number_tracker["unique_ids_combined"] = len(unique_ids)    
    
    # 2. Extract uniprot ids from IDs for AlphaFold downloads
    uids_set = set()
    for id in unique_ids:
        uid = extract_uid_from_recordID(id)
        if uid:
            uids_set.add(uid)
    uids = list(uids_set)
    print(f"Total unique uniprot UIDs for AlphaFold: {len(uids)}/{len(unique_ids)}\n")
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
    """Gets the unique IDs from all annotation files in the specified list.
    Args:
        annotation_files (list[str]): A list of paths to annotation CSV files.

    Returns:
        list[str]: A list of unique IDs.
    """
    all_ids = set()
    for file_path in annotation_files:
        ids = pd.read_csv(file_path)["id"].dropna().drop_duplicates().to_list()
        all_ids.update(ids)
    unique_ids = list(all_ids)
    return unique_ids
