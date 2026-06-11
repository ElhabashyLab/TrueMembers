from scripts.utils import *
import pandas as pd
import os
    
    
def get_all_structures_and_sequences(seed: list, bit_scores: list[float], number_tracker: dict, evcouplings_dir: str, foldseek_dir: str, alphafold_structures_dir: str, root_dir: str):
    unique_ids_with_regions, foldseek_webserver_ids_df = get_jackhmmer_and_foldseek_ids_with_regions(number_tracker, seed, bit_scores, evcouplings_dir, foldseek_dir)
    combined_unique_ids = combine_foldseek_and_jackhmmer_ids(unique_ids_with_regions, foldseek_webserver_ids_df)
    uids = list(combined_unique_ids["id"]) 
    number_tracker["unique_uniprot_uids_combined"] = len(uids)

    # Download structures from AlphaFold
    successful_uids = run_download_tasks(uids, alphafold_structures_dir)
    uids_with_region = combined_unique_ids[combined_unique_ids["id"].isin(successful_uids)] 
    number_tracker["successful_structure_downloads"] = len(successful_uids)
    
    # Download all records which have a successfully downloaded structure
    records = get_records_from_accessions(uids_with_region["id"].tolist())
    number_tracker["downloaded_records"] = len(records)
    
    # Keep only uids which have a uniprot record and a successfully downloaded structure
    record_ids = [extract_uid_from_recordID(record.id) for record in records]
    successful_uids_with_region = uids_with_region[uids_with_region["id"].isin(record_ids)]
    
    # Write records to fasta and csv with regions
    seed_containing_sequences_full_fasta_path = f"{root_dir}/data/seed_containing_sequences_full.fasta"
    seed_containing_sequences_with_region_path = f"{root_dir}/data/seed_containing_sequences_with_region.csv"
    seed_containing_sequences_full = write_fasta_from_records(records, seed_containing_sequences_full_fasta_path)
    successful_uids_with_region.to_csv(seed_containing_sequences_with_region_path, index=False)
    
    return seed_containing_sequences_full, successful_uids_with_region


def get_jackhmmer_and_foldseek_ids_with_regions(number_tracker: dict, seed: list, bit_scores: list[float], evcouplings_dir: str, foldseek_dir: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    annotation_files = get_annotation_files(evcouplings_dir, seed, bit_scores)
    number_tracker["annotation_files_jackhmmer"] = len(annotation_files)
    
    unique_ids_with_regions = get_unique_ids_with_regions_from_annotation_files(annotation_files)
    number_tracker["unique_ids_jackhmmer"] = len(unique_ids_with_regions)
    
    foldseek_webserver_output_path = f"{foldseek_dir}/foldseek_webserver_output.csv"
    foldseek_webserver_output = pd.read_csv(foldseek_webserver_output_path)
    foldseek_webserver_ids_df = foldseek_webserver_output[["id", "region_start", "region_end"]]
    return unique_ids_with_regions, foldseek_webserver_ids_df


def combine_foldseek_and_jackhmmer_ids(jackhmmer_ids_df: pd.DataFrame, foldseek_ids_df: pd.DataFrame) -> pd.DataFrame:    
    combined_unique_ids = (
        pd.concat([jackhmmer_ids_df, foldseek_ids_df])
        .groupby("id", as_index=False)
        .agg(
            region_start=("region_start", lambda x: round(x.mean())),
            region_end=("region_end", lambda x: round(x.mean())),
        )
        .astype({"region_start": "Int64", "region_end": "Int64"})
    )
    
    return combined_unique_ids
    
    
def get_annotation_files(evcouplings_dir: str, seed: list, bit_scores: list[float]) -> list[str]:
    protein_uid = seed[0]
    r_start = seed[1]
    r_end = seed[2]
    prefix = f"{protein_uid}_{r_start}-{r_end}"
    
    annotation_files = []
    for b in bit_scores:
        batch_key = f"_b{b}"
        batch_prefix = f"{prefix}{batch_key}"
        annotation_file = os.path.join(evcouplings_dir, f"output/{batch_prefix}/align/{batch_prefix}_annotation.csv")
        if os.path.exists(annotation_file):
            annotation_files.append(annotation_file)
        else:
            print(f"Warning: Annotation file {annotation_file} does not exist.")
    
    return annotation_files


def get_unique_ids_with_regions_from_annotation_files(annotation_files: list[str]) -> pd.DataFrame:
    """Gets the unique IDs with regions from all annotation files in the specified list.
    Args:
        annotation_files (list[str]): A list of paths to annotation CSV files.

    Returns:
        pd.DataFrame: A DataFrame containing unique IDs with regions (id, region_start, region_end).
    """
    dfs_ids_with_regions = []
    for file_path in annotation_files:
        ids = pd.read_csv(file_path)["id"].dropna().drop_duplicates()
        if ids.empty:
            continue
        df = pd.DataFrame(columns=["id", "region_start", "region_end"])
        df["id"] = ids.apply(lambda x: extract_uid_from_recordID(x)).copy()
        regions = ids.apply(lambda x: extract_region_from_recordID(x))
        df["region_start"] = regions.apply(lambda x: x[0])
        df["region_end"] = regions.apply(lambda x: x[1])
        dfs_ids_with_regions.append(df)
    
    if not dfs_ids_with_regions:
        return pd.DataFrame(columns=["id", "region_start", "region_end"])
    
    ids_with_regions = pd.concat(dfs_ids_with_regions, ignore_index=True)
    # remove duplicates and take the average of the region start and end for duplicate ids
    ids_with_regions = ids_with_regions.groupby("id").mean().astype({"region_start": "Int64", "region_end": "Int64"}).reset_index()
    
    return ids_with_regions
