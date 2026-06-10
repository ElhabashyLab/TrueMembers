from scripts.utils import *
import pandas as pd
import os
    
    
def get_all_structures_and_sequences(seed: list, bit_scores: list[float], number_tracker: dict, root_dir: str) -> dict:
    
    # Get annotation files
    annotation_files = get_annotation_files(root_dir, seed, bit_scores)
    number_tracker["annotation_files_jackhmmer"] = len(annotation_files)
    
    # Finding all unique IDs across annotation files
    unique_ids_with_regions = get_unique_ids_with_regions_from_annotation_files(annotation_files)
    number_tracker["unique_ids_jackhmmer"] = len(unique_ids_with_regions)
    
    # Get unique IDs from foldseek_webserver_output
    foldseek_webserver_output_path = f"{root_dir}/data/filtering/foldseek/foldseek_webserver_output.csv"
    foldseek_webserver_ids_df = pd.read_csv(foldseek_webserver_output_path)[["id"]]
    foldseek_webserver_ids_df["id"] = foldseek_webserver_ids_df["id"].apply(lambda x: extract_uid_from_recordID(x))
    
    # Combining all unique IDs
    combined_unique_ids = pd.concat([unique_ids_with_regions, foldseek_webserver_ids_df]).drop_duplicates()
    # remove duplicates and take the average of the region start and end for duplicate ids
    combined_unique_ids = (
        combined_unique_ids.groupby("id", as_index=False)
        .mean(numeric_only=True)
        .round()
        .astype({"id": str, "region_start": int, "region_end": int})
    )
    
    uids = list(combined_unique_ids["id"])
    number_tracker["unique_ids_combined"] = len(combined_unique_ids)  
    number_tracker["unique_uniprot_uids"] = len(uids)

    # Download structures from AlphaFold
    structures_output_dir = f"{root_dir}/data/alphafold_structures"
    successful_uids = run_download_tasks(uids, structures_output_dir)
    successful_uids_with_region = combined_unique_ids[combined_unique_ids["id"].isin(successful_uids)] 
    successful_uids_with_region.to_csv(f"{root_dir}/data/all_sequences_with_region.csv", index=False)
    number_tracker["successful_structure_downloads"] = len(successful_uids)
    
    # Download all records which have a successfully downloaded structure
    records = get_records_from_accessions(successful_uids_with_region["id"].tolist())
    number_tracker["downloaded_records"] = len(records)
    
    # Write records to fasta
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
    ids_with_regions = ids_with_regions.groupby("id").mean().reset_index().mean().astype({"id": str, "region_start": int, "region_end": int})
    
    return ids_with_regions
