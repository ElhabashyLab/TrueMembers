import logging
import os

import pandas as pd

from True_Members.scripts.types import JobStats
from True_Members.scripts.utils import (extract_uid_from_recordID, extract_region_from_recordID,
                           get_records_from_accessions, run_download_tasks, write_fasta_from_records)

logger = logging.getLogger(__name__)


def get_all_structures_and_sequences(seed: list, bit_scores: list[float], stats: JobStats, evcouplings_dir: str, foldseek_dir: str, shared_structures_dir: str, shared_sequences_dir: str, root_dir: str, local_fasta_path: str = "") -> tuple[str, pd.DataFrame]:

    unique_ids_with_regions, foldseek_webserver_ids_df = get_jackhmmer_and_foldseek_ids_with_regions(stats, seed, bit_scores, evcouplings_dir, foldseek_dir)
    combined_unique_ids_with_regions = combine_foldseek_and_jackhmmer_ids(unique_ids_with_regions, foldseek_webserver_ids_df)
    stats.ids_combined = len(combined_unique_ids_with_regions["id"])

    ids = combined_unique_ids_with_regions["id"].tolist()
    successful_ids = run_download_tasks(ids, shared_structures_dir)
    combined_unique_ids_with_regions = combined_unique_ids_with_regions[combined_unique_ids_with_regions["id"].isin(successful_ids)]
    stats.structure_downloads = len(successful_ids)

    records = get_records_from_accessions(combined_unique_ids_with_regions["id"].tolist(), cache_dir=shared_sequences_dir, local_fasta_path=local_fasta_path)
    logger.info("%d records downloaded from Uniprot.", len(records))
    stats.uniprot_records = len(records)

    record_ids = [extract_uid_from_recordID(record.id) for record in records]
    uids_with_region = combined_unique_ids_with_regions[combined_unique_ids_with_regions["id"].isin(record_ids)]

    uids_with_region_fasta_path = f"{root_dir}/uids_with_region.fasta"
    uids_with_region_path = f"{root_dir}/uids_with_region.csv"
    uids_with_region_full_fasta_path = write_fasta_from_records(records, uids_with_region_fasta_path)
    uids_with_region.to_csv(uids_with_region_path, index=False)

    return uids_with_region_full_fasta_path, uids_with_region


def get_jackhmmer_and_foldseek_ids_with_regions(stats: JobStats, seed: list, bit_scores: list[float], evcouplings_dir: str, foldseek_dir: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    annotation_files = get_annotation_files(evcouplings_dir, seed, bit_scores)
    stats.annotation_files_jackhmmer = len(annotation_files)

    unique_ids_with_regions = get_unique_ids_with_regions_from_annotation_files(annotation_files)
    stats.ids_jackhmmer = len(unique_ids_with_regions)

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
            logger.warning("Annotation file %s does not exist.", annotation_file)

    return annotation_files


def get_unique_ids_with_regions_from_annotation_files(annotation_files: list[str]) -> pd.DataFrame:
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
    ids_with_regions = (
        ids_with_regions
        .groupby("id")
        .agg(region_start=("region_start", "min"), region_end=("region_end", "max"))
        .astype({"region_start": "Int64", "region_end": "Int64"})
        .reset_index())

    return ids_with_regions
