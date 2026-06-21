import io
import logging
import os
import time
from pathlib import Path

import pandas as pd
import requests
import subprocess

from scripts.types import JobStats
from scripts.utils import extract_uid_from_recordID

logger = logging.getLogger(__name__)


def filter_foldseek_output(input_file: str, TM_SCORE_CUTOFF: float) -> "pd.DataFrame":
    df = pd.read_csv(input_file, sep="\t", header=None, names=["query", "target", "evalue", "qtmscore", "ttmscore"])
    filtered_df = df[df["qtmscore"] >= TM_SCORE_CUTOFF].copy()
    logger.info("Remaining structures after filtering of foldseek output: %d/%d", len(filtered_df), len(df))
    return filtered_df[["query", "target", "qtmscore"]]


def filter_usalign_output(input_file: str, TM_SCORE_CUTOFF: float) -> "pd.DataFrame":
    data_lines = []

    with open(input_file, 'r') as file:
        for line in file:
            if line.startswith('=='):
                break

            if not line.strip():
                continue

            if line.startswith('#'):
                line = line.lstrip('#')

            data_lines.append(line)

    raw_data = "".join(data_lines)

    df = pd.read_csv(io.StringIO(raw_data), sep=r'\s+')

    df = df.rename(columns={"PDBchain1": "target", "PDBchain2": "query", "TM2": "qtmscore"})
    df = df[["query", "target", "qtmscore"]]

    filtered_df = df[df["qtmscore"] >= TM_SCORE_CUTOFF].copy()

    filtered_df['target'] = filtered_df['target'].apply(lambda x: Path(x).stem)
    filtered_df['query'] = filtered_df['query'].apply(lambda x: Path(x).stem)

    logger.info("Remaining structures after filtering usalign output: %d/%d", len(filtered_df), len(df))
    return filtered_df


def run_foldseek(query_structure: str, target_structures: str, foldseek_dir: str) -> str:
    output_path = f"{foldseek_dir}/output_foldseek.csv"

    try:
        subprocess.run(
            [
                "foldseek",
                "easy-search",
                query_structure,
                target_structures,
                output_path,
                f"{foldseek_dir}/tmp",
                "--format-output",
                "query,target,evalue,qtmscore,ttmscore",
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        return output_path

    except subprocess.CalledProcessError as e:
        logger.error("Foldseek failed! Exit code: %s\nSTDOUT:\n%s\nSTDERR:\n%s",
                     e.returncode, e.stdout, e.stderr)
        raise


def run_usalign(query_structure: str, target_structures: str, target_list: str, usalign_dir: str) -> str:
    out_path = f"{usalign_dir}/usalign_results.txt"

    with open(out_path, "w") as f:
        try:
            subprocess.run(
                [
                    "USalign",
                    "-dir1", target_structures,
                    target_list,
                    query_structure,
                    "-outfmt", "2"
                ],
                stdout=f,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
        except subprocess.CalledProcessError as e:
            logger.error("USalign failed! Exit code: %s\nSTDERR:\n%s", e.returncode, e.stderr)
            raise

    return out_path


def create_target_list_for_usalign(target_structures_dir: str, usalign_dir: str) -> str:
    target_list_path = f"{usalign_dir}/target_list.txt"
    with open(target_list_path, "w") as f:
        for file in os.listdir(target_structures_dir):
            if file.endswith(".cif") or file.endswith(".pdb"):
                f.write(file + "\n")
    return target_list_path


def run_structural_similarity_filter(query_structure_roi: str, stats: JobStats, TM_SCORE_CUTOFF: float, structures_dir: str, foldseek_dir: str, usalign_dir: str) -> pd.DataFrame:

    foldseek_output = run_foldseek(query_structure_roi, structures_dir, foldseek_dir)

    target_list = create_target_list_for_usalign(structures_dir, usalign_dir)
    usalign_output = run_usalign(query_structure_roi, structures_dir, target_list, usalign_dir)

    filtered_ids_foldseek_df = filter_foldseek_output(foldseek_output, TM_SCORE_CUTOFF)
    filtered_ids_usalign_df = filter_usalign_output(usalign_output, TM_SCORE_CUTOFF)
    filtered_ids_foldseek_df.to_csv(f"{foldseek_dir}/filtered_foldseek_output_tm_{TM_SCORE_CUTOFF}.csv", index=False)
    filtered_ids_usalign_df.to_csv(f"{usalign_dir}/filtered_usalign_output_tm_{TM_SCORE_CUTOFF}.csv", index=False)

    merged = pd.merge(filtered_ids_usalign_df, filtered_ids_foldseek_df, on=["query", "target"], how="outer", suffixes=("_usalign", "_foldseek")).drop_duplicates().reset_index(drop=True)
    merged["qtmscore"] = merged["qtmscore_usalign"].combine_first(merged["qtmscore_foldseek"])
    total_filtered_df = merged.drop(columns=["query", "qtmscore_usalign", "qtmscore_foldseek"])
    total_filtered_df = total_filtered_df.rename(columns={"target": "id"})

    filtering_dir = os.path.dirname(usalign_dir)
    total_filtered_df.to_csv(f"{filtering_dir}/seed_containing_sequences_tm_{TM_SCORE_CUTOFF}.csv", index=False)
    stats.filtered_structures_foldseek = len(filtered_ids_foldseek_df)
    stats.filtered_structures_usalign = len(filtered_ids_usalign_df)
    stats.filtered_structures_total = len(total_filtered_df)

    return total_filtered_df


def run_foldseek_webserver(stats: JobStats, path_to_query_structure: str, root_dir: str, TM_SCORE_CUTOFF: float, out_dir: str):
    url = "https://search.foldseek.com/api/ticket"
    file = path_to_query_structure
    with open(file, "rb") as f:
        files = {"q": f}
        data = [
            ("mode", "tmalign"),
            ("database[]", "BFVD"),
            ("database[]", "afdb50"),
            ("database[]", "afdb-swissprot"),
            ("database[]", "afdb-proteome"),
            ("database[]", "bfmd"),
            ("database[]", "cath50"),
            ("database[]", "mgnify_esm30"),
            ("database[]", "pdb100"),
            ("database[]", "gmgcl_id"),
        ]
        response = requests.post(url, files=files, data=data).json()
        job_id = response["id"]

        df_results = pd.DataFrame()
        if wait_for_foldseek_webserver_job_completion(job_id):
            results = requests.get(
                f"https://search.foldseek.com/api/result/{job_id}/0").json()
            df_parsed_response = parse_foldseek_webserver_response(results)

            df_results = df_parsed_response[["target", "eval", "dbStartPos", "dbEndPos"]]
            df_results.rename(columns={"target": "id",
                                       "eval": "qtmscore",
                                       "dbStartPos": "region_start",
                                       "dbEndPos": "region_end"}, inplace=True)
            df_results = df_results[df_results["qtmscore"] >= TM_SCORE_CUTOFF]
            df_results["id"] = df_results["id"].apply(lambda x: extract_uid_from_recordID(x))

            df_results.to_csv(f"{out_dir}/foldseek_webserver_output.csv", index=False)
            stats.ids_foldseek_webserver = len(df_results)


def wait_for_foldseek_webserver_job_completion(job_id: str, max_wait_seconds: int = 3600) -> dict:
    start = time.time()
    while True:
        if time.time() - start > max_wait_seconds:
            raise TimeoutError(f"Foldseek webserver job {job_id} did not complete within {max_wait_seconds}s.")
        status = requests.get(
            f"https://search.foldseek.com/api/ticket/{job_id}", timeout=(10, 60)
        ).json()
        logger.info("Foldseek webserver status: %s", status)
        if status["status"] == "COMPLETE":
            return status
        if status["status"] == "ERROR":
            raise RuntimeError(status)
        time.sleep(5)


def parse_foldseek_webserver_response(data) -> pd.DataFrame:
    rows = []
    for db_result in data.get("results", []):
        db_name = db_result.get("db")
        for aln in db_result.get("alignments", []):
            for hsp_dict in aln:
                row = {"db": db_name}
                for key in hsp_dict.keys():
                    row[str(key)] = hsp_dict[str(key)]
                rows.append(row)

    return pd.DataFrame(rows)