import io
import tarfile
import time
import pandas as pd
from scripts.utils import *
import subprocess
import os
from pathlib import Path
import requests


def filter_foldseek_output(input_file: str, TM_SCORE_CUTOFF: float) -> "pd.DataFrame":
    """Reads the input file and filters for structures with TM-score >= TM_SCORE_CUTOFF."""
    df = pd.read_csv(input_file, sep="\t", header=None, names=["query", "target", "evalue", "qtmscore", "ttmscore"])
    filtered_df = df[df["qtmscore"] >= TM_SCORE_CUTOFF].copy()
    print(f"Remaining structures after filtering of foldseek output: {len(filtered_df)}/{len(df)}")
    return filtered_df[["query", "target", "qtmscore"]]


def filter_usalign_output(input_file: str, TM_SCORE_CUTOFF: float) -> "pd.DataFrame":
    """
    Parses alignment data output containing a Slurm job footer into a pandas DataFrame.
    """
    data_lines = []
    
    with open(input_file, 'r') as file:
        for line in file:
            # check for the job information footer
            if line.startswith('=='):
                break
            
            if not line.strip():
                continue
                
            # removing header
            if line.startswith('#'):
                line = line.lstrip('#')
                
            data_lines.append(line)
            
    # Combine the filtered lines into a single string
    raw_data = "".join(data_lines)
    
    # Read the string into a pandas DataFrame using a whitespace separator
    df = pd.read_csv(io.StringIO(raw_data), sep=r'\s+')
    
    df = df.rename(columns={"PDBchain1": "target", "PDBchain2": "query", "TM2": "qtmscore"})
    df = df[["query", "target", "qtmscore"]]
    
    filtered_df = df[df["qtmscore"] >= TM_SCORE_CUTOFF].copy()
    
    # clean up query and target columns to only hold uids
    filtered_df['target'] = filtered_df['target'].apply(lambda x: Path(x).stem)
    filtered_df['query'] = filtered_df['query'].apply(lambda x: Path(x).stem)
    
    print(f"Remaining structures after filtering usalign output: {len(filtered_df)}/{len(df)}")
    return filtered_df


def run_foldseek(query_structure: str, target_structures: str, root_dir: str) -> str:
    output_path = f"{root_dir}/data/filtering/foldseek/output_foldseek.csv"

    try:
        result = subprocess.run(
            [
                "foldseek",
                "easy-search",
                query_structure,
                target_structures,
                output_path,
                f"{root_dir}/data/filtering/foldseek/tmp",
                "--format-output",
                "query,target,evalue,qtmscore,ttmscore",
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        return output_path

    except subprocess.CalledProcessError as e:
        print("Foldseek failed!")
        print("Exit code:", e.returncode)
        print("STDOUT:")
        print(e.stdout)
        print("STDERR:")
        print(e.stderr)
        raise


def run_usalign(query_structure: str, target_structures: str, target_list: str, root_dir: str) -> str:
    out_path = f"{root_dir}/data/filtering/usalign/usalign_results.txt"

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
            print("USalign failed")
            print("Exit code:", e.returncode)
            print("STDERR:")
            print(e.stderr)
            raise

    return out_path
    
    
def create_target_list_for_usalign(target_structures_dir: str, root_dir: str) -> str:
    target_list_path = f"{root_dir}/data/filtering/usalign/target_list.txt"
    with open(target_list_path, "w") as f:
        for file in os.listdir(target_structures_dir):
            if file.endswith(".cif") or file.endswith(".pdb"):
                f.write(file + "\n")
    return target_list_path


def run_structural_similarity_filter(seed: list, number_tracker: dict, TM_SCORE_CUTOFF: float, root_dir: str) -> dict:
    query_structure_roi = f"{root_dir}/data/seed/{seed[0]}_roi.cif"
    target_structures_dir = f"{root_dir}/data/alphafold_structures"
    
    # run foldseek 
    foldseek_output = run_foldseek(query_structure_roi, target_structures_dir, root_dir)
    
    # run USalign
    target_list = create_target_list_for_usalign(target_structures_dir, root_dir)
    usalign_output = run_usalign(query_structure_roi, target_structures_dir, target_list, root_dir)
    
    # filter output
    filtered_ids_foldseek_df = filter_foldseek_output(foldseek_output, TM_SCORE_CUTOFF)
    filtered_ids_usalign_df = filter_usalign_output(usalign_output, TM_SCORE_CUTOFF)
    total_filtered_df = pd.merge(filtered_ids_usalign_df, filtered_ids_foldseek_df, on=["query", "target"], how="left").drop_duplicates().reset_index(drop=True)
    
    filtered_ids_foldseek_df.to_csv(f"{root_dir}/data/filtering/foldseek/filtered_foldseek_output_tm_{TM_SCORE_CUTOFF}.csv", index=False)
    filtered_ids_usalign_df.to_csv(f"{root_dir}/data/filtering/usalign/filtered_usalign_output_tm_{TM_SCORE_CUTOFF}.csv", index=False)
    total_filtered_df.to_csv(f"{root_dir}/data/filtering/structural_filtering_output_tm_{TM_SCORE_CUTOFF}.csv", index=False)
    
    # populate number tracker
    number_tracker["number_filtered_structures_foldseek"] = len(filtered_ids_foldseek_df)
    number_tracker["number_filtered_structures_usalign"] = len(filtered_ids_usalign_df)
    number_tracker["number_filtered_structures_total"] = len(total_filtered_df)
    
    return number_tracker


def run_foldseek_webserver(number_tracker: dict, path_to_query_structure: str, root_dir: str, TM_SCORE_CUTOFF: float):
    url = "https://search.foldseek.com/api/ticket"
    file = path_to_query_structure
    with open(file, "rb") as f:
        files = { "q": f }
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
            
            df_results = df_parsed_response[["query", "target", "eval", "dbStartPos", "dbEndPos"]]
            df_results.rename(columns={"target": "id",
                                       "eval": "qtmscore", # eval seems to be the TM-score
                                       "dbStartPos": "target_start",
                                       "dbEndPos": "target_end"}, inplace=True) 
            df_results = df_results[df_results["qtmscore"] >= TM_SCORE_CUTOFF]
            
            df_results.to_csv(f"{root_dir}/data/filtering/foldseek/foldseek_webserver_output.csv", index=False)
            number_tracker["number_filtered_structures_foldseek_webserver"] = len(df_results)
        return number_tracker
        

def wait_for_foldseek_webserver_job_completion(job_id: str) -> dict:
    while True:
        status = requests.get(
            f"https://search.foldseek.com/api/ticket/{job_id}"
        ).json()

        print(status)

        if status["status"] == "COMPLETE":
            return status

        if status["status"] == "ERROR":
            raise RuntimeError(status)

        time.sleep(5)
        
        
def parse_foldseek_webserver_response(data) -> pd.DataFrame:
    rows = []
    # for each db 
    for db_result in data.get("results", []):
        db_name = db_result.get("db")
        # for each alignment in db
        for aln in db_result.get("alignments", []):
            # for each hsp in alignment
            for hsp_dict in aln:
                row = {"db": db_name}
                # add all hsp_dict keys to row
                for key in hsp_dict.keys():
                    row[str(key)] = hsp_dict[str(key)]
                rows.append(row)

    return pd.DataFrame(rows)
                    