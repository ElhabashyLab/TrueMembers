import io
import pandas as pd
from scripts.utils import *
import subprocess
import os
from pathlib import Path


def filter_foldseek_output(input_file: str) -> "pd.DataFrame":
    """Reads the input file and filters for structures with TM-score >= 0.6."""
    df = pd.read_csv(input_file, sep="\t", header=None, names=["query", "target", "evalue", "qtmscore", "ttmscore"])
    filtered_df = df[df["qtmscore"] >= 0.6].copy()
    print(f"Remaining structures after filtering of foldseek output: {len(filtered_df)}/{len(df)}")
    return filtered_df[["query", "target", "qtmscore"]]


def filter_usalign_output(input_file: str) -> "pd.DataFrame":
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
    
    filtered_df = df[df["qtmscore"] >= 0.6].copy()
    
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


def run_structural_similarity_filter(seed: list, number_tracker: dict, root_dir: str) -> dict:
    # get correct roi of seed
    seed_uid = seed[0]
    seed_r_start = seed[1]
    seed_r_end = seed[2]
    run_download_tasks([seed_uid], output_dir=f"{root_dir}/data/seed/", max_workers=1)
    query_structure_full = f"{root_dir}/data/seed/{seed_uid}.cif"
    query_structure_roi = extract_region_of_protein(query_structure_full, (seed_r_start, seed_r_end))
    
    target_structures_dir = f"{root_dir}/data/alphafold_structures"
    
    # run foldseek 
    foldseek_output = run_foldseek(query_structure_roi, target_structures_dir, root_dir)
    
    # run USalign
    target_list = create_target_list_for_usalign(target_structures_dir, root_dir)
    usalign_output = run_usalign(query_structure_roi, target_structures_dir, target_list, root_dir)
    
    # filter output
    filtered_ids_foldseek_df = filter_foldseek_output(foldseek_output)
    filtered_ids_usalign_df = filter_usalign_output(usalign_output)
    total_filtered_df = pd.concat([filtered_ids_foldseek_df, filtered_ids_usalign_df]).drop_duplicates(subset=["query", "target"]).reset_index(drop=True)
    filtered_ids_foldseek_df.to_csv(f"{root_dir}/data/filtering/foldseek/filtered_foldseek_output_tm_0.6.csv", index=False)
    filtered_ids_usalign_df.to_csv(f"{root_dir}/data/filtering/usalign/filtered_usalign_output_tm_0.6.csv", index=False)
    total_filtered_df.to_csv(f"{root_dir}/data/filtering/structural_filtering_output_tm_0.6.csv", index=False)
    
    # populate number tracker
    number_tracker["number_filtered_structures_foldseek"] = len(filtered_ids_foldseek_df)
    number_tracker["number_filtered_structures_usalign"] = len(filtered_ids_usalign_df)
    number_tracker["number_filtered_structures_total"] = len(total_filtered_df)
    
    return number_tracker
    