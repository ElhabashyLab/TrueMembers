import sys
import os
import time
from pathlib import Path

import pandas as pd

from True_Members.JobScheduler import JobScheduler
from True_Members.scripts.utils import read_yaml_config

PACKAGE_DIR = Path(__file__).resolve().parent
TEMPLATE_DIR = PACKAGE_DIR / "template_files"


def main() -> None:
    if len(sys.argv) != 2:
        print(f"Usage: truemembers <path_to_config_file>", file=sys.stderr)
        sys.exit(1)

    config_path = sys.argv[1]
    config = read_yaml_config(config_path)

    seed = config["run"]["seed_uid"]
    is_convergence_run = config["global"]["convergence_run"]
    if not is_convergence_run:
        iterations = config["global"]["iterations"]
    else:
        iterations = None

    paths = set_up_dirs(seed, is_convergence_run, iterations)

    job_scheduler = JobScheduler(config,
                                 paths["GLOBAL_RUN_DIR"],
                                 paths["GLOBAL_DATA_DIR"],
                                 paths["SHARED_STRUCTURES_DIR"],
                                 paths["SHARED_SEQUENCES_DIR"],
                                 paths["SEEDS_PATH"],
                                 paths["USED_SEEDS_PATH"],
                                 paths["ALL_MEMBERS_PATH"],
                                 paths["SEQUENCE_IDENTITY_DIR"],
                                 paths["SEED_TO_MEMBERS_SEQUENCE_IDENTITY_PATH"],
                                 paths["RESULTS_DIR"],
                                 paths["JOB_CONFIGS_DIR"],
                                 paths["JOB_SLURM_SCRIPTS_DIR"],
                                 paths["JOB_CONFIG_TEMPLATE"],
                                 paths["EVC_CONFIG_TEMPLATE"],
                                 paths["SUBMIT_JOB_SLURM_TEMPLATE"])
    job_scheduler.start()
        
        
def set_up_dirs(seed, is_convergence_run, iterations):
    os.makedirs("runs", exist_ok=True)
    if is_convergence_run:
        run_dir = os.path.join("runs", f"convergence_run_seed_{seed}_{time.strftime('%Y-%m-%d_%H-%M-%S')}")
    else:
        run_dir = os.path.join("runs", f"iteration_{iterations}_run_seed_{seed}_{time.strftime('%Y-%m-%d_%H-%M-%S')}")
    os.makedirs(run_dir, exist_ok=True)
    
    paths = {
        "GLOBAL_RUN_DIR": run_dir,
        "GLOBAL_DATA_DIR": f"{run_dir}/data",
        "SHARED_STRUCTURES_DIR": f"{run_dir}/shared/structures",
        "SHARED_SEQUENCES_DIR": f"{run_dir}/shared/sequences",
        "JOB_CONFIGS_DIR": f"{run_dir}/job_configs",
        "JOB_SLURM_SCRIPTS_DIR": f"{run_dir}/slurm_scripts",
        "RESULTS_DIR": f"{run_dir}/final_true_members_results",
        "SEQUENCE_IDENTITY_DIR": f"{run_dir}/data/sequence_identity",
        
        "SEEDS_PATH": f"{run_dir}/data/seeds.csv",
        "USED_SEEDS_PATH": f"{run_dir}/data/used_seeds.csv",
        "ALL_MEMBERS_PATH": f"{run_dir}/data/All_Members_global.csv",
        "SEED_TO_MEMBERS_SEQUENCE_IDENTITY_PATH": f"{run_dir}/data/sequence_identity/seed_to_members_sequence_identity.csv",
        
        "JOB_CONFIG_TEMPLATE": str(TEMPLATE_DIR / "job_config_template.yaml"),
        "EVC_CONFIG_TEMPLATE": str(TEMPLATE_DIR / "evcouplings_config_template.yaml"),
        "SUBMIT_JOB_SLURM_TEMPLATE": str(TEMPLATE_DIR / "submit_true_members_job_template.sh")
        }
    
    # create directories
    for key in ["GLOBAL_RUN_DIR", "GLOBAL_DATA_DIR", "SHARED_STRUCTURES_DIR", "SHARED_SEQUENCES_DIR",
                "JOB_CONFIGS_DIR", "JOB_SLURM_SCRIPTS_DIR", "RESULTS_DIR", "SEQUENCE_IDENTITY_DIR"]:
        os.makedirs(paths[key], exist_ok=True)
    
    # create empty csv files for seeds, used seeds, all members and seed to members sequence identity
    create_csv_files(paths["SEEDS_PATH"], paths["USED_SEEDS_PATH"], paths["ALL_MEMBERS_PATH"], 
                     paths["SEED_TO_MEMBERS_SEQUENCE_IDENTITY_PATH"])
    
    return paths


def create_csv_files(seeds_path, used_seeds_path, all_members_path, seed_to_members_sequence_identity_path):
    """Creates empty seeds, used_seeds, all_members and seed_to_members_sequence_identity CSV files with their base schemas.
    HMM domain bool columns are not included here — they're added dynamically later based on the HMMs the user chose."""

    # seeds and all_members share the same base columns
    member_dtypes = {
        "id": "string",
        "qtmscore": "float64",
        "region_start": "Int64",
        "region_end": "Int64",
    }
    empty_members_df = pd.DataFrame({col: pd.Series(dtype=dt) for col, dt in member_dtypes.items()})
    empty_members_df.to_csv(seeds_path, index=False)
    empty_members_df.to_csv(all_members_path, index=False)

    used_seeds_dtypes = {
        "id": "string",
        "region_start": "Int64",
        "region_end": "Int64",
        "iteration": "Int64",
        "new_true_members": "Int64",
    }
    empty_used_seeds_df = pd.DataFrame({col: pd.Series(dtype=dt) for col, dt in used_seeds_dtypes.items()})
    empty_used_seeds_df.to_csv(used_seeds_path, index=False)
    
    # seed_to_members_sequence_identity is a TSV file with columns: seed_id, member_id, sequence_identity
    sequence_identity_dtypes = {
        "seed_id": "string",
        "member_id": "string",
        "sequence_identity": "float64",
    }
    empty_sequence_identity_df = pd.DataFrame({col: pd.Series(dtype=dt) for col, dt in sequence_identity_dtypes.items()})
    empty_sequence_identity_df.to_csv(seed_to_members_sequence_identity_path, index=False, sep="\t")
