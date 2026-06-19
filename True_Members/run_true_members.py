import time
import pandas as pd
from JobScheduler import JobScheduler
from scripts.utils import read_yaml_config
import sys
import os


def read_sys_args() -> str:
        args = sys.argv
        if len(args) < 2 or len(args) > 2:
            print("Usage: python run_true_members.py <path_to_config_file>")
            sys.exit(1)
        else:
            config_path = args[1]
            return config_path
        
        
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
        "JOB_CONFIGS_DIR": f"{run_dir}/job_configs",
        "JOB_SLURM_SCRIPTS_DIR": f"{run_dir}/slurm_scripts",
        "RESULTS_DIR": f"{run_dir}/final_true_members_results",
        
        "SEEDS_PATH": f"{run_dir}/data/seeds.csv",
        "USED_SEEDS_PATH": f"{run_dir}/data/used_seeds.csv",
        "ALL_MEMBERS_PATH": f"{run_dir}/data/All_Members_global.csv",
        
        "JOB_CONFIG_TEMPLATE": "template_files/job_config_template.yaml",
        "EVC_CONFIG_TEMPLATE": "template_files/evcouplings_config_template.yaml",
        "SUBMIT_JOB_SLURM_TEMPLATE": "template_files/submit_true_members_job_template.sh"
        }
    
    # create directories
    os.makedirs(paths["GLOBAL_RUN_DIR"], exist_ok=True)
    os.makedirs(paths["GLOBAL_DATA_DIR"], exist_ok=True)
    os.makedirs(paths["JOB_CONFIGS_DIR"], exist_ok=True)
    os.makedirs(paths["JOB_SLURM_SCRIPTS_DIR"], exist_ok=True)
    os.makedirs(paths["RESULTS_DIR"], exist_ok=True)
    
    # create empty csv files for seeds, used seeds and all members
    create_csv_files(paths["SEEDS_PATH"], paths["USED_SEEDS_PATH"], paths["ALL_MEMBERS_PATH"])
    
    return paths


def create_csv_files(seeds_path, used_seeds_path, all_members_path):
    """Creates empty seeds, used_seeds, and all_members CSV files with their base schemas.
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


# Before running check the config files and makes sure run_evcouplings is uncommented
if __name__ == "__main__":
    config_path = read_sys_args()
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
                                 paths["SEEDS_PATH"],
                                 paths["USED_SEEDS_PATH"],
                                 paths["ALL_MEMBERS_PATH"],
                                 paths["RESULTS_DIR"],
                                 paths["JOB_CONFIGS_DIR"],
                                 paths["JOB_SLURM_SCRIPTS_DIR"],
                                 paths["JOB_CONFIG_TEMPLATE"],
                                 paths["EVC_CONFIG_TEMPLATE"],
                                 paths["SUBMIT_JOB_SLURM_TEMPLATE"])
    job_scheduler.start()
    