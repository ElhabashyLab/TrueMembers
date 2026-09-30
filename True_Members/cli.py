import sys
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

    target = config["global"]["target_uid"]
    is_convergence_run = config["global"]["convergence_mode"]
    if not is_convergence_run:
        max_iterations = config["global"]["max_iterations"]
    else:
        max_iterations = None

    paths = set_up_dirs(target, is_convergence_run, max_iterations)

    job_scheduler = JobScheduler(config, paths)
    job_scheduler.start()
        
        
def set_up_dirs(target, is_convergence_run, max_iterations):
    runs_dir = Path("runs")
    runs_dir.mkdir(exist_ok=True)

    timestamp = time.strftime("%Y-%m-%d_%H-%M-%S")
    if is_convergence_run:
        run_dir = runs_dir / f"convergence_run_target_{target}_{timestamp}"
    else:
        run_dir = runs_dir / f"iteration_{max_iterations}_run_target_{target}_{timestamp}"

    data_dir = run_dir / "data"
    seqid_dir = data_dir / "sequence_identity"
    shared_structures_dir = run_dir / "shared" / "structures"
    shared_sequences_dir = run_dir / "shared" / "sequences"
    job_configs_dir = run_dir / "job_configs"
    job_slurm_scripts_dir = run_dir / "job_slurm_scripts"
    results_dir = run_dir / "results"

    members_path = data_dir / "members.csv"
    used_seeds_path = data_dir / "used_seeds.csv"
    all_members_path = data_dir / "all_members.csv"
    seqid_path = seqid_dir / "seed_to_members_sequence_identity.csv"

    job_config_template = TEMPLATE_DIR / "job_config_template.yaml"
    evc_config_template = TEMPLATE_DIR / "evcouplings_config_template.yaml"
    submit_job_template = TEMPLATE_DIR / "submit_true_members_job_template.sh"

    for d in (run_dir, data_dir, seqid_dir, shared_structures_dir,
              shared_sequences_dir, job_configs_dir, job_slurm_scripts_dir, results_dir):
        d.mkdir(parents=True, exist_ok=True)

    create_csv_files(members_path, used_seeds_path, all_members_path, seqid_path)

    return {
        "GLOBAL_RUN_DIR": run_dir,
        "GLOBAL_DATA_DIR": data_dir,
        "SHARED_STRUCTURES_DIR": shared_structures_dir,
        "SHARED_SEQUENCES_DIR": shared_sequences_dir,
        "JOB_CONFIGS_DIR": job_configs_dir,
        "JOB_SLURM_SCRIPTS_DIR": job_slurm_scripts_dir,
        "RESULTS_DIR": results_dir,
        "SEQUENCE_IDENTITY_DIR": seqid_dir,
        "MEMBERS_PATH": members_path,
        "USED_SEEDS_PATH": used_seeds_path,
        "ALL_MEMBERS_PATH": all_members_path,
        "SEED_TO_MEMBERS_SEQUENCE_IDENTITY_PATH": seqid_path,
        "JOB_CONFIG_TEMPLATE": job_config_template,
        "EVC_CONFIG_TEMPLATE": evc_config_template,
        "SUBMIT_JOB_SLURM_TEMPLATE": submit_job_template,
    }


def create_csv_files(members_path, used_seeds_path, all_members_path, seed_to_members_sequence_identity_path):
    """Creates empty members, used_seeds, all_members and seed_to_members_sequence_identity CSV files with their base schemas.
    HMM domain bool columns are not included here — they're added dynamically later based on the HMMs the user chose."""

    # members and all_members share the same base columns
    member_dtypes = {
        "id": "string",
        "qtmscore": "float64",
        "region_start": "Int64",
        "region_end": "Int64",
    }
    empty_members_df = pd.DataFrame({col: pd.Series(dtype=dt) for col, dt in member_dtypes.items()})
    empty_members_df.to_csv(members_path, index=False)
    empty_members_df.to_csv(all_members_path, index=False)

    used_seeds_dtypes = {
        "id": "string",
        "region_start": "Int64",
        "region_end": "Int64",
        "iteration": "Int64",
        "finishing_rank": "Int64",
        "new_members": "Int64",
    }
    empty_used_seeds_df = pd.DataFrame({col: pd.Series(dtype=dt) for col, dt in used_seeds_dtypes.items()})
    empty_used_seeds_df.to_csv(used_seeds_path, index=False)
    
    sequence_identity_dtypes = {
        "seed_id": "string",
        "hit_id": "string",
        "pident": "float64",
        "detected": "bool"
    }
    empty_sequence_identity_df = pd.DataFrame({col: pd.Series(dtype=dt) for col, dt in sequence_identity_dtypes.items()})
    empty_sequence_identity_df.to_csv(seed_to_members_sequence_identity_path, index=False)
