import sys
import time
from pathlib import Path

from src.JobScheduler import JobScheduler
from src.scripts.utils import read_yaml_config
from src.scripts.pipeline_db import init_db

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
        run_dir = runs_dir / f"iterations_run_{max_iterations}_target_{target}_{timestamp}"

    data_dir = run_dir / "data"
    shared_structures_dir = run_dir / "shared" / "structures"
    shared_sequences_dir = run_dir / "shared" / "sequences"
    job_configs_dir = run_dir / "job_configs"
    job_slurm_scripts_dir = run_dir / "job_slurm_scripts"
    results_dir = run_dir / "results"
    db_path = data_dir / "pipeline.db"

    evc_config_template = TEMPLATE_DIR / "evcouplings_config_template.yaml"
    submit_job_template = TEMPLATE_DIR / "submit_true_members_job_template.sh"

    for d in (run_dir, data_dir, shared_structures_dir, shared_sequences_dir,
              job_configs_dir, job_slurm_scripts_dir, results_dir):
        d.mkdir(parents=True, exist_ok=True)

    init_db(db_path)

    return {
        "GLOBAL_RUN_DIR": run_dir,
        "GLOBAL_DATA_DIR": data_dir,
        "SHARED_STRUCTURES_DIR": shared_structures_dir,
        "SHARED_SEQUENCES_DIR": shared_sequences_dir,
        "JOB_CONFIGS_DIR": job_configs_dir,
        "JOB_SLURM_SCRIPTS_DIR": job_slurm_scripts_dir,
        "RESULTS_DIR": results_dir,
        "EVC_CONFIG_TEMPLATE": evc_config_template,
        "SUBMIT_JOB_SLURM_TEMPLATE": submit_job_template,
        "DB_PATH": db_path,
    }