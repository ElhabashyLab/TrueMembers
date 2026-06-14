from JobScheduler import JobScheduler
from scripts.utils import read_yaml_config
import sys


def read_sys_args() -> str:
        args = sys.argv
        if len(args) < 2 or len(args) > 2:
            print("Usage: python run_true_members.py <path_to_config_file>")
            sys.exit(1)
        else:
            config_path = args[1]
            return config_path


# Before running check the config files and makes sure run_evcouplings is uncommented
if __name__ == "__main__":
    config_path = read_sys_args()
    config = read_yaml_config(config_path)
    
    GLOBAL_DATA_DIR = "data"
    SEEDS_PATH = "data/seeds.csv"
    USED_SEEDS_PATH = "data/used_seeds.csv"
    ALL_MEMBERS_PATH = "data/All_Members_global.csv"
    JOB_CONFIGS_DIR = "job_configs"
    JOB_SLURM_SCRIPTS_DIR = "slurm_scripts"
    RESULTS_DIR = "final_true_members_results"
    
    JOB_CONFIG_TEMPLATE = "template_files/true_members_config_template.yaml"
    EVC_CONFIG_TEMPLATE = "template_files/evcouplings_config_template.yaml"
    SUBMIT_JOB_SLURM_TEMPLATE = "template_files/submit_true_members_job_template.sh"
    
    job_scheduler = JobScheduler(config,
                                 GLOBAL_DATA_DIR,
                                 SEEDS_PATH,
                                 USED_SEEDS_PATH,
                                 ALL_MEMBERS_PATH,
                                 RESULTS_DIR,
                                 JOB_CONFIGS_DIR,
                                 JOB_SLURM_SCRIPTS_DIR,
                                 JOB_CONFIG_TEMPLATE,
                                 EVC_CONFIG_TEMPLATE,
                                 SUBMIT_JOB_SLURM_TEMPLATE)
    job_scheduler.start()
    