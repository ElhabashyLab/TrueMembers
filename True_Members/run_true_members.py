from JobScheluder import JobScheduler
import yaml
import sys

CONFIG_TEMPLATE = "evcouplings_config_template.yaml"
        
def read_sys_args() -> str:
        args = sys.argv
        if len(args) < 2 or len(args) > 2:
            print("Usage: python run_true_members.py <path_to_config_file>")
            sys.exit(1)
        else:
            config_path = args[1]
            return config_path
        
        
def read_true_members_config(config_path):
        with open(config_path, 'r') as file:
            config = yaml.safe_load(file)
        return config


if __name__ == "__main__":
    config_path = read_sys_args()
    config = read_true_members_config(config_path)
    
    seeds_path = "seed_list.csv"
    used_seeds_path = "used_seeds.csv"
    all_members_path = "All_Members_global.csv"
    
    starting_iteration = 1
    
    job_scheduler = JobScheduler(config, starting_iteration, seeds_path, used_seeds_path, all_members_path)
    job_scheduler.start()
    