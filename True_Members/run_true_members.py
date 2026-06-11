from scripts.utils import *
from scripts.evcouplings import *
from scripts.get_structures_and_sequences import *
from scripts.structural_filter import *
from scripts.hmm_domain_filter import *
from Job import Job

CONFIG_TEMPLATE = "evcouplings_config_template.yaml"
        
def read_sys_args() -> str:
        args = sys.argv
        if len(args) < 2 or len(args) > 2:
            print("Usage: python run_true_members.py <path_to_config_file>")
            sys.exit(1)
        else:
            config_path = args[1]
            return config_path

if __name__ == "__main__":
    config_path = read_sys_args()
    job = Job(config_path)
    job.run()
    