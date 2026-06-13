from Job import Job
import sys


def read_sys_args() -> tuple[str, str, str]:
        args = sys.argv
        if len(args) < 3 or len(args) > 4:
            print("Usage: python run_job.py <path_to_job_config_file> <path_to_evcouplings_config_template> <path_to_global_data_dir>")
            sys.exit(1)
        else:
            job_config_path = args[1]
            evcouplings_config_template = args[2]
            global_data_dir = args[3]
            return job_config_path, evcouplings_config_template, global_data_dir
        
        
if __name__ == "__main__":
    job_config_path, evcouplings_config_template, global_data_dir = read_sys_args()
    job = Job(job_config_path, evcouplings_config_template, global_data_dir)
    job.run()
    