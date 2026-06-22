from True_Members.Job import Job
import sys
import os


def read_sys_args() -> tuple[str, str]:
        args = sys.argv
        if len(args) < 3 or len(args) > 3:
            print("Usage: python run_job.py <path_to_job_config_file> <path_to_evcouplings_config_template>")
            sys.exit(1)
        else:
            job_config_path = args[1]
            evcouplings_config_template = args[2]
            return job_config_path, evcouplings_config_template
        
        
if __name__ == "__main__":
    job_config_path, evcouplings_config_template = read_sys_args()
    job = Job(job_config_path, evcouplings_config_template)
    job.run()
    os._exit(0) # stop finished job from hanging due to multiprocessing pool not closing properly
    