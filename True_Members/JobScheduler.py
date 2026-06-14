import os
import pandas as pd
from scripts.JobResultCollector import JobResultCollector
from scripts.utils import copy_file, read_yaml_config
from scripts.JobStatus import JobStatus
import yaml
import subprocess
import time



class JobScheduler:
    def __init__(self, config, global_data_dir,seed_list_path, used_seeds_path, all_members_path, results_dir, job_configs_dir, job_slurm_scripts_dir, job_config_template_path, evc_config_template_path, submit_job_slurm_template_path):
        self.start_config = config      
        self.MAX_ITERATIONS = self.start_config["global"]["iterations"]
        self.MAX_PARALLEL_JOBS = self.start_config["global"]["max_parallel_jobs"]
        self.CONVERGENCE_RUN = self.start_config["global"]["convergence_run"]
        self.ABSOLUTE_FLOOR = self.start_config["global"]["absolute_floor"]
        self.PATIENCE_WINDOW = self.start_config["global"]["patience_window"]
                
        self.seed_list_path = seed_list_path
        self.used_seeds_path = used_seeds_path
        self.all_members_path = all_members_path
        
        self.job_config_template_path = job_config_template_path
        self.evc_config_template_path = evc_config_template_path
        self.submit_job_slurm_template_path = submit_job_slurm_template_path

        self.global_data_dir = global_data_dir
        self.results_dir = results_dir
        self.job_configs_dir = job_configs_dir
        self.job_slurm_scripts_dir = job_slurm_scripts_dir
        os.makedirs(self.global_data_dir, exist_ok=True)
        os.makedirs(self.results_dir, exist_ok=True)
        os.makedirs(self.job_configs_dir, exist_ok=True)
        os.makedirs(self.job_slurm_scripts_dir, exist_ok=True)
        
        self.current_iteration = 1
        self.jobs: dict[int, dict] = {} # contains job_id as key and dict with keys "status" and "config_path" as value
        self.job_results_collector = JobResultCollector(self.seed_list_path, self.used_seeds_path, self.all_members_path)


    def start(self):
        while True:
            self.jobs = self.update_job_status()
            self.job_results_collector.get_results_from_finished_jobs(self.jobs)
            active_jobs_count = self.get_active_jobs_count()
            if self.current_iteration != 1:
                seed_list = pd.read_csv(f"{self.seed_list_path}")
                used_seeds = pd.read_csv(f"{self.used_seeds_path}")
                available_seeds = seed_list[~seed_list["id"].isin(used_seeds["id"])]
                if not self.CONVERGENCE_RUN and self.current_iteration > self.MAX_ITERATIONS:
                    print(f"Reached maximum number of iterations ({self.MAX_ITERATIONS}). Stopping pipeline.")
                    self.end_pipeline()
                    break
                elif available_seeds.empty:
                    print(f"No available seeds to run.")
                    if active_jobs_count > 0:
                        print(f"Waiting for running jobs to finish. Currently running {active_jobs_count} jobs.")
                        time.sleep(180) # waiting for 3 minutes before checking again
                        continue
                    else:
                        print(f"No running jobs. Stopping pipeline.")
                        self.end_pipeline()
                        break
                elif self.MAX_PARALLEL_JOBS <= active_jobs_count:
                    print(f"Currently running {active_jobs_count}/{self.MAX_PARALLEL_JOBS} jobs. Waiting for jobs to finish.")
                    time.sleep(180) # waiting for 3 minutes before checking again
                    continue
                else:
                    if self.CONVERGENCE_RUN:
                        converged, number_below_floor = self.compute_convergence(self.ABSOLUTE_FLOOR, self.PATIENCE_WINDOW)
                        if converged:
                            print(f"Number of True Members converged with {number_below_floor} consecutive jobs below absolute floor {self.ABSOLUTE_FLOOR}. Stopping pipeline.")
                            self.end_pipeline()
                            break
                        else:
                            print(f"Number of True Members did not converge with {number_below_floor} consecutive jobs below absolute floor {self.ABSOLUTE_FLOOR}. Continuing pipeline.")
                            print(f"Starting iteration {self.current_iteration} until convergence.")
                            self.start_job()
                            self.current_iteration += 1
                    else:
                        print(f"Starting iteration {self.current_iteration}/{self.MAX_ITERATIONS}.")
                        self.start_job()
                        self.current_iteration += 1
            else:
                print(f"Starting initial iteration {self.current_iteration} of {self.MAX_ITERATIONS}")
                job_id = self.start_job()
                while True:
                    print(f"Waiting for initial job to finish ...")
                    if self.check_if_job_finished(job_id):
                        print("First iteration finished.")
                        self.current_iteration += 1
                        break
                    else:                        
                        time.sleep(180) # waiting for 3 minutes before checking again
                        
    
    def update_job_status(self) -> dict[int, dict]:
        try:
            for job_id in self.jobs.keys():

                result = subprocess.run(
                    [
                        "sacct",
                        "-j", str(job_id),
                        "-X",
                        "--format=State",
                        "--noheader",
                        "--parsable2"
                    ],
                    capture_output=True,
                    text=True,
                    check=True
                )
                state = result.stdout.strip().split()[0] if result.stdout.strip() else None
                if not state:
                    continue  # job not yet registered in sacct DB  
                elif self.jobs[job_id]["status"] == JobStatus.HARVESTED:
                    continue  # already harvested
                elif state == "PENDING":
                    self.jobs[job_id]["status"] = JobStatus.PENDING
                elif state == "RUNNING":
                    self.jobs[job_id]["status"] = JobStatus.RUNNING
                elif state == "COMPLETED":
                    self.jobs[job_id]["status"] = JobStatus.COMPLETED
                elif state in ("FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL"):
                    self.jobs[job_id]["status"] = JobStatus.FAILED
                else:
                    print(f"Unknown job state '{state}' for job {job_id}. Keeping previous status {self.jobs[job_id]['status']}.")

        except subprocess.CalledProcessError as e:
            print(f"sacct error:\n{e.stderr}")
            
        return self.jobs        
        
        
    def check_if_job_finished(self, job_id: int) -> bool:
        self.jobs = self.update_job_status()
        finished_jobs =(JobStatus.COMPLETED, JobStatus.HARVESTED, JobStatus.FAILED)
        return self.jobs[job_id]["status"] in finished_jobs
        
    
    def get_active_jobs_count(self) -> int:
        active_states = [JobStatus.RUNNING, JobStatus.PENDING]
        return len([job_id for job_id in self.jobs.keys() if self.jobs[job_id]["status"] in active_states])
                
                
    def end_pipeline(self):
        running_jobs_count = self.get_active_jobs_count()
        if running_jobs_count > 0:
            while True:
                self.jobs = self.update_job_status()
                running_jobs_count = self.get_active_jobs_count()
                
                if running_jobs_count == 0:
                    print(f"All jobs finished. Ending pipeline.")
                    break
                else:
                    print(f"Waiting for {running_jobs_count} running jobs to finish before ending pipeline.")
                    time.sleep(180) # waiting for 3 minutes before checking again
        
        self.jobs = self.update_job_status()
        self.job_results_collector.get_results_from_finished_jobs(self.jobs)
        print(f"Saving final results to {self.results_dir}.")
        copy_file(self.seed_list_path, f"{self.results_dir}/final_True_Members_i_{self.current_iteration}.csv")
        copy_file(self.used_seeds_path, f"{self.results_dir}/used_seeds_iteration_i_{self.current_iteration}.csv")
        copy_file(self.all_members_path, f"{self.results_dir}/all_members_iteration_i_{self.current_iteration}.csv")


    def compute_convergence(self, absolute_floor: int, patience_window: int) -> tuple[bool, int]:
        used_seeds_df = pd.read_csv(self.used_seeds_path)
        finished_seeds_df = used_seeds_df.dropna(subset=["new_true_members"]).sort_values(by="iteration").copy()
        
        if len(used_seeds_df) < patience_window:
            # patience window not reached yet
            jobs_in_window = finished_seeds_df
            jobs_in_window["below_floor"] = jobs_in_window["new_true_members"] < absolute_floor
            number_of_jobs_below_floor = jobs_in_window["below_floor"].sum()
            return False, number_of_jobs_below_floor
        else:
            jobs_in_window = used_seeds_df.tail(patience_window)
            jobs_in_window["below_floor"] = jobs_in_window["new_true_members"] < absolute_floor
            number_of_jobs_below_floor = jobs_in_window["below_floor"].sum()
            converged = jobs_in_window["below_floor"].all()
            return bool(converged), number_of_jobs_below_floor


    def start_job(self) -> int:
        job_config_path, job_config = self.create_job_config()
        slurm_script_path = self.create_slurm_script(job_config_path, job_config)
        self.add_current_seed_to_used_seeds(job_config)
        # sbatch slurm_script_path in subprocess            
        try:
            result = subprocess.run(
                ["sbatch", slurm_script_path],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True)
            
            job_id = int(result.stdout.strip().split()[-1])
            self.jobs[job_id] = {"status": JobStatus.PENDING, "config_path": job_config_path}
            return job_id
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Failed to submit job {slurm_script_path} with error:\n{e.stderr}") from e
        
        
    def add_current_seed_to_used_seeds(self, job_config):
        seed_row = pd.DataFrame([{
            "id": job_config["run"]["seed_uid"],
            "region_start": job_config["run"]["seed_region_start"],
            "region_end": job_config["run"]["seed_region_end"],
            "iteration": self.current_iteration,
            "new_true_members": pd.NA
        }])
        used_seeds_df = pd.read_csv(self.used_seeds_path)
        used_seeds_df = pd.concat([used_seeds_df, seed_row], ignore_index=True)
        used_seeds_df.to_csv(self.used_seeds_path, index=False)
                
    
    def create_slurm_script(self, job_config_path, job_config):
        with open(self.submit_job_slurm_template_path, "r") as f:
            slurm_template = f.read()
            
        # replace placeholder with job_config_path
        seed_uid = job_config["run"]["seed_uid"]
        slurm_script = slurm_template.replace("<seed_uid>", seed_uid)
        slurm_script = slurm_script.replace("<iteration>", str(self.current_iteration))
        slurm_script = slurm_script.replace("<job_config_path>", job_config_path)
        slurm_script = slurm_script.replace("<ev_couplings_config_template_path>", self.evc_config_template_path)

        # save slurm script
        job_slurm_script_path = f"{self.job_slurm_scripts_dir}/submit_job_{seed_uid}_iteration_{self.current_iteration}.sh"
        with open(job_slurm_script_path, "w") as f:
            f.write(slurm_script)
            
        return job_slurm_script_path
    
    
    def create_job_config(self):
        job_config = read_yaml_config(self.job_config_template_path).copy()
        
        if self.current_iteration == 1:
            new_seed = self.start_config["run"]["seed_uid"]
            job_config["run"]["seed_uid"] = new_seed
            job_config["run"]["seed_region_start"] = self.start_config["run"]["seed_region_start"]
            job_config["run"]["seed_region_end"] = self.start_config["run"]["seed_region_end"]
            
        else:
            top_seed_row = self.get_highest_scoring_seed()
            new_seed = top_seed_row["id"].iloc[0]
            job_config["run"]["seed_uid"] = new_seed
            job_config["run"]["seed_region_start"] = top_seed_row["region_start"].iloc[0]
            job_config["run"]["seed_region_end"] = top_seed_row["region_end"].iloc[0]
            
        job_config["global"]["iteration"] = self.current_iteration
        job_config["global"]["root_dir"] = f"runs/true_members_run_{new_seed}_{time.strftime('%Y-%m-%d_%H-%M-%S')}"
        job_config["run"]["bit_scores"] = self.start_config["run"]["bit_scores"]
        job_config["run"]["tm_score_cutoff"] = self.start_config["run"]["tm_score_cutoff"]
        job_config["run"]["hmms_folder_path"] = self.start_config["run"]["hmms_folder_path"]
        job_config["run"]["CPU"] = self.start_config["run"]["CPU"]
        
        job_config_path = f"{self.job_configs_dir}/job_{new_seed}_iteration_{self.current_iteration}.yaml"
        with open(job_config_path, "w") as f:
            yaml.safe_dump(job_config, f, sort_keys=False)
            
        return job_config_path, job_config
                
                
    def get_highest_scoring_seed(self) -> pd.DataFrame:
        seeds_df = pd.read_csv(self.seed_list_path)
        used_seeds_df = pd.read_csv(self.used_seeds_path)
        available_seeds_df = seeds_df[~seeds_df["id"].isin(used_seeds_df["id"])]
        available_seeds_df = available_seeds_df.sort_values(by="qtmscore", ascending=False)
        top_seed_row = available_seeds_df.iloc[[0]]
        return top_seed_row
                