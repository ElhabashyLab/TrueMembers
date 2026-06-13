import os
import pandas as pd
from scripts.utils import copy_file, read_yaml_config
import yaml
import subprocess
import time


class JobScheduler:
    def __init__(self, config, global_data_dir,seed_list_path, used_seeds_path, all_members_path, results_dir, job_configs_dir, job_slurm_scripts_dir, job_config_template_path, evc_config_template_path, submit_job_slurm_template_path):
        self.start_config = config      
        self.ITERATIONS = self.start_config["global"]["iterations"]
        self.CONVERGENCE_RUN = self.start_config["global"]["convergence_run"]
        self.ABSOLUTE_FLOOR = self.start_config["global"]["absolute_floor"]
        self.PATIENCE_WINDOW = self.start_config["global"]["patience_window"]
        self.global_data_dir = global_data_dir
        self.seed_list_path = seed_list_path
        self.used_seeds_path = used_seeds_path
        self.all_members_path = all_members_path
        self.results_dir = results_dir
        self.job_configs_dir = job_configs_dir
        self.job_slurm_scripts_dir = job_slurm_scripts_dir
        self.job_config_template_path = job_config_template_path
        self.evc_config_template_path = evc_config_template_path
        self.submit_job_slurm_template_path = submit_job_slurm_template_path
        os.makedirs(self.results_dir, exist_ok=True)
        os.makedirs(self.job_configs_dir, exist_ok=True)
        os.makedirs(self.job_slurm_scripts_dir, exist_ok=True)
        self.current_iteration = 1
    
    def start(self):
        while True:
            if not self.current_iteration is 1:
                # slow down endless loop check used_seeds TODO
                seed_list = pd.read_csv(f"{self.seed_list_path}")
                if seed_list.empty:
                    print("Seed list is empty. Stopping pipeline.")
                    self.end_pipeline()
                    break
                if not self.CONVERGENCE_RUN and self.current_iteration > self.ITERATIONS:
                    print(f"Reached maximum number of iterations ({self.ITERATIONS}). Stopping pipeline.")
                    self.end_pipeline()
                    break
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
                        print(f"Starting iteration {self.current_iteration}/{self.ITERATIONS}.")
                        self.start_job()
                        self.current_iteration += 1
            else:
                print(f"Starting initial iteration {self.current_iteration} of {self.ITERATIONS}")
                self.start_job()
                used_seeds_df = pd.read_csv(self.used_seeds_path)
                while True:
                    if len(used_seeds_df) < 1:
                        print("Waiting for first iteration to finish ...")
                        time.sleep(1800) # waiting for 30 minutes before checking again
                    else:
                        print("First iteration finished. Starting next iteration.")
                        self.current_iteration += 1
                        break
                
                
    def end_pipeline(self):
        print(f"Pipeline finished. Saving final results to {self.results_dir}.")
        copy_file(self.seed_list_path, f"{self.results_dir}/final_True_Members_i_{self.current_iteration}.csv")
        copy_file(self.used_seeds_path, f"{self.results_dir}/used_seeds_iteration_i_{self.current_iteration}.csv")
        copy_file(self.all_members_path, f"{self.results_dir}/all_members_iteration_i_{self.current_iteration}.csv")


    def compute_convergence(self, absolute_floor: int, patience_window: int) -> tuple[bool, int]:
        used_seeds_df = pd.read_csv(self.used_seeds_path)
        used_seeds_df = used_seeds_df.sort_values(by="iteration")
        
        if len(used_seeds_df) < patience_window:
            # patience window not reached yet
            jobs_in_window = used_seeds_df
        else:
            jobs_in_window = used_seeds_df.tail(patience_window)
            
        converged = jobs_in_window["new_true_members"] < absolute_floor
        number_of_jobs_below_floor = converged.sum()
        converged = converged.all()
        
        return bool(converged), number_of_jobs_below_floor


    def start_job(self):
        job_config_path, job_config = self.create_job_config()
        slurm_script_path = self.create_slurm_script(job_config_path, job_config)
        # sbatch slurm_script_path in subprocess            
        try:
            result = subprocess.run(
                ["sbatch", slurm_script_path],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Failed to submit job {slurm_script_path} with error:\n{e.stderr}") from e
                
    
    def create_slurm_script(self, job_config_path, job_config):
        with open(self.submit_job_slurm_template_path, "r") as f:
            slurm_template = f.read()
            
        # replace placeholder with job_config_path
        seed_uid = job_config["run"]["seed_uid"]
        slurm_script = slurm_template.replace("<seed_uid>", seed_uid)
        slurm_script = slurm_script.replace("<iteration>", str(self.current_iteration))
        slurm_script = slurm_script.replace("<job_config_path>", job_config_path)
        slurm_script = slurm_script.replace("<ev_couplings_config_template_path>", self.evc_config_template_path)
        slurm_script = slurm_script.replace("<global_data_dir>", self.global_data_dir)

        # save slurm script
        job_slurm_script_path = f"{self.job_slurm_scripts_dir}/submit_job_{seed_uid}_iteration_{self.current_iteration}.sh"
        with open(job_slurm_script_path, "w") as f:
            f.write(slurm_script)
            
        return job_slurm_script_path
    
    
    def create_job_config(self):
        job_config = read_yaml_config(self.job_config_template_path).copy()
        job_config["global"]["iteration"] = self.current_iteration
        job_config["run"]["bit_scores"] = self.start_config["run"]["bit_scores"]
        job_config["run"]["tm_score_cutoff"] = self.start_config["run"]["tm_score_cutoff"]
        job_config["run"]["hmms_folder_path"] = self.start_config["run"]["hmms_folder_path"]
        job_config["run"]["CPU"] = self.start_config["run"]["CPU"]
        
        if self.current_iteration == 1:
            new_seed = self.start_config["run"]["seed_uid"]
            job_config["run"]["seed_uid"] = new_seed
            job_config["run"]["seed_region_start"] = self.start_config["run"]["seed_region_start"]
            job_config["run"]["seed_region_end"] = self.start_config["run"]["seed_region_end"]
            
        else:
            top_seed_row = self.get_highest_scoring_seed()
            new_seed = top_seed_row["id"]
            job_config["run"]["seed_uid"] = new_seed
            job_config["run"]["seed_region_start"] = top_seed_row["region_start"]
            job_config["run"]["seed_region_end"] = top_seed_row["region_end"]
        
        job_config_path = f"{self.job_configs_dir}/job_{new_seed}_iteration_{self.current_iteration}.yaml"
        with open(job_config_path, "w") as f:
            yaml.safe_dump(job_config, f, sort_keys=False)
            
        return job_config_path, job_config
                
                
    def get_highest_scoring_seed(self):
        seeds_df = pd.read_csv(self.seed_list_path)
        used_seeds_df = pd.read_csv(self.used_seeds_path)
        # filter out used seeds from seeds_df
        seeds_df = seeds_df[~seeds_df["id"].isin(used_seeds_df["id"])]
        seeds_df = seeds_df.sort_values(by="qtmscore", ascending=False)
        top_seed_row = seeds_df.iloc[0]
        return top_seed_row
                