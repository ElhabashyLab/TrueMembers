from scripts.JobStatus import JobStatus
import pandas as pd
import os   
from scripts.utils import read_yaml_config


class JobResultCollector:
    
    def __init__(self, seed_list_path: str, used_seeds_path: str, all_members_path: str):
        self.seed_list_path = seed_list_path
        self.used_seeds_path = used_seeds_path
        self.all_members_path = all_members_path
        
    
    def get_results_from_finished_jobs(self, jobs: dict[int, dict]) -> dict[int, dict]:
        finished_jobs = [job_id for job_id in jobs.keys() if jobs[job_id]["status"] == JobStatus.COMPLETED]
        for job_id in finished_jobs:
            
            # job_config_path looks like this: f"{self.job_configs_dir}/job_{new_seed}_iteration_{self.current_iteration}.yaml"
            job_config_path = jobs[job_id]["config_path"]
            job_config = read_yaml_config(job_config_path)
            seed_id_of_job = job_config["run"]["seed_uid"]
            job_root_dir = job_config["global"]["root_dir"]
            local_results_dir = f"{job_root_dir}/results"
            
            sorted_merged_seeds, number_new_true_members = self.merge_local_true_members_with_global_seeds(job_id, local_results_dir)
            sorted_merged_seeds.to_csv(self.seed_list_path, index=False)
            
            updated_used_seeds_df = self.add_current_seed_info_to_used_seeds(seed_id_of_job, number_new_true_members)
            updated_used_seeds_df.to_csv(self.used_seeds_path, index=False)
            
            merged_all_members_global = self.merge_local_all_members_with_global_all_members(job_id, local_results_dir)
            merged_all_members_global.to_csv(self.all_members_path, index=False)
            
            jobs[job_id]["status"] = JobStatus.HARVESTED
        return jobs
        
        
    def merge_local_all_members_with_global_all_members(self, job_id: int, local_results_dir: str) -> pd.DataFrame:
        global_all_members_df = pd.read_csv(os.path.join(self.all_members_path))
        local_all_members_df = pd.read_csv(f"{local_results_dir}/All_Members.csv")
        dfs_combined = pd.concat([global_all_members_df, local_all_members_df], ignore_index=True)
        
        # identify hmm domain columns
        bool_cols = [c for c in dfs_combined.columns if c not in ["id", "qtmscore", "region_start", "region_end"]]
        # aggregation rules for key columns
        agg_rules = {
            "qtmscore": "median", # robust
            "region_start": "min", # use lowest start region found
            "region_end": "max" # use highest end region found
        }
        # add aggregation rules for hmm domain columns
        for col in bool_cols:
            agg_rules[col] = "any" # some entries might be true

        global_all_members_mreged_df = dfs_combined.groupby("id", as_index=False).agg(agg_rules)
        return global_all_members_mreged_df 
    
    
    def add_current_seed_info_to_used_seeds(self, seed_id: str, number_new_true_members: int) -> pd.DataFrame:
        used_seeds_df = pd.read_csv(os.path.join(self.used_seeds_path))      
        used_seeds_df.loc[used_seeds_df["id"] == seed_id, "new_true_members"] = number_new_true_members
        return used_seeds_df
        
    
    def merge_local_true_members_with_global_seeds(self, job_id: int, local_results_dir: str) -> tuple[pd.DataFrame, int]:  
        global_seeds_df = pd.read_csv(self.seed_list_path)
        local_true_members_df = pd.read_csv(f"{local_results_dir}/True_Members.csv")
        dfs_combined = pd.concat([global_seeds_df, local_true_members_df], ignore_index=True)
        
        # identify hmm domain columns (all entries are true)
        bool_cols = [c for c in dfs_combined.columns if c not in ["id", "qtmscore", "region_start", "region_end"]]

        # aggregation rules for key columns
        agg_rules = {
            "qtmscore": "median", # robust
            "region_start": "min", # use lowest start region found
            "region_end": "max" # use highest end region found
        }
        # aggregation rules for hmm domain columns
        for col in bool_cols:
            agg_rules[col] = "all" # all entries should be true

        merged_seeds = dfs_combined.groupby("id", as_index=False).agg(agg_rules)
        sorted_merged_seeds = merged_seeds.sort_values("qtmscore", ascending=False)
        
        # compute number of new true members found in this run
        seeds_length_before_merge = len(global_seeds_df)
        seeds_length_after_merge = len(sorted_merged_seeds)
        number_new_true_members = seeds_length_after_merge - seeds_length_before_merge
        
        return sorted_merged_seeds, number_new_true_members
    