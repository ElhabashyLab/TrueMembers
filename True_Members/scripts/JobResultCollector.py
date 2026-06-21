import logging
import os

import pandas as pd

from scripts.JobStatus import JobStatus
from scripts.utils import read_yaml_config

logger = logging.getLogger(__name__)


class JobResultCollector:

    def __init__(self, seed_list_path: str, used_seeds_path: str, all_members_path: str):
        self.seed_list_path = seed_list_path
        self.used_seeds_path = used_seeds_path
        self.all_members_path = all_members_path

    def get_results_from_finished_jobs(self, jobs: dict[int, dict]) -> dict[int, dict]:
        finished_jobs = [job_id for job_id in jobs.keys() if jobs[job_id]["status"] == JobStatus.COMPLETED]
        for job_id in finished_jobs:
            try:
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

                logger.info("Collected results from job %d (seed %s): %d new true members.",
                            job_id, seed_id_of_job, number_new_true_members)
                jobs[job_id]["status"] = JobStatus.HARVESTED
            except Exception as e:
                logger.error("Error while collecting results from job %d: %s", job_id, e)
        return jobs

    def merge_local_all_members_with_global_all_members(self, job_id: int, local_results_dir: str) -> pd.DataFrame:
        global_all_members_df = pd.read_csv(os.path.join(self.all_members_path))
        local_all_members_df = pd.read_csv(f"{local_results_dir}/All_Members.csv")
        dfs_combined = pd.concat([global_all_members_df, local_all_members_df], ignore_index=True)

        bool_cols = [c for c in dfs_combined.columns if c not in ["id", "qtmscore", "region_start", "region_end"]]
        agg_rules = {
            "qtmscore": "median",
            "region_start": "min",
            "region_end": "max"
        }
        for col in bool_cols:
            agg_rules[col] = "any"

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

        bool_cols = [c for c in dfs_combined.columns if c not in ["id", "qtmscore", "region_start", "region_end"]]

        agg_rules = {
            "qtmscore": "median",
            "region_start": "min",
            "region_end": "max"
        }
        for col in bool_cols:
            agg_rules[col] = "all"

        merged_seeds = dfs_combined.groupby("id", as_index=False).agg(agg_rules)
        sorted_merged_seeds = merged_seeds.sort_values("qtmscore", ascending=False)

        seeds_length_before_merge = len(global_seeds_df)
        seeds_length_after_merge = len(sorted_merged_seeds)
        number_new_true_members = seeds_length_after_merge - seeds_length_before_merge

        return sorted_merged_seeds, number_new_true_members