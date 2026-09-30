from itertools import chain
import logging
import os

import pandas as pd

from Bio import SeqIO

from src.scripts.SequenceIdentityCalculator import SequenceIdentityCalculator
from src.scripts.JobStatus import JobStatus
from src.scripts.utils import read_yaml_config, get_records_from_accessions

logger = logging.getLogger(__name__)


class JobResultCollector:

    def __init__(self, 
                 seed_list_path: str, 
                 used_seeds_path: str, 
                 all_members_path: str, 
                 sequence_identity_dir: str, 
                 seed_to_hits_sequence_identity_path: str,
                 shared_sequences_dir: str):
        self.seed_list_path = seed_list_path
        self.used_seeds_path = used_seeds_path
        self.all_members_path = all_members_path
        self.sequence_identity_dir = sequence_identity_dir
        self.seed_to_hits_sequence_identity_path = seed_to_hits_sequence_identity_path
        self.shared_sequences_dir = shared_sequences_dir


    def get_results_from_finished_jobs(self, jobs: dict[int, dict]) -> dict[int, dict]:
        finished_jobs = [job_id for job_id in jobs.keys() if jobs[job_id]["status"] == JobStatus.COMPLETED]
        for job_id in finished_jobs:
            try:
                job_config_path = jobs[job_id]["config_path"]
                job_config = read_yaml_config(job_config_path)
                seed_id_of_job = job_config["run"]["seed_uid"]
                job_root_dir = job_config["global"]["root_dir"]
                local_results_dir = f"{job_root_dir}/results"

                sorted_merged_seeds, new_members = self.merge_local_true_members_with_global_seeds(local_results_dir)
                sorted_merged_seeds.to_csv(self.seed_list_path, index=False)

                updated_used_seeds_df = self.add_current_seed_info_to_used_seeds(seed_id_of_job, len(new_members))
                updated_used_seeds_df.to_csv(self.used_seeds_path, index=False)

                merged_all_members_global = self.merge_local_all_members_with_global_all_members(local_results_dir)
                merged_all_members_global.to_csv(self.all_members_path, index=False)
                
                self._compute_sequence_identity_to_seeds(updated_used_seeds_df, new_members, seed_id_of_job)
                
                logger.info("Collected results from job %d (seed %s): %d new true members.",
                            job_id, seed_id_of_job, len(new_members))
                jobs[job_id]["status"] = JobStatus.HARVESTED
            except Exception as e:
                logger.error("Error while collecting results from job %d: %s", job_id, e)
        return jobs
    

    def _compute_sequence_identity_to_seeds(self, used_seeds_df: pd.DataFrame, new_members_df: pd.DataFrame, seed_id_of_job: str):
        current_seed_fasta_path, new_members_fasta_path, previous_used_seeds_fasta_path = self._update_and_generate_fasta_files_for_sequence_identity(used_seeds_df, new_members_df, seed_id_of_job)
        sequence_identity_calculator = SequenceIdentityCalculator(work_dir=self.sequence_identity_dir)
        sequence_identity_calculator.compute_sequence_identity(
            output_csv_path=self.seed_to_hits_sequence_identity_path,
            current_seed_fasta_path=current_seed_fasta_path,
            new_members_fasta_path=new_members_fasta_path,
            previous_used_seeds_fasta_path=previous_used_seeds_fasta_path,
        )
        

    def _update_and_generate_fasta_files_for_sequence_identity(self, used_seeds_df: pd.DataFrame, new_members_df: pd.DataFrame, seed_id_of_job: str):
        # create fasta file for seed (temp)
        seed_record = get_records_from_accessions([seed_id_of_job], cache_dir=self.shared_sequences_dir)
        current_seed_fasta_path = os.path.join(self.sequence_identity_dir, "current_seed.fasta")
        SeqIO.write(seed_record, current_seed_fasta_path, "fasta")
        
        # create fasta file for new detected members (temp)
        new_members = new_members_df["id"].tolist()
        new_member_records = get_records_from_accessions(new_members, cache_dir=self.shared_sequences_dir)
        new_members_fasta_path = os.path.join(self.sequence_identity_dir, "new_members.fasta")
        SeqIO.write(new_member_records, new_members_fasta_path, "fasta")
        
        # fasta of previous used seeds (temp because small)
        previous_used_seeds = used_seeds_df[used_seeds_df["id"] != seed_id_of_job]
        previous_used_seed_records = get_records_from_accessions(previous_used_seeds["id"].tolist(), cache_dir=self.shared_sequences_dir)
        previous_used_seeds_fasta_path = os.path.join(self.sequence_identity_dir, "previous_used_seeds.fasta")
        SeqIO.write(previous_used_seed_records, previous_used_seeds_fasta_path, "fasta")
        
        return current_seed_fasta_path, new_members_fasta_path, previous_used_seeds_fasta_path


    def merge_local_all_members_with_global_all_members(self, local_results_dir: str) -> pd.DataFrame:
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


    def merge_local_true_members_with_global_seeds(self, local_results_dir: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        global_seeds_df = pd.read_csv(self.seed_list_path)
        local_true_members_df = pd.read_csv(f"{local_results_dir}/True_Members.csv")
        new_members = local_true_members_df[~local_true_members_df["id"].isin(global_seeds_df["id"])]
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

        return sorted_merged_seeds, new_members
    