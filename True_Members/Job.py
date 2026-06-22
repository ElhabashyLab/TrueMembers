import logging
import os
import shutil
from dataclasses import asdict

import pandas as pd

from scripts.evcouplings import run_evcouplings
from scripts.get_structures_and_sequences import get_all_structures_and_sequences
from scripts.hmm_domain_filter import get_hmms_with_ga_cutoffs, run_hmm_domain_filter
from scripts.structural_filter import run_foldseek_webserver, run_structural_similarity_filter
from scripts.types import JobStats
from scripts.utils import extract_region_of_protein, read_yaml_config, run_download_tasks

logger = logging.getLogger(__name__)


class Job:
    def __init__(self, true_members_config_path: str, ev_couplings_config_template: str, use_existing_root_dir: bool = False):
        self.true_members_config = read_yaml_config(true_members_config_path)
        self.CONFIG_TEMPLATE = ev_couplings_config_template
        self.stats = JobStats()

        try:
            self.GLOBAL_SEED_structure_roi = self.true_members_config["global"]["roi_structure_of_initial_seed"]
            self.iteration = self.true_members_config["global"]["iteration"]
            self.root_dir = self.true_members_config["global"]["root_dir"]
            self.SEED_uid = self.true_members_config["run"]["seed_uid"]
            self.SEED_region_start = self.true_members_config["run"]["seed_region_start"]
            self.SEED_region_end = self.true_members_config["run"]["seed_region_end"]
            self.SEED = [self.SEED_uid, self.SEED_region_start, self.SEED_region_end]
            self.BIT_SCORES = self.true_members_config["run"]["bit_scores"]
            self.TM_SCORE_CUTOFF = self.true_members_config["run"]["tm_score_cutoff"]
            self.HMMs_DIR = self.true_members_config["run"]["hmms_folder_path"]
            self.CPUS = self.true_members_config["run"]["CPU"]
        except KeyError as e:
            logger.error("Missing key in config file: %s", e)
            raise

        self._setup_logging()
        self.setup_directories(use_existing_root_dir)

    def _setup_logging(self):
        root_logger = logging.getLogger()
        root_logger.setLevel(logging.INFO)
        root_logger.handlers.clear()
        formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s")

        fh = logging.FileHandler(os.path.join(self.root_dir, "job.log"))
        fh.setFormatter(formatter)
        root_logger.addHandler(fh)

        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        root_logger.addHandler(ch)

    def setup_directories(self, use_existing_root_dir: bool):
        self.evcouplings_dir = f"{self.root_dir}/evcouplings"
        self.usalign_dir = f"{self.root_dir}/filtering/usalign"
        self.foldseek_dir = f"{self.root_dir}/filtering/foldseek"
        self.hmm_domain_filtering_dir = f"{self.root_dir}/filtering/hmm_domain_filtering"
        self.alphafold_structures_dir = f"{self.root_dir}/alphafold_structures"
        self.seed_dir = f"{self.root_dir}/seed"
        self.results_dir = f"{self.root_dir}/results"

        if not use_existing_root_dir:
            directories = [
                self.evcouplings_dir,
                self.usalign_dir,
                self.foldseek_dir,
                self.hmm_domain_filtering_dir,
                self.alphafold_structures_dir,
                self.seed_dir,
                self.results_dir
            ]

            for directory in directories:
                if os.path.exists(directory):
                    shutil.rmtree(directory)
                os.makedirs(directory, exist_ok=True)

    def run(self):
        logger.info("Starting job for seed %s, iteration %s", self.SEED_uid, self.iteration)

        run_evcouplings(self.SEED, self.BIT_SCORES, self.CPUS, self.CONFIG_TEMPLATE, self.root_dir)

        run_download_tasks([self.SEED_uid], output_dir=self.seed_dir, max_workers=1)
        query_structure_full = f"{self.seed_dir}/{self.SEED_uid}.cif"
        current_seed_roi = extract_region_of_protein(query_structure_full, (self.SEED_region_start, self.SEED_region_end))
        run_foldseek_webserver(self.stats, current_seed_roi, self.root_dir, self.TM_SCORE_CUTOFF, self.foldseek_dir)

        uids_with_region_full_fasta_path, uids_with_region = get_all_structures_and_sequences(
            self.SEED,
            self.BIT_SCORES,
            self.stats,
            self.evcouplings_dir,
            self.foldseek_dir,
            self.alphafold_structures_dir,
            self.root_dir)

        seed_containing_sequences_with_tm = run_structural_similarity_filter(self.GLOBAL_SEED_structure_roi,
                                         self.stats,
                                         self.TM_SCORE_CUTOFF,
                                         self.alphafold_structures_dir,
                                         self.foldseek_dir,
                                         self.usalign_dir)

        hmms_with_ga_cutoffs = get_hmms_with_ga_cutoffs(self.HMMs_DIR)
        hmm_search_results = run_hmm_domain_filter(
            self.stats,
            self.CPUS,
            uids_with_region_full_fasta_path,
            seed_containing_sequences_with_tm,
            self.HMMs_DIR,
            self.hmm_domain_filtering_dir,
            self.root_dir,
            self.TM_SCORE_CUTOFF,
            hmms_with_ga_cutoffs)

        all_members, true_members = self.combine_results(uids_with_region, seed_containing_sequences_with_tm, hmm_search_results, hmms_with_ga_cutoffs)
        self.print_and_save_results(all_members, true_members)

    def combine_results(self, uids_with_region: pd.DataFrame, seed_containing_sequences_with_tm: pd.DataFrame, hmm_search_results: pd.DataFrame, hmms_with_ga_cutoffs: dict):
        combined_results = uids_with_region.merge(seed_containing_sequences_with_tm, on="id", how="left")
        combined_results = combined_results.dropna(subset=["qtmscore"])
        combined_results = pd.merge(combined_results, hmm_search_results, on="id", how="left")
        expected_hmms = list(hmms_with_ga_cutoffs.keys())
        for col in expected_hmms:
            if col not in combined_results.columns:
                combined_results[col] = False
        combined_results[expected_hmms] = combined_results[expected_hmms].fillna(False).astype(bool)

        all_members = combined_results
        self.stats.all_members = len(all_members)
        true_members = all_members[all_members[expected_hmms].all(axis=1)]
        self.stats.true_members = len(true_members)
        return all_members, true_members

    def print_and_save_results(self, all_members: pd.DataFrame, true_members: pd.DataFrame):
        stats_dict = asdict(self.stats)
        for key, value in stats_dict.items():
            logger.info("%s: %s", key, value)

        df_number_tracker = pd.DataFrame(list(stats_dict.items()), columns=["Stage", "Numbers"])
        df_number_tracker.to_csv(f"{self.results_dir}/number_tracker.csv", index=False)

        all_members.to_csv(f"{self.results_dir}/All_Members.csv", index=False)
        true_members.to_csv(f"{self.results_dir}/True_Members.csv", index=False)

        logger.info("Results saved to %s", self.results_dir)
        