import time
import os
import shutil
from scripts.utils import *
from scripts.evcouplings import *
from scripts.get_structures_and_sequences import *
from scripts.structural_filter import *
from scripts.hmm_domain_filter import *
from scripts.evcouplings import *

class Job:
    def __init__(self, true_members_config_path: str, ev_couplings_config_template: str, use_existing_root_dir: bool = False):
        self.true_members_config = read_yaml_config(true_members_config_path)
        self.CONFIG_TEMPLATE = ev_couplings_config_template
        self.number_tracker = {}
        
        try:
            self.iteration = self.true_members_config["global"]["iteration"]
            self.root_dir = self.true_members_config["global"]["root_dir"]
            self.SEED_uid = self.true_members_config["run"]["seed_uid"]
            self.SEED_region_start = self.true_members_config["run"]["seed_region_start"]
            self.SEED_region_end = self.true_members_config["run"]["seed_region_end"]
            self.SEED = [self.SEED_uid, self.SEED_region_start, self.SEED_region_end]
            self.BIT_SCORES = self.true_members_config["run"]["bit_scores"]
            self.TM_SCORE_CUTOFF = self.true_members_config["run"]["tm_score_cutoff"]
            self.HMMs_DIR = self.true_members_config["run"]["hmms_folder_path"] # holds HMMs specifies by user
            self.CPUS = self.true_members_config["run"]["CPU"]
        except KeyError as e:
            print(f"Missing key in config file: {e}")
            raise
        
        self.setup_directories(use_existing_root_dir)
        
    
    def setup_directories(self, use_existing_root_dir: bool):
            
        self.evcouplings_dir = f"{self.root_dir}/evcouplings"
        self.usalign_dir = f"{self.root_dir}/filtering/usalign"
        self.foldseek_dir = f"{self.root_dir}/filtering/foldseek"
        self.hmm_domain_filtering_dir = f"{self.root_dir}/filtering/hmm_domain_filtering"
        self.alphafold_structures_dir = f"{self.root_dir}/alphafold_structures"
        self.seed_dir = f"{self.root_dir}/seed"
        self.results_dir = f"{self.root_dir}/results"
        
        # create directories
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
                # clear directory if it already exists
                if os.path.exists(directory):
                    shutil.rmtree(directory)
                os.makedirs(directory, exist_ok=True)
    
    
    def run(self):          
        # Run evcouplings align and couplings stage
        # run_evcouplings(self.SEED, self.BIT_SCORES, self.CPUS, self.CONFIG_TEMPLATE, self.root_dir)
        
        # Run foldseek webserver search with ROI of SEED and filter output with TM_SCORE_CUTOFF
        run_download_tasks([self.SEED_uid], output_dir=self.seed_dir, max_workers=1)
        query_structure_full = f"{self.seed_dir}/{self.SEED_uid}.cif"
        query_structure_roi = extract_region_of_protein(query_structure_full, (self.SEED_region_start, self.SEED_region_end))
        run_foldseek_webserver(self.number_tracker, query_structure_roi, self.root_dir, self.TM_SCORE_CUTOFF, self.foldseek_dir)
        
        # Get all fasta sequences and structures of evcouplings and foldseek output
        uids_with_region_full_fasta_path, uids_with_region = get_all_structures_and_sequences(
            self.SEED,
            self.BIT_SCORES,
            self.number_tracker,
            self.evcouplings_dir,
            self.foldseek_dir,
            self.alphafold_structures_dir,
            self.root_dir)
        
        # Run structural filter with USalign and foldseek
        seed_containing_sequences_with_tm = run_structural_similarity_filter(query_structure_roi,
                                         self.number_tracker, 
                                         self.TM_SCORE_CUTOFF,
                                         self.alphafold_structures_dir,
                                         self.foldseek_dir,
                                         self.usalign_dir)
        
        # Run hmmsearch domain filter
        hmms_with_ga_cutoffs = get_hmms_with_ga_cutoffs(self.HMMs_DIR)
        hmm_search_results = run_hmm_domain_filter(
            self.number_tracker,
            self.CPUS,
            uids_with_region_full_fasta_path,
            seed_containing_sequences_with_tm,
            self.HMMs_DIR,
            self.hmm_domain_filtering_dir,
            self.root_dir,
            self.TM_SCORE_CUTOFF,
            hmms_with_ga_cutoffs)
        
        # combine results of structural and domain filter
        all_members, true_members = self.combine_results(uids_with_region, seed_containing_sequences_with_tm, hmm_search_results, hmms_with_ga_cutoffs)
        self.print_and_save_results(all_members, true_members)
        
        
    def combine_results(self, uids_with_region: pd.DataFrame, seed_containing_sequences_with_tm: pd.DataFrame, hmm_search_results: pd.DataFrame, hmms_with_ga_cutoffs: dict):
        combined_results = uids_with_region.merge(seed_containing_sequences_with_tm, on="id", how="left")
        combined_results = combined_results.dropna(subset=["qtmscore"])
        combined_results = pd.merge(combined_results, hmm_search_results, on="id", how="left")
        expected_hmms = list(hmms_with_ga_cutoffs.keys())
        combined_results[expected_hmms] = combined_results[expected_hmms].fillna(False).astype(bool)
        
        all_members = combined_results
        self.number_tracker["all_members"] = len(all_members)
        true_members = all_members[all_members[expected_hmms].all(axis=1)]
        self.number_tracker["true_members"] = len(true_members)
        return all_members, true_members
            
        
    def print_and_save_results(self, all_members: pd.DataFrame, true_members: pd.DataFrame):
        # Print number of true members at each stage in the pipeline and save results in results directory
        for key, value in self.number_tracker.items():
            print(f"{key}: {value}")
            
        df_number_tracker = pd.DataFrame(list(self.number_tracker.items()), columns=["Stage", "Numbers"])
        df_number_tracker.to_csv(f"{self.results_dir}/number_tracker.csv", index=False)

        # save to local results
        all_members.to_csv(f"{self.results_dir}/All_Members.csv", index=False)
        true_members.to_csv(f"{self.results_dir}/True_Members.csv", index=False)
        