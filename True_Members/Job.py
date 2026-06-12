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
    def __init__(self, true_members_config, ev_couplings_config_template, global_data_dir, root_dir = None):
        self.true_members_config = true_members_config
        self.global_data_dir = global_data_dir
        self.CONFIG_TEMPLATE = ev_couplings_config_template
        self.number_tracker = {}
        
        try:
            self.iteration = self.true_members_config["global"]["iteration"]
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
        
        self.setup_directories(root_dir)
        
    
    def setup_directories(self, root_dir):
        if root_dir is not None:
            self.root_dir = root_dir
        else:
            self.root_dir = f"runs/true_members_run_{self.SEED[0]}_{time.strftime('%Y-%m-%d_%H-%M-%S')}"
            
        self.evcouplings_dir = f"{self.root_dir}/data/evcouplings"
        self.usalign_dir = f"{self.root_dir}/data/filtering/usalign"
        self.foldseek_dir = f"{self.root_dir}/data/filtering/foldseek"
        self.hmm_domain_filtering_dir = f"{self.root_dir}/data/filtering/hmm_domain_filtering"
        self.alphafold_structures_dir = f"{self.root_dir}/data/alphafold_structures"
        self.seed_dir = f"{self.root_dir}/data/seed"
        self.results_dir = f"{self.root_dir}/data/results"
        
        # create directories
        if root_dir is None:
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
        #run_evcouplings(self.SEED, self.BIT_SCORES, self.CPUS, self.CONFIG_TEMPLATE, self.root_dir)
        
        # Run foldseek webserver search with ROI of SEED and filter output with TM_SCORE_CUTOFF
        run_download_tasks([self.SEED_uid], output_dir=self.seed_dir, max_workers=1)
        query_structure_full = f"{self.seed_dir}/{self.SEED_uid}.cif"
        query_structure_roi = extract_region_of_protein(query_structure_full, (self.SEED_region_start, self.SEED_region_end))
        run_foldseek_webserver(self.number_tracker, query_structure_roi, self.root_dir, self.TM_SCORE_CUTOFF, self.foldseek_dir)
        
        # Get all structures and full fasta of all sequences of evcouplings output
        seed_containing_sequences_full_fasta_path, seed_containing_sequences_with_region_path = get_all_structures_and_sequences(
            self.SEED,
            self.BIT_SCORES,
            self.number_tracker,
            self.evcouplings_dir,
            self.foldseek_dir,
            self.alphafold_structures_dir,
            self.root_dir)
        
        # Run structural filter with USalign and fildseek
        seed_roi_path = query_structure_roi
        run_structural_similarity_filter(seed_roi_path,
                                         self.number_tracker, 
                                         self.TM_SCORE_CUTOFF,
                                         self.alphafold_structures_dir,
                                         self.foldseek_dir,
                                         self.usalign_dir,
                                         seed_containing_sequences_with_region_path)
        
        # Run hmmsearch with HMMs in HMMs_dir on all_sequences_full.fasta and filter true members
        run_hmm_domain_filter(self.number_tracker,
                              self.CPUS,
                              seed_containing_sequences_full_fasta_path,
                              seed_containing_sequences_with_region_path,
                              self.HMMs_DIR,
                              self.hmm_domain_filtering_dir,
                              self.root_dir)
        
        self.print_and_save_results()
        
        
    def print_and_save_results(self):
        # Print number of true members at each stage in the pipeline and save results in results directory
        for key, value in self.number_tracker.items():
            print(f"{key}: {value}")
            
        df_number_tracker = pd.DataFrame(list(self.number_tracker.items()), columns=["Stage", "Numbers"])
        df_number_tracker.to_csv(f"{self.results_dir}/number_tracker.csv", index=False)

        # save to local results
        copy_file(f"{self.root_dir}/data/filtering/True_Members.csv", f"{self.results_dir}/True_Members.csv")
        copy_file(f"{self.root_dir}/data/filtering/All_Members.csv", f"{self.results_dir}/All_Members.csv")
        copy_file(f"{self.root_dir}/data/filtering/hmm_domain_filtering/hmmsearch_results.domtblout", f"{self.results_dir}/hmmsearch_results.domtblout")
        
        # save to global data (seeds, used seeds, all members)
        sorted_merged_seeds, new_true_members = self.merge_true_members_with_global_seeds()
        sorted_merged_seeds.to_csv(os.path.join(self.global_data_dir, "seed_list.csv"), index=False)
        updated_used_seeds_df = self.add_current_seed_to_used_seeds(new_true_members)
        updated_used_seeds_df.to_csv(os.path.join(self.global_data_dir, "used_seeds.csv"), index=False)
        merged_all_members_global = self.merge_all_members_with_global_all_members()
        merged_all_members_global.to_csv(os.path.join(self.global_data_dir, "All_Members_global.csv"), index=False)
        
        
    def merge_all_members_with_global_all_members(self) -> pd.DataFrame:
        all_members_df = pd.read_csv(f"{self.results_dir}/All_Members.csv")
        all_members_global_df = pd.read_csv(os.path.join(self.global_data_dir, "All_Members_global.csv"))
        
        dfs_combined = pd.concat([all_members_global_df, all_members_df], ignore_index=True)
        # identify hmm domain columns (all entries are true)
        bool_cols = [c for c in dfs_combined.columns if c not in ["id", "qtmscore", "region_start", "region_end"]]

        # aggregation rules for key columns
        agg_rules = {
            "qtmscore": "median",
            "region_start": "min", # use lowest start region found
            "region_end": "max" # use highest end region found
        }
        # add aggregation rules for hmm domain columns
        for col in bool_cols:
            agg_rules[col] = "any" # some entries might be true

        merged_all_members_global = dfs_combined.groupby("id", as_index=False).agg(agg_rules)
        return merged_all_members_global 
        
    
    def add_current_seed_to_used_seeds(self, new_true_members: int) -> pd.DataFrame:
        used_seeds_df = pd.read_csv(os.path.join(self.global_data_dir, "used_seeds.csv"))
        current_seed_row = pd.DataFrame([{
            "id": self.SEED[0], 
            "region_start": self.SEED[1],
            "region_end": self.SEED[2],
            "iteration": self.iteration,
            "new_true_members": new_true_members
            
        }])
        updated_used_seeds_df = pd.concat([used_seeds_df, current_seed_row], ignore_index=True)
        return updated_used_seeds_df
        
    
    def merge_true_members_with_global_seeds(self) -> tuple[pd.DataFrame, int]:  
        true_members_df = pd.read_csv(f"{self.results_dir}/True_Members.csv")
        seeds_df = pd.read_csv(os.path.join(self.global_data_dir, "seed_list.csv"))
        
        
        dfs_combined = pd.concat([seeds_df, true_members_df], ignore_index=True)
        # identify hmm domain columns (all entries are true)
        bool_cols = [c for c in dfs_combined.columns if c not in ["id", "qtmscore", "region_start", "region_end"]]

        # aggregation rules for key columns
        agg_rules = {
            "qtmscore": "median",
            "region_start": "min", # use lowest start region found
            "region_end": "max" # use highest end region found
        }
        # add aggregation rules for hmm domain columns
        for col in bool_cols:
            agg_rules[col] = "all" # all entries should be true

        merged_seeds = dfs_combined.groupby("id", as_index=False).agg(agg_rules)
        sorted_merged_seeds = merged_seeds.sort_values("qtmscore", ascending=False)
        
        # compute number of new true members found in this iteration
        seeds_length_before_merge = len(seeds_df)
        seeds_length_after_merge = len(sorted_merged_seeds)
        new_true_members = seeds_length_after_merge - seeds_length_before_merge
        
        return sorted_merged_seeds, new_true_members
                