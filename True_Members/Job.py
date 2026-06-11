import time
import yaml
import os
import shutil
from scripts.utils import *
from scripts.evcouplings import *
from scripts.get_structures_and_sequences import *
from scripts.structural_filter import *
from scripts.hmm_domain_filter import *
from scripts.evcouplings import *

class Job:
    def __init__(self, true_members_config: str, root_dir = None):
        self.true_members_config = self.read_true_members_config(true_members_config)
        
        self.CONFIG_TEMPLATE = "evcouplings_config_template.yaml"
        self.number_tracker = {}
        
        self.SEED_uid = self.true_members_config["global"]["seed_uid"]
        self.SEED_region_start = self.true_members_config["global"]["seed_region_start"]
        self.SEED_region_end = self.true_members_config["global"]["seed_region_end"]
        self.SEED = [self.SEED_uid, self.SEED_region_start, self.SEED_region_end]
        self.BIT_SCORES = self.true_members_config["global"]["bit_scores"]
        self.TM_SCORE_CUTOFF = self.true_members_config["global"]["tm_score_cutoff"]
        self.HMMs_DIR = self.true_members_config["global"]["hmms_folder_path"] # holds HMMs specifies by user
        self.CPUS = self.true_members_config["global"]["CPU"]
        
        self.setup_directories(root_dir)
        
        
    def read_true_members_config(self, config_path):
        with open(config_path, 'r') as file:
            config = yaml.safe_load(file)
        return config
    
    
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
        # run_evcouplings(self.SEED, self.BIT_SCORES, self.CPUS, self.CONFIG_TEMPLATE, self.root_dir)
        
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

        copy_file(f"{self.root_dir}/data/filtering/True_Members.csv", f"{self.results_dir}/True_Members.csv")
        copy_file(f"{self.root_dir}/data/filtering/structures_with_seed_and_not_all_domains.csv", f"{self.results_dir}/structures_with_seed_and_not_all_domains.csv")
        copy_file(f"{self.root_dir}/data/filtering/hmm_domain_filtering/hmmsearch_results.domtblout", f"{self.results_dir}/hmmsearch_results.domtblout")
                