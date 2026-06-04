from scripts.utils import *
from scripts.setup_evcouplings import *
from scripts.get_structures_and_sequences import *
from scripts.structural_filter import *
from scripts.hmm_domain_filter import *
import time 


# IMPORTANT PARAMETERS TO SET BEFORE RUNNING:
# evcouplings, foldseek, USalign, hmmpress and hmmsearch should be installed and in Path

CPUS = 16
BIT_SCORES = [0.8]
SEED= ["Q5SK37", 507, 569] # [uid,region_start,region_end]
HMMs_DIR = "STAC_hmm_domains" # should contain GA cutoff
CONFIG_TEMPLATE = "evcouplings_config_template.yaml"
number_tracker = {} # dictionary to keep track of true members at each stage in the pipeline


# directory setup for each run:

# root_dir
    # data
        # evcouplings
            # contains config file and output
        # filtering
            # usalgin_dir
            # foldseek_dir
            # structural_filtering_output.csv
            # True_Members.csv
            # structures_containing_seeds_and_not_al_domains.csv
            # hmm_domain_filtering_dir
        # alphafold_structures
            # contains structures of potential true members
        # seed
            # contains structure of seed
        # HMMs
            # contains HMMs
        # all_sequences_full.fasta


def setup_directories():
    root_dir = f"runs/true_members_run_{SEED[0]}_{time.strftime('%Y-%m-%d_%H-%M-%S')}"
    directories = [
        f"{root_dir}/data/evcouplings",
        f"{root_dir}/data/filtering/usalign",
        f"{root_dir}/data/filtering/foldseek",
        f"{root_dir}/data/filtering/hmm_domain_filtering",
        f"{root_dir}/data/alphafold_structures",
        f"{root_dir}/data/seed",
        f"{root_dir}/data/HMMs",
    ]
    for directory in directories:
        # clear directory if it already exists
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.makedirs(directory, exist_ok=True)
    return root_dir
        
root_dir = setup_directories()


# Run evcouplings align and couplings stage
run_evcouplings(SEED, BIT_SCORES, CPUS, CONFIG_TEMPLATE, root_dir)


# Get all structures and full fasta of all sequences of evcouplings output
EVCOUPLINGS_OUT_DIR = f"{root_dir}/data/evcouplings",
number_tracker = get_all_structures_and_sequences(SEED, BIT_SCORES, number_tracker, root_dir)


# Run USalgin, foldseek and filter output
number_tracker = run_structural_similarity_filter(SEED, number_tracker, root_dir)


# Run hmmsearch with CbrA related HMMs on fasta of STAC containing sequences and filter to get STAC containing CbrAs and STAC containing proteins without CbrA domains
seed_containing_sequences_ids_path = f"{root_dir}/data/filtering/structural_filtering_output_tm_0.6.csv"
number_tracker = run_hmm_domain_filter(number_tracker, CPUS, seed_containing_sequences_ids_path, HMMs_DIR, root_dir)


# Print number of true members at each stage in the pipeline
for key, value in number_tracker.items():
    print(f"{key}: {value}")
    