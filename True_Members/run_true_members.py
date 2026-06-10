from scripts.utils import *
from scripts.evcouplings import *
from scripts.get_structures_and_sequences import *
from scripts.structural_filter import *
from scripts.hmm_domain_filter import *
import time 

CONFIG_TEMPLATE = "evcouplings_config_template.yaml"
number_tracker = {} # dictionary to keep track of true members at each stage in the pipeline

# output directory structure in runs directory:

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
        # results
            # number_tracker.csv
            # True_Members.csv
            # structures_with_seed_and_not_all_domains.csv
        
def read_sys_args():
        args = sys.argv
        if len(args) < 2 or len(args) > 2:
            print("Usage: python run_true_members.py <path_to_config_file>")
            sys.exit(1)
        config_path = sys.argv[1] if len(sys.argv) > 1 else print("Please provide the path to the config file as a command line argument.")
        config = read_true_members_config(config_path)
        return config
        
def read_true_members_config(config_path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

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
        f"{root_dir}/data/results"
    ]
    for directory in directories:
        # clear directory if it already exists
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.makedirs(directory, exist_ok=True)
    return root_dir

if __name__ == "__main__":
    
    config = read_sys_args()
    
    SEED_uid = config["global"]["seed_uid"]
    SEED_region_start = config["global"]["seed_region_start"]
    SEED_region_end = config["global"]["seed_region_end"]
    SEED = [SEED_uid, SEED_region_start, SEED_region_end]
    BIT_SCORES = config["global"]["bit_scores"]
    TM_SCORE_CUTOFF = config["global"]["tm_score_cutoff"]
    HMMs_DIR = config["global"]["hmms_folder_path"]
    CPUS = config["global"]["CPU"]
    
    root_dir = setup_directories()

    # Run evcouplings align and couplings stage
    run_evcouplings(SEED, BIT_SCORES, CPUS, CONFIG_TEMPLATE, root_dir)
    
    # Run foldseek webserver search and filter output with TM_SCORE_CUTOFF
    # get correct roi of seed
    seed_uid = SEED[0]
    seed_r_start = SEED[1]
    seed_r_end = SEED[2]
    run_download_tasks([seed_uid], output_dir=f"{root_dir}/data/seed/", max_workers=1)
    query_structure_full = f"{root_dir}/data/seed/{seed_uid}.cif"
    query_structure_roi = extract_region_of_protein(query_structure_full, (seed_r_start, seed_r_end))

    number_tracker = run_foldseek_webserver(number_tracker, query_structure_roi, root_dir, TM_SCORE_CUTOFF)

    # Get all structures and full fasta of all sequences of evcouplings output
    EVCOUPLINGS_OUT_DIR = f"{root_dir}/data/evcouplings"
    number_tracker = get_all_structures_and_sequences(SEED, BIT_SCORES, number_tracker, root_dir)

    # Run USalgin, foldseek and filter output
    number_tracker = run_structural_similarity_filter(SEED, number_tracker, TM_SCORE_CUTOFF, root_dir)

    # Run hmmsearch with CbrA related HMMs on fasta of STAC containing sequences and filter to get STAC containing CbrAs and STAC containing proteins without CbrA domains
    seed_containing_sequences_ids_path = f"{root_dir}/data/filtering/structural_filtering_output_tm_0.6.csv"
    number_tracker = run_hmm_domain_filter(number_tracker, CPUS, seed_containing_sequences_ids_path, HMMs_DIR, root_dir)

    # Print number of true members at each stage in the pipeline and save results in results directory
    for key, value in number_tracker.items():
        print(f"{key}: {value}")
        
    df_number_tracker = pd.DataFrame(list(number_tracker.items()), columns=["Stage", "Numbers"])
    df_number_tracker.to_csv(f"{root_dir}/data/results/number_tracker.csv", index=False)

    copy_file(f"{root_dir}/data/filtering/True_Members.csv", f"{root_dir}/data/results/True_Members.csv")
    copy_file(f"{root_dir}/data/filtering/structures_with_seed_and_not_all_domains.csv", f"{root_dir}/data/results/structures_with_seed_and_not_all_domains.csv")
    copy_file(f"{root_dir}/data/filtering/usalign/filtered_usalign_output_tm_0.6.csv", f"{root_dir}/data/results/filtered_usalign_output_tm_0.6.csv")
    copy_file(f"{root_dir}/data/filtering/hmm_domain_filtering/hmmsearch_results.domtblout", f"{root_dir}/data/results/hmmsearch_results.domtblout")
    