import shutil
import subprocess
import pandas as pd
from scripts.utils import *
import sys
from Bio import SearchIO
import os
        

def run_hmm_domain_filter(number_tracker: dict, cpus: int, seed_containing_sequences_ids_path: str, HMMs_DIR: str, root_dir: str) -> dict:
    # press hmms in HMMs dir
    domain_validator_hmm_path = press_hmms(HMMs_DIR, root_dir)
    
    # grab fasta records of all seed containing sequences
    all_sequences_fasta_path = f"{root_dir}/data/all_sequences_full.fasta"    
    seed_containing_sequences_ids = pd.read_csv(seed_containing_sequences_ids_path)["target"].tolist()
    seed_containing_sequences_fasta_path = filter_fasta_file(all_sequences_fasta_path, seed_containing_sequences_ids, f"{root_dir}/data/filtering/hmm_domain_filtering/seed_containing_sequences.fasta")
    
    # run hmmsearch to check if desired domains are present in the seed containing sequences
    hmmsearch_domtblout_output = run_hmmsearch(domain_validator_hmm_path, seed_containing_sequences_fasta_path, f"{root_dir}/data/filtering/hmm_domain_filtering/hmmsearch_results", cpus)
    
    # filter hmmsearch results to get seed containing proteins with domains and seed containing proteins without domains
    hmms_with_ga_cutoffs = get_GA_cutoffs_from_hmms(HMMs_DIR)
    hmmsearch_results = parse_hmmsearch_results(hmmsearch_domtblout_output, hmms_with_ga_cutoffs)
    hmms_as_list = list(hmms_with_ga_cutoffs.keys())
    number_tracker = filter_hmmsearch_results(hmmsearch_results, seed_containing_sequences_ids, number_tracker, hmms_as_list, root_dir)
    
    return number_tracker
    
    
def parse_hmmsearch_results(hmmsearch_domtblout_output: str, hmms_with_ga_cutoffs: dict) -> list[dict]:
    data = [] # consists of {"id": sequence_id, "hmm": hmm_name, "contained": is_hit_valid}
    for qresult in SearchIO.parse(hmmsearch_domtblout_output, "hmmsearch3-domtab"):
        hmm_name = qresult.id        
        for hit in qresult.hits:
            sequence_id = hit.id
            cutoff = hmms_with_ga_cutoffs.get(hmm_name, 25.0)
            is_hit_valid = any(hsp.bitscore >= cutoff for hsp in hit.hsps)
            data.append({"id": extract_uid_from_recordID(sequence_id), "hmm": hmm_name, "contained": is_hit_valid})
            
    return data


def filter_hmmsearch_results(data: list[dict], seed_containing_sequences_ids: list[str], number_tracker: dict, expected_hmms: list[str], root_dir: str) -> dict:
    if data:
        df = pd.DataFrame(data)
        seed_df = df.pivot_table(
            index="id", 
            columns="hmm", 
            values="contained", 
            fill_value=False
        ).astype(bool)
        
        seed_df = seed_df.reindex(columns=expected_hmms, fill_value=False)
        
        # add sequences which contain seed have no domain hit (not contained in hmmsearch results)
        seed_df = seed_df.reindex(seed_containing_sequences_ids, fill_value=False)
        
        # Filter seed containing sequences in contains_all_domains and not_contains_all_domains
        seed_df["contains_all_domains"] = seed_df.all(axis=1)
        structures_with_seed_and_domains = seed_df[seed_df["contains_all_domains"]].reset_index()
        number_tracker["structures_with_seed_and_domains"] = len(structures_with_seed_and_domains)
        
        structures_with_seed_and_not_all_domains = seed_df[~seed_df["contains_all_domains"]].reset_index()
        number_tracker["structures_with_seed_and_not_all_domains"] = len(structures_with_seed_and_not_all_domains)
        
        # Save output
        structures_with_seed_and_domains.to_csv(f"{root_dir}/data/filtering/True_Members.csv", index=False, header=True)
        structures_with_seed_and_not_all_domains.to_csv(f"{root_dir}/data/filtering/structures_with_seed_and_not_all_domains.csv", index=False, header=True)
    
    else:
        print("No valid domains passed hits found based on the E-value cutoff.")
        number_tracker["structures_with_seed_and_domains"] = 0
        number_tracker["structures_with_seed_and_not_all_domains"] = len(seed_containing_sequences_ids)
    
    return number_tracker


def get_GA_cutoffs_from_hmms(HMMs_DIR: str) -> dict:
    GA_cutoffs = {}
    for filename in os.listdir(HMMs_DIR):
        if filename.endswith(".hmm") and filename != "domain_validator.hmm":
            with open(os.path.join(HMMs_DIR, filename), "r") as f:
                name = None
                ga_domain_score = None
                for line in f:
                    clean_line = line.strip()
                    if clean_line.startswith("NAME"):
                        name = parse_name(clean_line)
                    if clean_line.startswith("GA"):
                        ga_domain_score = parse_ga(clean_line)
                    if name and ga_domain_score is not None:
                        GA_cutoffs[name] = ga_domain_score
                        break
    return GA_cutoffs


def parse_name(clean_line: str):
    parts = clean_line.split()
    if len(parts) >= 2:
        return parts[1]
    else:
        return None


def parse_ga(clean_line: str):
    parts = clean_line.rstrip(';').split()
    if len(parts) >= 3:
        domain_cutoff = parts[2].rstrip(';')
        return float(domain_cutoff)
    else:
        return None
    

def run_hmmsearch(hmm_path: str, fasta_path: str, output_prefix: str, cpus: int) -> str:
    try:
        print(f"Running hmmsearch with HMM: {hmm_path} on fasta: {fasta_path}...")
        subprocess.run(
            ["hmmsearch",
                        "--cpu", 
                        str(cpus),
                        "--tblout",
                        f"{output_prefix}.tbl",
                        "--domtblout",
                        f"{output_prefix}.domtblout",
                        hmm_path,
                        fasta_path],
            check=True, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            text=True)
        
        return f"{output_prefix}.domtblout"
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"hmmsearch failed with error:\n{e.stderr}") from e
    
    
def press_hmms(hmm_dir: str, root_dir: str) -> str:
    output_hmm_path = f"{root_dir}/data/HMMs/domain_validator.hmm"
    hmm_files = [f for f in os.listdir(hmm_dir) if f != "domain_validator.hmm" and f.endswith(".hmm")]

    with open(output_hmm_path, 'wb') as outfile:
        for hmm_file in hmm_files:
            with open(os.path.join(hmm_dir, hmm_file), 'rb') as infile:
                shutil.copyfileobj(infile, outfile)
    
    try:
        subprocess.run(
            ["hmmpress", str(output_hmm_path)],
            check=True,                 # Raises an exception if hmmpress fails
            stdout=subprocess.PIPE,     # Captures standard output
            stderr=subprocess.PIPE,     # Captures standard error
            text=True                   # Returns strings instead of bytes
        )
    except subprocess.CalledProcessError as e:
        # Provide a highly descriptive error if the HMMER tool fails
        raise RuntimeError(f"hmmpress failed with error:\n{e.stderr}") from e
    
    return output_hmm_path
                