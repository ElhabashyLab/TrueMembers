import shutil
import subprocess
import pandas as pd
from scripts.utils import *
import sys
from Bio import SearchIO
import os


def run_hmm_domain_filter(number_tracker: dict, cpus: int, seed_containing_sequences_full_fasta: str, seed_containing_sequences_ids: list[str], HMMs_DIR: str, hmm_domain_filtering_dir: str, root_dir: str):
    # press hmms of HMMs_DIR in hmm_domain_filtering_dir
    domain_validator_hmm_path = press_hmms(hmm_domain_filtering_dir, HMMs_DIR)

    # run hmmsearch to check if desired domains are present in the seed containing sequences
    hmmsearch_domtblout_output = run_hmmsearch(domain_validator_hmm_path, seed_containing_sequences_full_fasta, f"{hmm_domain_filtering_dir}/hmmsearch_results", cpus)

    # filter hmmsearch results to get seed containing proteins with domains and seed containing proteins without domains
    hmms_with_ga_cutoffs = get_hmms_with_ga_cutoffs(HMMs_DIR)
    hmmsearch_results = parse_hmmsearch_results(hmmsearch_domtblout_output, hmms_with_ga_cutoffs)
    hmms_as_list = list(hmms_with_ga_cutoffs.keys())
    filter_hmmsearch_results(hmmsearch_results, seed_containing_sequences_ids, number_tracker, hmms_as_list, root_dir)


def parse_hmmsearch_results(hmmsearch_domtblout_output: str, hmms_with_ga_cutoffs: dict) -> list[dict]:
    data = [] # consists of {"id": sequence_id, "hmm": hmm_name, "contained": is_hit_valid}
    for qresult in SearchIO.parse(hmmsearch_domtblout_output, "hmmsearch3-domtab"):
        hmm_name = qresult.id
        for hit in qresult.hits:
            sequence_id = hit.id
            cutoff = hmms_with_ga_cutoffs.get(hmm_name)
            if cutoff is None:
                valid_hit = any(hsp.evalue <= 1e-3 for hsp in hit.hsps)
            else:
                valid_hit = any(hsp.bitscore >= cutoff for hsp in hit.hsps)
            data.append({"id": extract_uid_from_recordID(sequence_id), "hmm": hmm_name, "contained": valid_hit})

    return data


def filter_hmmsearch_results(data: list[dict], seed_containing_sequences_ids: list[str], number_tracker: dict, expected_hmms: list[str], root_dir: str):
    if data:
        df = pd.DataFrame(data)
        hmm_search_results = df.pivot_table(
            index="id",
            columns="hmm",
            values="contained",
            fill_value=False
        ).astype(bool)

        # fill hits for hmms which were not detected with false
        hmm_search_results = hmm_search_results.reindex(columns=expected_hmms, fill_value=False)

        # add sequences which contain seed but have no domain hit (not contained in hmmsearch results)
        seed_containing_structures = hmm_search_results.reindex(seed_containing_sequences_ids, fill_value=False)

        # Filter seed containing sequences in contains_all_domains and not_contains_all_domains
        seed_containing_structures["contains_all_domains"] = seed_containing_structures.all(axis=1)
        true_members = seed_containing_structures[seed_containing_structures["contains_all_domains"]].reset_index()
        number_tracker["structures_with_seed_and_domains"] = len(true_members)
        seed_containing_structures = seed_containing_structures.reset_index()
        number_tracker["structures_with_seed_and_not_all_domains"] = len(seed_containing_structures)

        # Save output
        true_members.to_csv(f"{root_dir}/data/filtering/True_Members.csv", index=False, header=True)
        seed_containing_structures.to_csv(f"{root_dir}/data/filtering/structures_with_seed_and_not_all_domains.csv", index=False, header=True)

    else:
        print("No valid domains found in hmmsearch result.")
        number_tracker["structures_with_seed_and_domains"] = 0
        number_tracker["structures_with_seed_and_not_all_domains"] = len(seed_containing_sequences_ids)


def get_hmms_with_ga_cutoffs(HMMs_DIR: str) -> dict:
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
                    elif clean_line.startswith("GA"):
                        ga_domain_score = parse_ga(clean_line)
                    elif clean_line.startswith("HMM ") or clean_line == "HMM":
                        break
                if name:
                    GA_cutoffs[name] = ga_domain_score
                else:
                    print(f"Warning: Could not find NAME line in {filename}. Skipping this HMM.")
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


def press_hmms(hmm_domain_filtering_dir: str, HMMs_DIR: str) -> str:
    domain_validator = f"{hmm_domain_filtering_dir}/domain_validator.hmm"
    hmm_files = [f for f in os.listdir(HMMs_DIR) if f.endswith(".hmm")]

    # copy all HMMs into domain_validator.hmm
    with open(domain_validator, 'wb') as outfile:
        for hmm_file in hmm_files:
            with open(os.path.join(HMMs_DIR, hmm_file), 'rb') as infile:
                shutil.copyfileobj(infile, outfile)

    try:
        subprocess.run(
            ["hmmpress", str(domain_validator)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"hmmpress failed with error:\n{e.stderr}") from e

    return domain_validator
