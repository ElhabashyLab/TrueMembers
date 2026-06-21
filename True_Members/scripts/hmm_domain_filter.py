import logging
import os
import shutil
import subprocess

import pandas as pd
from Bio import SearchIO

from scripts.types import JobStats
from scripts.utils import extract_uid_from_recordID, get_records_from_fasta, write_fasta_from_records

logger = logging.getLogger(__name__)


def run_hmm_domain_filter(stats: JobStats, cpus: int, uids_with_region_fasta_path: str, seed_containing_sequences: pd.DataFrame, HMMs_DIR: str, hmm_domain_filtering_dir: str, root_dir: str, TM_SCORE_CUTOFF: float, hmms_with_ga_cutoffs: dict) -> pd.DataFrame:
    domain_validator_hmm_path = press_hmms(hmm_domain_filtering_dir, HMMs_DIR)

    records = get_records_from_fasta(uids_with_region_fasta_path)
    seed_containing_records = [record for record in records if extract_uid_from_recordID(record.id) in seed_containing_sequences["id"].tolist()]
    seed_containing_sequences_full_fasta = f"{root_dir}/filtering/seed_containing_sequences_tm_{TM_SCORE_CUTOFF}.fasta"
    write_fasta_from_records(seed_containing_records, seed_containing_sequences_full_fasta)

    hmmsearch_domtblout_output = run_hmmsearch(domain_validator_hmm_path, seed_containing_sequences_full_fasta, f"{hmm_domain_filtering_dir}/hmmsearch_results", cpus)

    hmmsearch_results = parse_hmmsearch_results(hmmsearch_domtblout_output, hmms_with_ga_cutoffs)
    stats.hmm_search_results = len(hmmsearch_results)
    return hmmsearch_results


def parse_hmmsearch_results(hmmsearch_domtblout_output: str, hmms_with_ga_cutoffs: dict) -> pd.DataFrame:
    data = []
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

    hmm_search_results = pd.DataFrame()
    if data:
        df = pd.DataFrame(data)
        hmm_search_results = df.pivot_table(
            index="id",
            columns="hmm",
            values="contained",
            fill_value=False
        ).astype(bool)

        expected_hmms = list(hmms_with_ga_cutoffs.keys())
        hmm_search_results = hmm_search_results.reindex(columns=expected_hmms, fill_value=False).reset_index()
    else:
        logger.info("No valid domains found in hmmsearch result.")

    return hmm_search_results


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
                    logger.warning("Could not find NAME line in %s. Skipping this HMM.", filename)
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
        logger.info("Running hmmsearch with HMM: %s on fasta: %s...", hmm_path, fasta_path)
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