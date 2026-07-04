"""
The goal of this script is to create a summary of a TrueMembers pipeline run.
It summarizes the results of all jobs into a dict and saves it as JSON.

Usage:
    python create_run_summary.py <run_dir> <output_json_path>

The <run_dir> is a pipeline run directory (e.g. runs/run_Q5SK37_2026-06-27_...).

Data sources:
  - run_dir/final_true_members_results/Used_Seeds_i_*  → iteration, new_true_members
  - run_dir/runs/job_{seed_id}_iteration_{iteration}_*/results/number_tracker.csv
    → all pipeline stage counts
"""

import csv
import json
import os
import sys


def find_final_all_members_file(run_dir):
    """Return the path to the All_Members_i_* CSV file under final_true_members_results."""
    results_dir = os.path.join(run_dir, "final_true_members_results")
    if not os.path.isdir(results_dir):
        raise FileNotFoundError(f"Missing directory: {results_dir}")

    for entry in os.listdir(results_dir):
        if entry.startswith("All_Members_i_"):
            return os.path.join(results_dir, entry)

    raise FileNotFoundError(f"No All_Members_i_* file found in {results_dir}")


def find_final_true_members_file(run_dir):
    """Return the path to the True_Members_i_* CSV file under final_true_members_results."""
    results_dir = os.path.join(run_dir, "final_true_members_results")
    if not os.path.isdir(results_dir):
        raise FileNotFoundError(f"Missing directory: {results_dir}")

    for entry in os.listdir(results_dir):
        if entry.startswith("True_Members_i_"):
            return os.path.join(results_dir, entry)

    raise FileNotFoundError(f"No True_Members_i_* file found in {results_dir}")


def find_used_seeds_file(run_dir):
    """Return the path to the Used_Seeds_i_* CSV file under final_true_members_results."""
    results_dir = os.path.join(run_dir, "final_true_members_results")
    if not os.path.isdir(results_dir):
        raise FileNotFoundError(f"Missing directory: {results_dir}")

    for entry in os.listdir(results_dir):
        if entry.startswith("Used_Seeds_i_"):
            return os.path.join(results_dir, entry)

    raise FileNotFoundError(f"No Used_Seeds_i_* file found in {results_dir}")


def find_job_directory(run_dir, seed_id, iteration):
    """Return the job directory matching job_{seed_id}_iteration_{iteration}_* under run_dir/runs/."""
    runs_dir = os.path.join(run_dir, "runs")
    if not os.path.isdir(runs_dir):
        raise FileNotFoundError(f"Missing directory: {runs_dir}")

    prefix = f"job_{seed_id}_iteration_{iteration}_"
    for entry in os.listdir(runs_dir):
        if entry.startswith(prefix):
            return os.path.join(runs_dir, entry)

    raise FileNotFoundError(
        f"No job directory matching '{prefix}*' found in {runs_dir}"
    )


def read_number_tracker(job_dir):
    """Read number_tracker.csv and return a dict of {stage: value}."""
    path = os.path.join(job_dir, "results", "number_tracker.csv")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Missing number_tracker.csv at {path}")

    data = {}
    with open(path, newline="") as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            if len(row) < 2:
                continue
            stage = row[0].strip()
            try:
                data[stage] = int(row[1])
            except ValueError:
                data[stage] = row[1].strip()
    return data


def count_csv_lines(file_path):
    """
    Counts the number of data rows in a CSV file, excluding the header.
    Efficient for large files because it doesn't load the file into memory.
    """
    with open(file_path, 'r', encoding='utf-8') as f:
        # Skip the header line safely
        next(f, None) 
        
        # Count the remaining lines efficiently
        return sum(1 for line in f)


def build_summary(run_dir):
    """Iterate over Used_Seeds_i_* rows, collect per-job data, return summary dict."""
    used_seeds_path = find_used_seeds_file(run_dir)
    final_true_members_path = find_final_true_members_file(run_dir)
    final_all_members_path = find_final_all_members_file(run_dir)

    summary = {}
    with open(used_seeds_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            seed_id = row["id"].strip()
            iteration = int(row["iteration"])
            new_true_members_raw = row.get("new_true_members", "0").strip()
            try:
                new_true_members = int(float(new_true_members_raw))
            except ValueError:
                new_true_members = 0

            job_dir = find_job_directory(run_dir, seed_id, iteration)
            tracker = read_number_tracker(job_dir)

            summary[seed_id] = {
                "iteration": iteration,
                "new_true_members": new_true_members,
                "annotation_files_jackhmmer": tracker.get("annotation_files_jackhmmer", 0),
                "ids_jackhmmer": tracker.get("ids_jackhmmer", 0),
                "ids_foldseek_webserver": tracker.get("ids_foldseek_webserver", 0),
                "ids_combined": tracker.get("ids_combined", 0),
                "structure_downloads": tracker.get("structure_downloads", 0),
                "uniprot_records": tracker.get("uniprot_records", 0),
                "filtered_structures_foldseek": tracker.get("filtered_structures_foldseek", 0),
                "filtered_structures_usalign": tracker.get("filtered_structures_usalign", 0),
                "filtered_structures_total": tracker.get("filtered_structures_total", 0),
                "hmm_search_results": tracker.get("hmm_search_results", 0),
                "all_members": tracker.get("all_members", 0),
                "true_members": tracker.get("true_members", 0),
            }
            
        # add final true members count
        summary["final_true_members_count"] = count_csv_lines(final_true_members_path)
        # add final all members count
        summary["final_all_members_count"] = count_csv_lines(final_all_members_path)

    return summary


def main():
    if len(sys.argv) != 3:
        print("Usage: python create_run_summary.py <run_dir> <output_json_path>", file=sys.stderr)
        sys.exit(1)

    run_dir = sys.argv[1]
    output_path = sys.argv[2]

    if not os.path.isdir(run_dir):
        print(f"Error: run_dir not found: {run_dir}", file=sys.stderr)
        sys.exit(1)

    summary = build_summary(run_dir)
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(summary, f, indent=4)

    print(f"Summary written to {output_path}")


if __name__ == "__main__":
    main()
    