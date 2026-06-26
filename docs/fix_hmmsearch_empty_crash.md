# Fix: hmmsearch Crash When No Structurally Filtered Sequences Exist

## Problem

In `True_Members/scripts/hmm_domain_filter.py:run_hmm_domain_filter()` (line 15-27), when the structural filter finds **zero** sequences (i.e., no homolog passes the TM-score cutoff), the pipeline crashes. The crash occurs because:

1. `seed_containing_sequences` is an empty DataFrame
2. `seed_containing_records` (line 19) becomes an empty list
3. `write_fasta_from_records` (line 21) writes an **empty FASTA file**
4. `hmmsearch` on an empty FASTA either returns a non-zero exit code (crashing with `RuntimeError` at line 125) or produces a domtblout without valid hits that `SearchIO.parse` cannot handle

## Root Cause Location

- **`hmm_domain_filter.py:15-27`** — `run_hmm_domain_filter()` unconditionally runs hmmsearch even when `seed_containing_sequences` is empty
- **`hmm_domain_filter.py:105-125`** — `run_hmmsearch()` raises RuntimeError on non-zero exit from hmmsearch (which happens on empty input)

## Fix Strategy

Add an early-return guard in `run_hmm_domain_filter()`: if `seed_containing_sequences` is empty, log a warning, set `stats.hmm_search_results = 0`, and return an empty DataFrame without invoking hmmsearch.

## Implementation Steps

### 1. Modify `run_hmm_domain_filter()` in `True_Members/scripts/hmm_domain_filter.py`

Add an early-return check after line 19 (after filtering records by seed_containing_sequences):

```python
def run_hmm_domain_filter(stats: JobStats, cpus: int, uids_with_region_fasta_path: str, seed_containing_sequences: pd.DataFrame, HMMs_DIR: str, hmm_domain_filtering_dir: str, root_dir: str, TM_SCORE_CUTOFF: float, hmms_with_ga_cutoffs: dict) -> pd.DataFrame:
    domain_validator_hmm_path = press_hmms(hmm_domain_filtering_dir, HMMs_DIR)

    records = get_records_from_fasta(uids_with_region_fasta_path)
    seed_containing_records = [record for record in records if extract_uid_from_recordID(record.id) in seed_containing_sequences["id"].tolist()]

    if not seed_containing_records:
        logger.warning("No structurally filtered sequences to run hmmsearch on. Skipping HMM domain filtering.")
        stats.hmm_search_results = 0
        return pd.DataFrame()

    seed_containing_sequences_full_fasta = f"{root_dir}/filtering/seed_containing_sequences_tm_{TM_SCORE_CUTOFF}.fasta"
    # ... rest of existing code unchanged
```

Indent the remaining lines (21-27) under the `if` branch. The function returns an empty DataFrame when there are no sequences, which `combine_results()` in `Job.py` handles gracefully (it already merges with the result and fills missing HMM columns with `False`).

### 2. Verify `parse_hmmsearch_results()` handles empty domtblout safely (already correct)

`parse_hmmsearch_results()` at line 43-58 correctly initializes `hmm_search_results = pd.DataFrame()` and only populates it if `data` is non-empty. With the early return added above, this code path won't be reached for empty inputs, but it serves as a secondary safety net.

### 3. Verify `combine_results()` in `Job.py` handles empty HMM results

`Job.py:129-143` — `combine_results()` does a left merge with `hmm_search_results` (line 132), then fills missing HMM columns with `False` (lines 134-136). This works correctly with an empty DataFrame returned from `run_hmm_domain_filter()`. When all HMM columns are False, `true_members` will be empty (line 141), which is the expected behavior when no sequences pass structural filtering.

## Files Changed

- `True_Members/scripts/hmm_domain_filter.py` — Add early return to `run_hmm_domain_filter()`

## Testing

- Unit test: mock `seed_containing_sequences` as empty DataFrame, verify `run_hmm_domain_filter()` returns empty DataFrame and sets `stats.hmm_search_results = 0` without calling `subprocess.run`
- Integration test: run a job with a TM-score cutoff of 1.0 (guarantees no structural matches), verify pipeline completes without crashing and produces empty True_Members.csv