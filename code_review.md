# TrueMembers — Code Review

## Quick Summary

The codebase has a sound architectural design (scheduler/worker pattern for HPC, clean pipeline steps) and good HTTP client practices. The main weaknesses are: two confirmed bugs, lack of tests, mutable global state via dict, print-based logging, race conditions on CSV files, and inconsistent type usage. Overall the code is functional and the design is solid — these issues are fixable with targeted refactoring.

---

## What's Good

### Architecture
- **Clean separation of concerns:** `JobScheduler` handles orchestration, `Job` handles execution, `JobResultCollector` handles merging — each class has a clear, narrow responsibility.
- **Template-based configuration:** Using YAML templates for both evcouplings and job configs is clean and avoids hardcoding.
- **Soft failure handling:** Failed jobs' seeds are automatically freed from `used_seeds` so they can be retried. This is robust for HPC environments where transient failures (node crash, OOM) are common.
- **Convergence mode:** Well-designed with `absolute_floor` + `patience_window` — avoids stopping on a single noisy iteration.

### HTTP & Networking
- **Connection pooling via `_build_session()`:** Using a single `requests.Session` with `HTTPAdapter` and `Retry` for backoff is professional-grade. Most codebases miss this.
- **ThreadPoolExecutor for parallel downloads:** Appropriate for I/O-bound tasks.
- **Progress logging with adaptive intervals:** `_log_interval()` scales from 1-item batches to larger ones — good UX for both small and large jobs.

### Data Processing
- **HMM domain filter pivot table:** The `parse_hmmsearch_results()` function that parses domtblout into a boolean pivot table is elegant and efficient.
- **GA cutoff extraction from HMM files:** Proper use of HMMER's gathering thresholds (`parse_ga()`) instead of arbitrary E-value cutoffs.
- **Region extraction from record IDs:** The `/<start>-<end>` parsing in `extract_region_from_recordID()` handles the evcouplings naming convention correctly.

### Robustness
- **`os._exit(0)` in `run_job.py`:** Explicitly handles the multiprocessing pool cleanup issue — pragmatic and effective.
- **Graceful handling of missing annotation files:** `get_annotation_files()` warns but doesn't crash if a bitscore batch produced no output.

---

## Bugs Confirmed

### 1. Duplicate evcouplings call (`evcouplings.py:9-11` and `:14-23`)
Lines 9-11 call `subprocess.run` without `capture_output`, `check`, or error handling — this runs evcouplings and dumps output to stdout with no error checking. Lines 14-23 then re-run it with proper error handling. The first call is dead code that doubles execution time.

**Fix:** Delete lines 9-11.

### 2. Foldseek webserver uses wrong ROI (`Job.py:86` / `structural_filter.py` via caller)
In `Job.py:run()`, the Foldseek webserver search downloads the **current seed's** structure and extracts its ROI, but the RHS of line 86 (`run_foldseek_webserver(...)`) then passes this to the webserver. The AC also confirmed: it should use `self.SEED`'s ROI (the job's current seed), not `self.GLOBAL_SEED_structure_roi` (the initial pipeline seed). Wait — re-reading the code, the webserver call on line 88 of Job.py uses `query_structure_roi` which is extracted from the current seed's downloaded structure (lines 85-86). But line 93 (`run_structural_similarity_filter`) passes `self.GLOBAL_SEED_structure_roi` — the initial seed's ROI. This means:
- **The webserver search (for breadth) correctly uses the current seed's ROI** — this is fine.
- **The structural filter (Foldseek local + USalign) uses the initial seed's ROI** — this means ALL jobs compare against the same reference structure regardless of which seed they're processing.

This is a **semantic inconsistency**. The structural filter should compare candidate structures against the current job's seed, not the pipeline's initial seed. Otherwise, a seed with qtmscore=0.6 against the initial seed gets added to the pool, but then job for THAT seed also compares everything against the initial seed — so the structural filter never re-anchors on the new seed.

**Fix:** Extract the ROI from the current seed's structure and pass THAT to `run_structural_similarity_filter`. The `GLOBAL_SEED_structure_roi` config field should be renamed to make clear what it represents.

### 3. Race condition on global CSV files (`JobResultCollector.py`)
`get_results_from_finished_jobs()` reads CSV files with `pd.read_csv()`, modifies them, and writes them back. Between the read and the write, another finished job's results could be merged in (the scheduler calls this once per loop iteration after updating all job statuses — this handles batch updates, but concurrent collection isn't the concern here). The real issue: if the scheduler process crashes mid-write, CSV corruption can occur. The fix is atomic writes: write to a `.tmp` file, then rename.

### 4. Potential KeyError on missing annotation files
`get_jackhmmer_and_foldseek_ids_with_regions()` at line 45-46 counts annotation files, but if `annotation_files` is empty, `get_unique_ids_with_regions_from_annotation_files()` returns an empty DataFrame and the code continues — but then `run_download_tasks([], ...)` runs with an empty list, and later operations on `combined_unique_ids_with_regions` may be semantically wrong (though they won't crash). The pipeline should detect this case and fail the job gracefully.

---

## Design Issues

### 1. `number_tracker` dict as mutable global state
The `number_tracker` dict is passed by reference through 5+ function calls (`Job.run() → run_evcouplings → run_foldseek_webserver → get_all_structures_and_sequences → get_jackhmmer_and_foldseek_ids_with_regions → etc.`). Each function mutates it via `number_tracker["some_key"] = value`. This makes it extremely hard to trace state changes and debug.
**Fix:** Make each function return its own stats dict, then merge them at the top level. Or use a dataclass.

```python
@dataclass
class JobStats:
    annotation_files: int = 0
    ids_jackhmmer: int = 0
    ids_foldseek_webserver: int = 0
    ids_combined: int = 0
    structure_downloads: int = 0
    uniprot_records: int = 0
    hmm_search_results: int = 0
    all_members: int = 0
    true_members: int = 0
```

### 2. Print-based logging
All status messages use `print()`. On Slurm, `print()` goes to the job's stdout file, which is fine but doesn't provide:
- Timestamps for each event
- Severity levels (INFO vs WARNING vs ERROR)
- Structured output for downstream parsing

**Fix:** Replace `print()` with Python's `logging` module. Slurm jobs should log to both console (stdout) and a file.

```python
import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(), logging.FileHandler(f"{root_dir}/job.log")]
)
logger = logging.getLogger(__name__)
```

### 3. No tests
There are zero test files. The test_data directory exists with some sample annotation CSVs, but no test runner or test functions. Given the pipeline's complexity (7+ steps per job, merging logic, convergence computation), bugs in edge cases are highly likely.

**Fix:** Add:
- Unit tests for `combine_results()` (merge logic)
- Unit tests for `parse_hmmsearch_results()` (domtblout parsing)
- Unit tests for `get_unique_ids_with_regions_from_annotation_files()`
- Unit tests for `compute_convergence()` (boundary cases)
- Integration test: run a job with the test data (already has annotation CSVs)
- Mock external tools (evcouplings, foldseek, hmmsearch) for fast unit tests

### 4. CSV files as database
Using CSV files for global state (`seeds.csv`, etc.) works for small datasets but has problems:
- No atomicity (crash during write = corruption)
- No concurrent access safety
- Schema is implicit (HMM columns appear dynamically)
- Full file read/write for incremental updates

**Fix:** Consider SQLite with WAL mode. The migration is straightforward since all operations are pandas `read_csv` / `to_csv` — replace with `pd.read_sql` / `df.to_sql`.

### 5. Snake-case vs camel-case inconsistency
- Class names: PascalCase (`JobScheduler`, `JobResultCollector`) — correct
- Function names: snake_case (`run_evcouplings`, `download_file`) — correct
- Variable names: UPPER_CASE for pseudo-constants (`HMMs_DIR`, `TM_SCORE_CUTOFF`) — inconsistent with Python conventions. These are constructor parameters, not module-level constants. Use `self.hmm_dir` instead of `self.HMMs_DIR`.

### 6. Hardcoded paths in template files
`submit_true_members_job_template.sh` contains hardcoded paths:
```bash
cd /ptmp/arowi/STAC_CbrA_Aggregation/true_members_pipeline_runs
```
This makes the pipeline non-portable. The `cd` should be a template variable like `<pipeline_root>` that gets substituted by the scheduler.

### 7. `evcouplings.py` first subprocess call has no error handling
Even after removing the duplicate (bug #1), the remaining call on lines 14-23 uses `capture_output=True`, which means errors are captured but also means no progress output is visible during the long jackhmmer run (minutes to hours). Consider using `subprocess.PIPE` for stderr only and letting stdout stream to the terminal, or use Popen with polling.

### 8. `get_records_from_accessions` filters with regex
```python
if record.id not in ids:
```
This is linear search — O(n²). For large result sets (hundreds of sequences), this can be slow. Convert `ids` to a set.

### 9. `parse_ga()` parsing is fragile
```python
parts = clean_line.rstrip(';').split()
if len(parts) >= 3:
    domain_cutoff = parts[2].rstrip(';')
```
The `.rstrip(';')` on the whole line is redundant (split already handles whitespace), and the second `.rstrip(';')` on `parts[2]` is not needed since there's no `;` at position 2 in a GA line. The GA line format is `GA  25.0  25.0;` — fields are space-separated, not semicolon-delimited. The semicolon is only a line terminator at position 4. This works by accident because `split()` splits on whitespace.

**Fix:** Cleaner parsing with regex.

---

## Style & Minor Issues

### `Job.py`
| Location | Issue | Fix |
|----------|-------|-----|
| L3-L9 | Wildcard imports (`from scripts.utils import *`) | Import only what's used |
| L7 | Duplicate import of evcouplings | Remove one |
| L84 | Comment refers to `run_download_tasks` but actually downloads and extracts seed structure | Split into clearly named helper method |
| L86 | Variable `query_structure_roi` name is misleading — it's the current seed's ROI, not the initial seed's | Rename to `current_seed_roi` |

### `JobScheduler.py`
| Location | Issue | Fix |
|----------|-------|-----|
| L57 | `self.jobs: dict[int, dict]` — nested dict type hint | Use TypedDict or NamedTuple for inner dict |
| L142-143 | `time.sleep(180)` — magic number | Extract to config or constant |
| L228 | `top_seed_row["id"].iloc[0]` — fragile indexing | Use `.iat[0, 0]` or check for empty |
| L257 | `new_seed` reused across two branches | Separate variable names for clarity |

### `structural_filter.py`
| Location | Issue | Fix |
|----------|-------|-----|
| L18 | `filter_foldseek_output`: unused parameter `evalue` is selected but `qtmscore` comparison uses integer column index | Hardcoded column order is fragile if foldseek output format changes — use named column access |
| L55 | `filter_usalign_output`: renames `TM2` to `qtmscore` after reading | Add comment explaining why TM2 (not TM1) is used |
| L158-196 | `run_foldseek_webserver` query structure is uploaded but annotation CSV is loaded from disk | Race condition if evcouplings hasn't finished writing |

### `utils.py`
| Location | Issue | Fix |
|----------|-------|-----|
| L14 | Unused import `from Bio.SeqRecord import SeqRecord` | Remove |
| L80 | `get_records_from_fasta` could use generator to avoid loading all into memory | `SeqIO.parse` returns iterator — `list()` forces full load |

### `hmm_domain_filter.py`
| Location | Issue | Fix |
|----------|-------|-----|
| L5 | `import sys` is unused | Remove |
| L128-131 | `press_hmms` concatenates HMM files by bytes | If HMM files have trailing newline issues, hmmpress may fail. Use text concatenation with explicit newlines. |

---

## Recommended Priority Fix List (Implementation Order)

| Priority | Fix | File(s) | Effort |
|----------|-----|---------|--------|
| **P0** | Remove duplicate evcouplings call (Bug #1) | `evcouplings.py:9-11` | 1 min |
| **P0** | Fix structural filter to use current seed ROI (Bug #2) | `Job.py:93`, `structural_filter.py` | 30 min |
| **P1** | Replace `number_tracker` dict with dataclass (Design #1) | `Job.py`, all scripts | 1 hr |
| **P1** | Replace print() with logging (Design #2) | All files | 30 min |
| **P1** | Atomic CSV writes (Bug #3) | `JobResultCollector.py` | 15 min |
| **P2** | Add unit tests for merge/filter/convergence logic (Design #3) | New test files | 3 hr |
| **P2** | Fix variable naming conventions (Style #5) | `Job.py`, `JobScheduler.py` | 30 min |
| **P2** | De-hardcode paths in Slurm template (Design #6) | `submit_true_members_job_template.sh` | 10 min |
| **P3** | Migrate global state to SQLite (Design #4) | `JobResultCollector.py`, `JobScheduler.py` | 3 hr |
| **P3** | Fix import style (wildcard imports → explicit) | All files | 15 min |