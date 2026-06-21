# TrueMembers — Testing Strategy

## Current State

- **No automated tests exist.** The only validation is a Jupyter notebook (`tests.ipynb`) that manually runs `Job.run()` with pre-computed evcouplings output and `run_evcouplings()` commented out.
- **Test data available:** Two evcouplings annotation CSVs (JackHMMER output) — one large (~5869 entries, bitscore 0.6), one small (~19 entries, bitscore 0.8, UniRef90-clustered).
- **No structures or sequences in test data:** CIFs and FASTAs need to be downloaded from AlphaFold/UniProt APIs at test time.
- **External tools not available locally:** evcouplings, foldseek, USalign, hmmsearch, hmmpress — all only on HPC via modules.

## What Can We Test? (Layered Approach)

```
Layer 3: End-to-end (HPC)     ← manual, one-off, full pipeline on real cluster
Layer 2: Integration (local)  ← mock external tools, use real annotation CSVs
Layer 1: Unit tests (local)   ← pure Python logic, no external deps
```

### Layer 1: Unit Tests — Pure Python Logic

These run anywhere (no HPC, no external tools, no network). They test individual functions with deterministic inputs.

| Test | Function(s) | What it validates | Test data needed |
|------|------------|-------------------|-----------------|
| **ID parsing** | `extract_uid_from_recordID()` | Parses `tr\|A0A1A9EVS2\|A0A1A9EVS2_9GAMM/553-621` → `A0A1A9EVS2`. Handles edge cases: UPI prefixes, bare accessions, malformed IDs. | None (inline strings) |
| **Region parsing** | `extract_region_from_recordID()` | Parses `/553-621` → `(553, 621)`. Handles missing regions. | None |
| **Annotation file processing** | `get_unique_ids_with_regions_from_annotation_files()` | Given real and mock annotation CSVs, extracts unique IDs with min/max regions. Dedup logic (same ID at multiple bitscores → wider region). | Real `*_annotation.csv` from test_data |
| **JackHMMER + Foldseek ID merging** | `combine_foldseek_and_jackhmmer_ids()` | Merges two DataFrames on `id`, averages regions. Handles empty inputs, non-overlapping sets. | Synthetic DataFrames |
| **Structure filtering** | `filter_foldseek_output()`, `filter_usalign_output()` | Parses foldseek TSV and USalign output, filters by TM-score. Handles empty results, edge-case column names. | Small fake TSV files or inline strings |
| **HMM domain parsing** | `get_hmms_with_ga_cutoffs()`, `parse_name()`, `parse_ga()` | Correctly extracts NAME and GA cutoffs from HMM files. Handles HMMs without GA lines, empty directories. | HMM files in `hmm_domain_filters/` (already exist) |
| **HMM results parsing** | `parse_hmmsearch_results()` | Parses a domtblout file into a boolean pivot table. Handles: no hits, all-positive hits, mixed hits, HMMs not present in output. | Small fake domtblout strings |
| **Result combination** | `Job.combine_results()` | Left-join logic: sequences without structural match are dropped, missing HMM columns are filled with False. All-HMM-domain check. | Synthetic DataFrames |
| **Convergence computation** | `JobScheduler.compute_convergence()` | Boundary cases: fewer jobs than patience window, all below floor, exactly at floor, one above floor, empty used_seeds. | Synthetic used_seeds DataFrames |
| **Seed selection** | `JobScheduler.get_highest_scoring_seed()` | Picks highest qtmscore unused seed. Handles: all used, single seed, tie scores, empty seeds. | Synthetic seeds + used_seeds DataFrames |
| **Result merging** | `JobResultCollector.merge_local_true_members_with_global_seeds()` | Dedup by id, median qtmscore, min region_start, max region_end. Computes new_true_members count. Handles: no overlap, full overlap, partial overlap. | Synthetic DataFrames |
| **All members merging** | `JobResultCollector.merge_local_all_members_with_global_all_members()` | Aggregation: any=True for HMM columns. Handles new HMM columns appearing mid-pipeline. | Synthetic DataFrames |
| **Job config creation** | `JobScheduler.create_job_config()` | Template YAML populated correctly for iteration 1 vs iteration N. All expected keys present. | Template YAML files (already exist) |
| **Slurm script generation** | `JobScheduler.create_slurm_script()` | Template placeholders replaced correctly. Output script is valid bash. | Template shell script (already exists) |

**Total:** ~13 unit test groups. All can run in <1 second each. No external dependencies beyond pandas and pytest.

### Layer 2: Integration Tests — Mocked Pipeline Run

Test the pipeline end-to-end locally by mocking only external tool calls, using real test annotation CSVs.

**Approach:** Use `unittest.mock.patch` to replace:
- `subprocess.run` (evcouplings, foldseek, USalign, hmmsearch, hmmpress, sbatch, sacct)
- `requests.get/post` (AlphaFold API, UniProt API, Foldseek webserver)
- File I/O for external tool outputs (provide pre-computed result files)

**Integration test scenarios:**

| Scenario | What it tests |
|----------|--------------|
| **Happy path: single job** | `Job.run()` completes all steps. Uses real annotation CSVs, mocks all external tools. Verifies output CSVs exist with correct schema. |
| **Happy path: multi-iteration scheduler** | `JobScheduler.start()` with mocked `sbatch` and `sacct`. Jobs complete instantly. Verifies seeds are selected correctly, convergence works. |
| **Failed job recovery** | Mock `sacct` reporting FAILED. Verifies seed is freed from used_seeds and can be retried. |
| **Convergence stops pipeline** | Mock `new_true_members` decreasing below floor. Verifies pipeline stops at correct iteration. |
| **Max iteration stop** | Non-convergence mode. Verifies pipeline stops at `MAX_ITERATIONS`. |
| **All seeds used** | Mock seed pool exhausted. Verifies pipeline waits for running jobs then terminates. |
| **Missing annotation file** | One bitscore produces no annotation CSV. Verifies warning logged, pipeline continues. |
| **Empty JackHMMER results** | All annotation files empty. Verifies pipeline fails gracefully (no crash). |

**Key implementation detail:** The `use_existing_root_dir` parameter on `Job.__init__` already supports reusing pre-existing directories. For integration tests, set up a test directory with:
1. Pre-populated `evcouplings/output/` (copy from `test_data/`)
2. Mocked `alphafold_structures/` (provide a few small CIF files)
3. Mocked foldseek output files
4. Mocked USalign output files
5. Mocked hmmsearch domtblout output

### Layer 3: End-to-End — HPC Validation

Manual process, run before releases or major changes:
1. Submit a single `run_job.py` with a known seed on the cluster
2. Verify all pipeline stages produce expected outputs
3. Check that `True_Members.csv` contains the known family members
4. Run a 2-iteration pipeline and verify convergence behavior

## Test Framework Choice: pytest

- `pytest` is the standard Python testing framework
- `pytest-mock` provides clean mocking (thin wrapper over `unittest.mock`)
- `pytest-cov` for coverage reporting
- Test files in `tests/` directory (not `test_data/`)

## Proposed File Structure

```
True_Members/
├── tests/
│   ├── __init__.py
│   ├── conftest.py                    ← shared fixtures (test data paths, mock helpers)
│   ├── test_utils.py                  ← unit tests for utils.py
│   ├── test_structural_filter.py      ← unit tests for filter_foldseek_output, filter_usalign_output
│   ├── test_hmm_domain_filter.py      ← unit tests for parse_hmmsearch_results, get_hmms_with_ga_cutoffs
│   ├── test_get_structures.py         ← unit tests for combine_foldseek_and_jackhmmer_ids, get_annotation_files
│   ├── test_job.py                    ← unit tests for Job.combine_results
│   ├── test_job_scheduler.py          ← unit tests for compute_convergence, get_highest_scoring_seed, create_job_config
│   ├── test_job_result_collector.py   ← unit tests for merge logic
│   ├── test_integration_job.py        ← integration: single Job.run() with mocks
│   └── test_integration_scheduler.py  ← integration: JobScheduler.start() with mocks
├── test_data/                         ← existing, unchanged (annotation CSVs used by tests)
├── test_data/mock_outputs/            ← pre-computed mock outputs for integration tests
│   ├── foldseek_output.csv
│   ├── usalign_results.txt
│   └── hmmsearch_results.domtblout
└── ...
```

## Shared Fixtures (`conftest.py`)

```python
import pytest
from pathlib import Path

TEST_DATA = Path(__file__).parent.parent / "test_data"

@pytest.fixture
def annotation_csv_small():
    """The small (19-entry) annotation CSV."""
    return str(TEST_DATA / "output_small" / "Q5SK37_507-569_b0.8" / "align" / "Q5SK37_507-569_b0.8_annotation.csv")

@pytest.fixture
def annotation_csv_large():
    """The large (~5869-entry) annotation CSV."""
    return str(TEST_DATA / "output" / "Q5SK37_507-569_b0.6" / "align" / "Q5SK37_507-569_b0.6_annotation.csv")

@pytest.fixture
def hmms_dir():
    return str(Path(__file__).parent.parent / "hmm_domain_filters")

@pytest.fixture
def mock_subprocess_run(mocker):
    """Mock subprocess.run for external tool calls."""
    return mocker.patch("subprocess.run")

@pytest.fixture
def mock_requests(mocker):
    """Mock requests.get/post for API calls."""
    mock_get = mocker.patch("requests.get")
    mock_post = mocker.patch("requests.post")
    return mock_get, mock_post
```

## Implementation Priority (Build Order)

| Order | Test file | Rationale |
|-------|-----------|-----------|
| **1** | `test_utils.py` | Simplest unit tests, no mocking needed (just string parsing). Validates the test harness. |
| **2** | `test_hmm_domain_filter.py` | Uses existing HMM files. Tests pure parsing logic. |
| **3** | `test_get_structures.py` | Uses real annotation CSVs from test_data. Tests the core ID collection logic. |
| **4** | `test_structural_filter.py` | Needs small fake TSV strings (easy to create inline). |
| **5** | `test_job.py` | Tests `combine_results` with synthetic DataFrames. |
| **6** | `test_job_scheduler.py` | Tests convergence and seed selection — pure logic. |
| **7** | `test_job_result_collector.py` | Tests merge logic — pure pandas operations. |
| **8** | `test_integration_job.py` | End-to-end Job run with mocks. Depends on all previous tests passing. |
| **9** | `test_integration_scheduler.py` | End-to-end scheduler with mocks. Most complex. |

## Running Tests

```bash
# All tests
pytest tests/ -v

# Unit tests only (fast)
pytest tests/ -v --ignore=tests/test_integration_job.py --ignore=tests/test_integration_scheduler.py

# With coverage
pytest tests/ -v --cov=True_Members --cov-report=term-missing

# Single test file
pytest tests/test_utils.py -v
```

## What NOT to Test

- **evcouplings subprocess itself** — tested by the evcouplings maintainers. We only test our config generation.
- **foldseek/USalign/hmmsearch correctness** — tested by their maintainers. We test our output parsing.
- **Network downloads** — tested by AlphaFold/UniProt APIs. We mock these.
- **Slurm job scheduling** — tested by the Slurm system. We only mock `sbatch`/`sacct`.

## Risks & Limitations

1. **Mock fidelity:** Integration tests mock external tools. If a tool changes its output format, mocks won't catch it. Mitigation: run Layer 3 (E2E on HPC) before releases.
2. **Test data freshness:** Annotation CSVs are snapshots of JackHMMER output. If evcouplings/UniProt DB changes, results differ. Mitigation: unit tests don't depend on specific annotation content — only on format and edge cases.
3. **No structure/sequence files in test data:** Tests that need CIF or FASTA files must mock the download functions. Mitigation: provide a small (1-2) test CIF files in `test_data/` so structural filter parsing can be tested with real data.