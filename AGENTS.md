# TrueMembers — Project Knowledge Map

## Purpose

TrueMembers automates the identification of all family members (homologs) of a given protein seed. It uses an iterative pipeline that combines JackHMMER-based sequence search, structural similarity filtering, and HMM domain validation, converging when few or no new family members are discovered.

## Architecture Overview

```
submit_run_true_memberes.sh    ← User runs `sbatch submit_run_true_memberes.sh` from repo root
    │  cd $SLURM_SUBMIT_DIR ensures job runs in project root
    │
    ▼
truemembers <config.yaml>      ← Console script (True_Members/cli.py:main), sets up dirs, starts JobScheduler
    │
    ▼
JobScheduler.py                ← Iteration loop: creates configs, submits Slurm jobs, collects results
    │
    ├──► JobResultCollector    ← Merges per-job results into global CSVs (seeds, used_seeds, all_members)
    │
    └──► (sbatch) → run_job.py → Job.py   ← A single job on a Slurm node
              │
              ├── evcouplings.py            ← Runs JackHMMER align stage with multiple bitscores
              ├── Foldseek webserver        ← Remote structure search (initial breadth)
              ├── get_structures_and_sequences.py  ← Downloads AlphaFold CIFs + UniProt FASTAs
              ├── structural_filter.py      ← Foldseek local + USalign → TM-score ≥ cutoff
              └── hmm_domain_filter.py      ← hmmsearch with domain HMMs → domain presence
```

**Important:** The `truemembers` command and all spawned Slurm jobs MUST be run from the project root directory. The CLI uses relative paths for template files (`True_Members/template_files/`), output directories (`runs/`), and config values like `hmms_folder_path` — all resolved relative to CWD.

## Key Classes

| Class | File | Role |
|-------|------|------|
| `JobScheduler` | `JobScheduler.py` | Master loop: monitors job status via `sacct`, spawns new jobs from seed pool, enforces max parallel jobs, handles convergence |
| `Job` | `Job.py` | Orchestrates a single pipeline run for one seed: evcouplings → structure download → structural filter → domain filter → combine results |
| `JobResultCollector` | `scripts/JobResultCollector.py` | Reads a finished job's CSV outputs, merges True Members into global seeds, All Members into global corpus, updates used_seeds with new_true_members count |
| `JobStatus` | `scripts/JobStatus.py` | Enum: PENDING, RUNNING, COMPLETED, HARVESTED, FAILED |
| `JobStats` | `scripts/types.py` | Dataclass tracking per-job pipeline stage counts (annotation files, IDs found, downloads, filter results, etc.). Replaced the ad-hoc `number_tracker` dict. |

## Pipeline Flow (per iteration)

### 1. Seed Selection
- Iteration 1: uses the seed from `run_true_members_config.yaml`
- Iteration 2+: picks the highest `qtmscore` seed from `seeds.csv` that has NOT been used yet

### 2. Single Job Run (`Job.run()`)

#### 2a. JackHMMER via evcouplings (`evcouplings.py`)
- Configures an evcouplings YAML with the seed's UniProt ID and residue region
- Creates one "batch" per bitscore value (e.g., `_b0.2`, `_b0.4`, …)
- Runs `evcouplings --yolo` to execute the `align` stage only (jackhmmer search)
- Output: per-bitscore annotation CSVs listing all homologs with sequence regions

#### 2b. Foldseek Webserver Search (`structural_filter.py:158-196`)
- Extracts the ROI structure of the **current job's seed**
- Posts a search job to the Foldseek webserver API
- Polls for completion, then downloads results
- Filters by E-value (labelled as TM-score equivalent) ≥ cutoff
- Output: `foldseek_webserver_output.csv`

#### 2c. Collect IDs & Download (`get_structures_and_sequences.py`)
- Collects unique (id, region_start, region_end) from JackHMMER annotation files + Foldseek webserver results
- Downloads AlphaFold CIF structures for all unique IDs via the EBI API
- Downloads UniProt/UniParc FASTA records for successfully downloaded structures
- Removes structures without a UniProt entry
- Output: `uids_with_region.fasta`, `uids_with_region.csv`, CIF files in `alphafold_structures/`

#### 2d. Structural Similarity Filter (`structural_filter.py`)
- **Foldseek local:** runs `foldseek easy-search` of seed ROI against all downloaded CIFs
- **USalign:** runs `USalign` for pairwise TM-score alignment of the seed ROI against all CIFs
- Merges Foldseek + USalign results, keeps sequences with `qtmscore ≥ TM_SCORE_CUTOFF` (default 0.6)
- Output: filtered DataFrame of structurally similar sequences

#### 2e. HMM Domain Filter (`hmm_domain_filter.py`)
- Presses all domain HMMs from the configured HMMs directory into one concatenated HMM
- Extracts GA (gathering threshold) cutoffs from HMM headers
- Runs `hmmsearch` on the FASTA of structurally filtered sequences
- Parses domtblout → pivot table: one boolean column per HMM, True if bitscore ≥ GA cutoff
- True Members = sequences where ALL HMM domains are present

#### 2f. Combine Results (`Job.py:combine_results`)
- Merges all sequences from JackHMMER/foldseek with structural filter results (left join)
- Drops any sequence that failed structural filter (no qtmscore)
- Merges with HMM domain filter results (left join, fills missing with False)
- True Members = those where all HMM domain columns are True
- Output: `results/All_Members.csv`, `results/True_Members.csv`, `results/number_tracker.csv`

### 3. Result Collection (`JobResultCollector`)
- Merges job's True Members into global `seeds.csv` (deduplicating by id, using median qtmscore, min/max region)
- Merges job's All Members into global `All_Members_global.csv`
- Updates `used_seeds.csv` with `new_true_members` count
- Failed jobs: their seed is removed from used_seeds so it can be retried

### 4. Convergence (optional)
- Convergence run mode (`convergence_run: True`): pipeline runs until N consecutive jobs each find fewer than `absolute_floor` new members
- Non-convergence mode: runs exactly `iterations` iterations

## Data Schemas

### `seeds.csv`
| id | qtmscore | region_start | region_end | {HMM_domain_1} | {HMM_domain_2} | … |
|----|----------|--------------|------------|----------------|----------------|---|

HMM domain columns are dynamically added based on the HMMs in the configured folder.

### `used_seeds.csv`
| id | region_start | region_end | iteration | new_true_members |
|----|--------------|------------|-----------|-----------------|

### `All_Members_global.csv`
Same schema as `seeds.csv`, but contains ALL sequences that passed structural + domain filters (not just True Members).

## Configuration

### Master Config (`run_true_members_config.yaml` — at repo root)
- `global.iterations` — max iterations in non-convergence mode
- `global.max_parallel_jobs` — Slurm job concurrency limit
- `global.convergence_run` — enable convergence mode
- `global.absolute_floor` — threshold for "new members per job" convergence
- `global.patience_window` — consecutive jobs below floor to declare convergence
- `run.seed_uid` — initial UniProt accession
- `run.seed_region_start/end` — residue region of interest
- `run.bit_scores` — list of JackHMMER bitscore thresholds
- `run.tm_score_cutoff` — TM-score cutoff for structural filter
- `run.structural_filter_tool` — which structural filter to use: `usalign`, `foldseek`, or `both`
- `run.hmms_folder_path` — directory containing `.hmm` files (path relative to project root)
- `run.CPU` — cores per job

### Per-Job Config (`job_config_template.yaml`)
Generated dynamically by `JobScheduler.create_job_config()`. Populated with current seed, iteration number, and a unique `root_dir`.

## External Dependencies

- **jackhmmer** (via evcouplings 0.2.1) — iterative sequence search
- **Foldseek** — structural alignment (local binary + webserver API)
- **USalign** — pairwise TM-score computation
- **hmmsearch/hmmpress** (HMMER) — domain detection
- **Slurm** — job scheduling on HPC cluster
- **UniProt REST API** — sequence downloads
- **AlphaFold DB API** — structure downloads

## Code Conventions

- **Imports:** Use explicit imports (no wildcard `from x import *`). Group imports: stdlib, third-party, local.
- **Logging:** Use `logging.getLogger(__name__)` module-level loggers. The root logger is configured by `Job._setup_logging()` (per-job: `job.log`) and `JobScheduler._setup_logging()` (per-pipeline: `pipeline.log`). Use `logger.info()` for progress, `logger.warning()` for recoverable issues, `logger.error()` for failures.
- **Stats:** Use `JobStats` dataclass (`scripts/types.py`) instead of ad-hoc dicts. Pass by reference to pipeline functions; each function sets the fields it's responsible for.
- **Config:** All pipeline configuration flows through YAML files. Master config (`run_true_members_config.yaml`) → JobScheduler → per-job config (dynamically generated from `job_config_template.yaml`) → Job.

## Testing

- **Unit tests** (`tests/test_*.py`): Pure Python logic, no external deps or network. Test string parsing, DataFrame operations, HMM parsing, convergence math.
- **Integration tests** (`tests/test_integration_*.py`): Mock external tool calls (`subprocess.run`, `requests`) and use real annotation CSVs from `test_data/`.
- **E2E:** Manual HPC validation before releases.
- Run: `pytest tests/ -v`
- Coverage: `pytest tests/ -v --cov=True_Members --cov-report=term-missing`
- See `testing_strategy.md` for detailed test plan.

## Repo File Layout

```
TrueMembers/                       ← Project root
├── True_Members/                  ← Python package (installed via pip)
│   ├── cli.py                     ← Console script entry point (truemembers)
│   ├── Job.py
│   ├── JobScheduler.py
│   ├── run_job.py
│   ├── scripts/
│   │   ├── JobResultCollector.py
│   │   ├── JobStatus.py
│   │   ├── types.py
│   │   └── utils.py
│   ├── template_files/
│   │   ├── job_config_template.yaml
│   │   ├── evcouplings_config_template.yaml
│   │   └── submit_true_members_job_template.sh
│   ├── test_data/
│   └── tests/
├── hmm_domain_filters/            ← Domain HMM profiles (.hmm files)
├── run_true_members_config.yaml   ← Master pipeline config
├── submit_run_true_memberes.sh    ← Slurm submission script (user runs: sbatch submit_run_true_memberes.sh)
├── hpc_environment.yml            ← Conda environment for HPC
├── pyproject.toml                 ← Package metadata + entry point
├── README.md
├── AGENTS.md
└── runs/                          ← Pipeline output (created at runtime)
```

## Directory Layout (per pipeline run)

```
runs/run_<seed>_<timestamp>/
├── data/
│   ├── seeds.csv
│   ├── used_seeds.csv
│   └── All_Members_global.csv
├── shared/
│   ├── structures/               ← All AlphaFold CIFs (cached across jobs)
│   └── sequences/                ← All UniProt FASTA records (cached across jobs)
├── job_configs/
├── slurm_scripts/
├── final_true_members_results/
├── pipeline.log                  ← Scheduler log
├── runs/
│   └── job_<seed>_iteration_<N>_<timestamp>/
│       ├── evcouplings/output/<seed>_<region>_b<score>/align/
│       ├── filtering/foldseek/
│       ├── filtering/usalign/
│       ├── filtering/hmm_domain_filtering/
│       ├── results/
│       └── job.log               ← Per-job log
└── <initial_seed>.cif_roi.cif    ← Initial seed ROI (downloaded to shared dir)
```

## Shared Directory Architecture

To avoid redundant downloads across iterations, each pipeline run creates `shared/structures/` and `shared/sequences/` directories. All jobs in the same pipeline run read from and write to these shared directories:

- **Structures:** `download_alphafold_structure()` checks `os.path.exists()` before downloading, so subsequent jobs skip files already fetched. Atomic writes (`.tmp` → rename) prevent corruption on shared filesystems.
- **Sequences:** Records are downloaded from UniProt on demand per job (they are small text files and the API is rate-limited, so per-ID caching is less critical).
- **Thread-safety:** `os.rename()` is atomic on NFS and Lustre. If two jobs download the same structure simultaneously, the last atomic rename wins — both produce valid results.