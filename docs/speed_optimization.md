# TrueMembers — Speed Optimization Analysis

## Primary Bottlenecks (Ranked by Impact)

### 1. Redundant AlphaFold Structure Downloads (HIGH impact)

**Problem:** Each job downloads AlphaFold CIF structures independently into its own `alphafold_structures/` directory. In practice, different jobs (iterations) discover overlapping sets of sequence IDs. The same structures are downloaded repeatedly.

**Example:** If job 1 discovers 500 IDs and job 2 discovers 400 IDs with 300 overlap, 300 structures are downloaded twice.

**Solution: Per-pipeline-run shared `structures/` directory**

Create a shared directory inside the pipeline run root:

```
runs/run_<seed>_<timestamp>/
├── shared/
│   ├── structures/          ← All downloaded AlphaFold CIFs
│   └── sequences/           ← All downloaded UniProt FASTA records
├── data/
├── runs/
└── ...
```

Implementation changes:

- Add `shared_structures_dir` and `shared_sequences_dir` to the paths dict in `run_true_members.py:set_up_dirs()`
- Pass these paths to `JobScheduler`, which passes them into `create_job_config()`
- In `get_structures_and_sequences.py`, modify `run_download_tasks()` to check the shared directory first
- `download_alphafold_structure()` already has a `os.path.exists(save_path)` check — it just needs to point at the shared dir
- Same for FASTA records: write to shared `sequences/` dir, check before downloading
- **Thread-safety:** Since Slurm jobs run on separate nodes and the shared directory is on a shared filesystem (NFS/Lustre), add a per-file lock or use atomic writes (write to temp file, then `os.rename()` which is atomic on most filesystems)

**Estimated speedup:** For a typical 5-iteration run with 50% ID overlap across jobs: ~40-60% reduction in download time. If a job spends 30 min downloading structures, this saves ~12-18 min per job.

### 2. Redundant UniProt Sequence Downloads (HIGH impact)

**Problem:** Same as above — each job independently downloads FASTA records from UniProt for every ID it discovers.

**Solution:** Shared `sequences/` directory with cached FASTA records.

- Maintain a combined `all_records.fasta` in the shared directory
- Before downloading, check if ID is already present
- Append new records to the combined FASTA, write per-job filtered FASTAs on demand
- Alternatively: one FASTA file per accession (more filesystem-friendly for parallel access)

**Estimated speedup:** Similar to structure downloads — 40-60% reduction in sequence download time.

### 3. Foldseek Webserver Search Runs on Every Job (MEDIUM-HIGH impact) DO NOT IMPLEMENT since the foldseek webserversearch does not take that long and should search with the current seed's ROI, not the initial seed's ROI

**Problem:** `Job.run()` always runs `run_foldseek_webserver()`, which posts a job to the remote Foldseek API and polls for completion. This is an inherently slow network operation (can take minutes) and is currently bugged — it uses the initial seed's ROI instead of the current job's seed ROI.

**Solution A (fix + optimize):** As you noted, it should use the current seed's ROI. Fixed, this webserver call provides valuable breadth by searching huge databases (BFVD, afdb50, etc.). However, it can still be expensive.

**Solution B (cache per seed):** Cache webserver results per seed ID in the shared directory. If two different iterations happen to use the same seed (which shouldn't happen, since used_seeds blocks re-use), skip the call.

**Solution C (optionalize):** Make the webserver search configurable. For some families, the JackHMMER results alone may be sufficient, and the structural filter from local Foldseek is enough.

### 4. JackHMMER Runs per Bitscore Independently (MEDIUM impact) DO NOT IMPLEMENT the different biscore runs are the core of the sequence search approach

**Problem:** `run_evcouplings()` creates N batches for N bitscores. Each batch runs jackhmmer independently against the same database. This means N independent searches against UniProt for the same seed.

**Analysis:** This is inherent to the approach — each bitscore gives a different sensitivity/specificity tradeoff, and the union captures a broader range of homologs. However, jackhmmer itself is iterative and already produces an MSA — re-running from scratch at different bitscores discards the intermediate alignments.

**Potential optimization:** Modify evcouplings to run a single jackhmmer search at the lowest bitscore, then post-filter the resulting alignment at each bitscore threshold. This requires understanding evcouplings internals — the `reuse_alignment: True` setting in the template already helps, but running jackhmmer from scratch at each bitscore still rebuilds the HMM iteratively.

**Estimated speedup:** If implemented, could reduce jackhmmer time by 50-75% (one run vs N runs). This is the hardest to implement but highest potential for CPU-bound speedup.

### 5. Sequential vs Parallel External Tool Calls (MEDIUM impact) DO NOT IMPLEMENT the order of the steps is important. Though 1 and 2 could run in parallel this would not gian a lot of speed. Also downloading structures and sequences needs to be in this order because it first grabs all structures and then downloads sequences for the successfully downloaded structures. So these steps are also dependent on each other. (check this in the code)

**Problem:** In a single job, the pipeline runs sequentially:

1. evcouplings (CPU-bound, but uses multiprocessing internally)
2. Foldseek webserver (network I/O)
3. Download structures (network I/O)
4. Download sequences (network I/O)
5. Foldseek local + USalign (CPU-bound)
6. hmmsearch (CPU-bound)

Steps 1 and 2 could run in parallel since they are independent operations on different IDs. Similarly, steps 3 and 4 could be combined/concurrent.

**Solution:**

- Run evcouplings and Foldseek webserver in parallel (they use different resources)
- Combine structure + sequence downloads into a single concurrent phase: for each ID, download structure and sequence in parallel

**Estimated speedup:** 20-30% reduction in total job wall time by overlapping CPU and network phases.

### 6. Foldseek Local + USalign Run Redundantly (LOW-MEDIUM impact) (actually remove the foldseek call but keep the code so the user could use it if they want to by f.e. specifying this in the config. we could implement this in a later stage.)

**Problem:** Both tools compute structural similarity of ALL downloaded structures against the seed ROI. They serve complementary roles (Foldseek is fast but approximate; USalign is slower but more accurate), but USalign runs on ALL structures even if Foldseek already rejected most.

**Solution:** Use Foldseek as a pre-filter: only run USalign on structures that passed Foldseek (or structures where Foldseek score is close to the cutoff). This is partially done via the merge logic, but USalign still runs on everything.

### 7. CSV Files as Global State (LOW impact, correctness not speed) DO NOT IMPLEMENT and check this again. i think there are no race conditions since only the jobresultscollector in the main loop runs its methods seuqentially

**Problem:** `seeds.csv`, `used_seeds.csv`, and `All_Members_global.csv` are read and written by both the scheduler process and the collector. Parallel job completions can cause race conditions when writing these files. While not a speed issue directly, data corruption from races can cause incorrect convergence decisions, wasting compute.

**Solution:** Use SQLite with WAL mode for the global state. WAL supports concurrent readers and a single writer, which matches the pipeline's access pattern. SQLite also eliminates the need to read/write entire CSVs for incremental updates.

## Shared Directory Implementation Plan (Recommended)

This is your primary speed improvement. Here's a concrete plan:

### New directory layout

```
runs/run_<seed>_<timestamp>/
├── shared/
│   ├── structures/           ← *.cif files
│   └── sequences/
│       ├── records.fasta     ← combined FASTA (optional)
│       └── by_id/            ← individual FASTA per accession
├── data/                     ← (unchanged)
└── runs/                     ← (unchanged)
```

### Code changes required

1. **`run_true_members.py`:**
   - Add `shared_structures_dir` and `shared_sequences_dir` to paths dict
   - Create these directories in `set_up_dirs()`

2. **`JobScheduler.py`:**
   - Pass shared paths into `create_job_config()` so they end up in each job's config YAML

3. **`job_config_template.yaml`:**
   - Add two new fields: `global.shared_structures_dir` and `global.shared_sequences_dir`

4. **`Job.py`:**
   - Read new shared paths from config
   - Remove `self.alphafold_structures_dir` in favor of `self.shared_structures_dir`
   - Pass shared paths to `get_all_structures_and_sequences()`

5. **`get_structures_and_sequences.py`:**
   - `download_alphafold_structure()` already checks `os.path.exists(save_path)` — just point `save_path` at shared dir
   - For thread safety on shared filesystems: write to a temp path, then rename atomically
   - For sequence downloads, maintain individual FASTA files in `shared/sequences/by_id/`

6. **`utils.py`:**
   - Add helper: `def get_or_download_structure(uid, shared_dir, session) -> str`
   - Add atomic write: write to `.tmp` file, then `os.rename(tmp_path, final_path)`

### Thread-safety notes

Since Slurm jobs run on different compute nodes and share a network filesystem (likely Lustre or NFS on your HPC), file-level concurrency is the concern:

- `os.rename()` is atomic on POSIX filesystems (including NFS v3+ and Lustre)
- Write each structure/sequence to `shared/structures/<uid>.cif.tmp`, then rename to `shared/structures/<uid>.cif`
- If two jobs try to download the same file simultaneously, the last rename "wins" — both produce correct results
- For the records FASTA, use per-ID files instead of a single combined file to avoid append conflicts

## Summary of Recommendations (Priority Order)

| Priority | Change                                           | Effort | Estimated Speedup         |
| -------- | ------------------------------------------------ | ------ | ------------------------- |
| **P0**   | Shared structures directory (per-run)            | Medium | 40-60% download reduction |
| **P0**   | Shared sequences directory (per-run)             | Medium | 40-60% download reduction |
| **P1**   | Parallelize evcouplings + Foldseek webserver     | Medium | 20-30% wall time          |
| **P1**   | Fix Foldseek webserver to use current seed's ROI | Small  | Correctness, not speed    |
| **P2**   | USalign only on Foldseek pre-passed structures   | Small  | 30-50% USalign time       |
| **P2**   | SQLite for global state (correctness)            | Medium | Prevents corruption       |
| **P3**   | Single jackhmmer run + post-filter by bitscore   | Large  | 50-75% jackhmmer time     |
