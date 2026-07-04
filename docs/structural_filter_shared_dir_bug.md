# Bug: Structural Filter Processes All Structures in Shared Directory

**Date:** 2026-07-04
**Status:** Confirmed

## Summary

`run_structural_similarity_filter()` in `True_Members/scripts/structural_filter.py` runs Foldseek and/or USalign against **all** structures in the shared directory (`shared_structures_dir`), not just the structures relevant to the current job. The actual filtering to the job's IDs (`uids_with_region`) happens only *after* the expensive structural alignment computation has already completed on all structures.

## Root Cause

### The Call Site (`True_Members/Job.py:105-112`)

```python
seed_containing_sequences_with_tm = run_structural_similarity_filter(
    self.GLOBAL_SEED_structure_roi,
    self.stats,
    self.TM_SCORE_CUTOFF,
    self.shared_structures_dir,       # ← FULL shared dir, not filtered
    self.foldseek_dir,
    self.usalign_dir,
    uids_with_region,                 # ← job-specific IDs, only used post-hoc
    self.STRUCTURAL_FILTER_TOOL)
```

### Three Problems in `structural_filter.py`

1. **Foldseek searches all structures** (lines 121, 139):

   ```python
   foldseek_output = run_foldseek(query_structure_roi, structures_dir, foldseek_dir)
   ```

   `structures_dir` is the full `shared_structures_dir`, which accumulates `.cif` files across ALL jobs and ALL iterations. Foldseek `easy-search` receives a directory and searches every structure in it.

2. **USalign target list includes all structures** (lines 109-115, 122, 148):

   ```python
   def create_target_list_for_usalign(target_structures_dir, usalign_dir):
       target_list_path = f"{usalign_dir}/target_list.txt"
       with open(target_list_path, "w") as f:
           for file in os.listdir(target_structures_dir):
               if file.endswith(".cif") or file.endswith(".pdb"):
                   f.write(file + "\n")
       return target_list_path
   ```

   This lists every `.cif`/`.pdb` file in the full shared directory, not just the job's subset.

3. **Post-hoc filtering is too late** (lines 160-162):

   ```python
   total_filtered_df = total_filtered_df[
       total_filtered_df["id"].isin(uids_with_region["id"])
   ]
   ```

   This filter discards results for structures the job didn't discover, but only *after* expensive TM-score alignment has already been computed for all of them.

## Impact

| Aspect | Detail |
|--------|--------|
| **Correctness** | Not affected — the post-hoc filter at lines 160-162 discards extraneous results, so the final output is correct |
| **Performance** | Severely degraded — the structural alignment is the most expensive step in the pipeline (Foldseek/USalign TM-score computation), and it runs against O(N_jobs × N_structures_in_shared_dir) structures instead of O(N_job's own structures) |
| **Scaling** | Gets worse over iterations as `shared_structures_dir` accumulates more structures from previous jobs |

### Concrete example:
- Iteration 1: shared dir has 500 structures, job found 400 → alignment runs on 500 (25% waste)
- Iteration 3: shared dir has 1500 structures, job found 300 → alignment runs on 1500 (400% waste)

## Fix Plan

Two changes are needed:

### 1. Create a job-specific sub-directory with only the relevant structures

In `Job.run()`, before calling `run_structural_similarity_filter()`, create a temp directory (e.g., `<root_dir>/filtering/target_structures/`) containing symlinks to only the `.cif` files whose IDs appear in `uids_with_region["id"]`. Pass this directory instead of `self.shared_structures_dir`.

### 2. Update `run_structural_similarity_filter()` to accept a pre-filtered directory

Either:
- **Option A:** Keep the function signature the same but change `Job.py` to pass a filtered directory, OR
- **Option B:** Modify the function to accept `uids_with_region` and internally create the symlink directory before calling Foldseek/USalign

**Recommended: Option B** — keeps the filtering logic close to where it matters and avoids leaking temporary-directory responsibilities to the caller.

### Files to modify

| File | Change |
|------|--------|
| `True_Members/scripts/structural_filter.py` | `run_structural_similarity_filter()`: create a temp symlink dir from `uids_with_region["id"]` + `structures_dir`, use it for alignment |
| `True_Members/scripts/structural_filter.py` | `create_target_list_for_usalign()`: accept `uids_with_region` as an additional parameter to generate a filtered target list (or remove it if using symlink approach) |
| `True_Members/Job.py` | No change needed if Option B is used |

### Alternative (simpler, USalign-only)

For USalign, a simpler fix is to change `create_target_list_for_usalign` to accept the list of IDs directly and only write matching files to the target list:

```python
def create_target_list_for_usalign(uids_with_region, structures_dir, usalign_dir):
    target_list_path = f"{usalign_dir}/target_list.txt"
    id_set = set(uids_with_region["id"])
    with open(target_list_path, "w") as f:
        for uid in id_set:
            cif_path = f"{structures_dir}/{uid}.cif"
            if os.path.exists(cif_path):
                f.write(f"{uid}.cif\n")
    return target_list_path
```

This fixes USalign but not Foldseek, because Foldseek's `easy-search` takes a directory, not a list file. For Foldseek, the symlink directory approach is still needed.