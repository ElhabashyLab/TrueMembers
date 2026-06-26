# Fix: Limit Failed Seed Retries to One Attempt

## Problem

In `True_Members/JobScheduler.py:remove_seeds_of_failed_jobs_from_used_seeds()` (line 184-216), when a job fails, its seed is unconditionally removed from `used_seeds.csv` (line 202), making it available for re-selection in future iterations. There is **no limit** on how many times a failed seed can be retried.

This means a seed that consistently fails (e.g., due to a structural issue with the protein, a corrupted input, or an external service outage) will be retried indefinitely, wasting cluster resources and never allowing the pipeline to converge.

## Root Cause Location

- **`JobScheduler.py:184-216`** — `remove_seeds_of_failed_jobs_from_used_seeds()` frees the seed regardless of how many times it has already been freed/retried
- **`JobScheduler.py:349-355`** — `get_highest_scoring_seed()` selects the highest-scoring available seed; a freed seed that was previously highest-scoring will be selected again immediately

## Fix Strategy

Track which seeds have already been granted one retry in an in-memory `set` on the `JobScheduler` instance. When a job fails:

1. If the seed has **not** used its retry yet: add it to the set, remove from `used_seeds` (free it for retry) — same as current behavior
2. If the seed has **already** used its retry: keep it in `used_seeds`, mark it as exhausted (set `new_true_members = 0`), log a warning, do **not** free it

This limits each seed to at most 1 retry (2 total attempts: initial + 1 retry).

## Data Flow

```
Job fails for seed X
    │
    ├─► X not in self._seeds_retried?  →  add X to set, remove from used_seeds → seed is freed for retry
    │
    └─► X already in self._seeds_retried?  →  set new_true_members=0, keep in used_seeds, log warning → seed is permanently exhausted
```

## Implementation Steps

### 1. Add retry tracker to `JobScheduler.__init__()`

In `JobScheduler.py`, add a new instance attribute in `__init__()` (after line 56):

```python
        self._seeds_retried: set[str] = set()
```

### 2. Modify `remove_seeds_of_failed_jobs_from_used_seeds()`

Replace the current logic (lines 184-216) with the following:

```python
    def remove_seeds_of_failed_jobs_from_used_seeds(self):
        failed_job_ids = [
            job_id for job_id in self.jobs.keys()
            if self.jobs[job_id]["status"] == JobStatus.FAILED
        ]
        if not failed_job_ids:
            return

        used_seeds_df = pd.read_csv(self.used_seeds_path)

        for job_id in failed_job_ids:
            try:
                job_config_path = self.jobs[job_id]["config_path"]
                job_config = read_yaml_config(job_config_path)
                failed_iteration = job_config["global"]["iteration"]
                seed_id_of_job = job_config["run"]["seed_uid"]

                if seed_id_of_job in self._seeds_retried:
                    logger.warning(
                        "Seed '%s' (iteration %d, job %d) failed for the second time. "
                        "Keeping it in used_seeds — no more retries.",
                        seed_id_of_job, failed_iteration, job_id
                    )
                    used_seeds_df.loc[used_seeds_df["iteration"] == failed_iteration, "new_true_members"] = 0
                    self.jobs[job_id]["status"] = JobStatus.HARVESTED
                    continue

                self._seeds_retried.add(seed_id_of_job)
                before = len(used_seeds_df)
                used_seeds_df = used_seeds_df[used_seeds_df["iteration"] != failed_iteration]
                removed = before - len(used_seeds_df)

                if removed:
                    logger.info("Freed seed '%s' (iteration %d) from used_seeds after job %d failed (retry 1/1).",
                                seed_id_of_job, failed_iteration, job_id)
                else:
                    logger.warning("No used_seeds row found for failed job %d (seed '%s', iteration %d).",
                                   job_id, seed_id_of_job, failed_iteration)
            except Exception as e:
                logger.error("Failed to handle failed job %d: %s", job_id, e)
            finally:
                if self.jobs[job_id]["status"] != JobStatus.HARVESTED:
                    self.jobs[job_id]["status"] = JobStatus.HARVESTED

        used_seeds_df.to_csv(self.used_seeds_path, index=False)
```

Key changes:
- Check `seed_id_of_job in self._seeds_retried` before freeing
- On first failure: add to `_seeds_retried`, then free as before
- On second failure: set `new_true_members=0` to mark it exhausted, keep it in used_seeds
- Only set `HARVESTED` status in `finally` if not already set inside the loop

### 3. Edge Cases Covered

| Scenario | Behavior |
|----------|----------|
| Seed fails once | Freed from used_seeds, can be retried |
| Seed fails twice | Permanently kept in used_seeds with `new_true_members=0`, never selected again |
| Seed fails, retry succeeds | Normal completion; `new_true_members` updated with real count, seed stays in used_seeds |
| Seed fails, retry fails (different iteration) | Caught by `_seeds_retried` set, treated as second failure |
| Multiple different seeds fail | Each tracked independently in `_seeds_retried` |

### 4. State Persistence Note

`_seeds_retried` is an in-memory set. This is correct for the expectation: the retry limit is per pipeline run. If the entire pipeline is restarted from scratch, seeds start fresh (0 retries). This is desirable — a transient HPC outage might cause many initial failures that should be retried in a fresh run.

## Files Changed

- `True_Members/JobScheduler.py` — Add `_seeds_retried` attribute and modify `remove_seeds_of_failed_jobs_from_used_seeds()`

## Testing

- Unit test: mock `used_seeds.csv` with a failing seed, call the method twice, verify the seed is freed on first call but kept on second call
- Integration test: simulate two consecutive failures for the same seed, verify it's not selected for a third job