# Integration Plan: Local UniProt DB for Sequence Fetching

## Overview

Replace round-trips to `rest.uniprot.org` with direct filesystem reads from the already-available UniProt FASTA database at `/mpcdf/soft/.../uniprot_2026_4.fasta`. REST API remains as fallback for any accession not found locally.

## Design

### New File: `True_Members/scripts/uniprot_local.py`

A `UniprotLocalDB` class responsible for indexing the FASTA once and extracting individual sequences by accession.

```
UniprotLocalDB
├── __init__(fasta_path, index_dir: str | None)
├── _build_index() → writes {fasta_path}.idx text file
├── _load_index() → dict[accession, (byte_offset, byte_length)]  
├── _lookup_offset(accession) → byte_offset | None
└── get_sequence(accession) → SeqRecord | None
```

#### Index Format

After the first call, a `.idx` text file is created alongside the FASTA with one line per sequence:

```
P00533\t123456789\n
P04637\t123458012\n
...
```

This file is sorted by accession for binary search (future optimization). Initially, the entire index is loaded into a Python dict. For a FASTA with ~2M sequences, the dict is ~200MB — well within HPC node memory (64–256GB typical).

#### FASTA Header Parsing

Use the existing `extract_uid_from_recordID()` from `scripts/utils.py:181` — it already handles the `>sp|P00533|EGFR_HUMAN ...` format via the `_UNIPROT_ACCESSION_RE` regex. No new parsing logic needed.

#### Sequence Extraction

```python
def get_sequence(self, accession: str) -> SeqRecord | None:
    offset = self._lookup_offset(accession)
    if offset is None:
        return None
    with open(self.fasta_path) as f:
        f.seek(offset)
        # Read lines until next '>' header
        lines = []
        for line in f:
            if line.startswith('>') and lines:
                break
            lines.append(line)
    return next(SeqIO.parse(StringIO(''.join(lines)), 'fasta'))
```

### Changes to Existing Files

#### 1. `run_true_members_config.yaml` — add field

```yaml
run:
  seed_uid: "P00533"
  # ... existing fields ...
  uniprot_fasta_path: "/mpcdf/soft/RHEL_9/packages.norpm/x86_64/evcouplings/0.2.1/data/DB/uniprot/uniprot_2026_4.fasta"
```

Set to empty string `""` to disable local DB and use REST API only. Do NOT use `null`/`None` — empty string is the sentinel.

#### 2. `True_Members/template_files/job_config_template.yaml` — add field

```yaml
run:
  uniprot_fasta_path: {{ UNIPROT_FASTA_PATH }}
```

The JobScheduler fills this from the master config.

#### 3. `True_Members/scripts/utils.py` — modify `get_records_from_accessions`

Add `local_fasta_path` parameter:

```python
def get_records_from_accessions(accessions, max_workers=8, cache_dir=None,
                                 batch_size=100, use_batching=True,
                                 local_fasta_path=""):
```

New flow:

```
get_records_from_accessions(accessions, ..., local_fasta_path=path)
  │
  ├── _collect_cached_records() → cached + uncached
  │
  ├── if local_fasta_path (non-empty):
  │     local_records, still_missing = _download_from_local_db(
  │         uncached, local_fasta_path, cache_dir)
  │   else:
  │     local_records = []
  │     still_missing = uncached
  │
  ├── if still_missing and use_batching:
  │     batched, missing = _download_records_batched(still_missing, ...)
  │     if missing:
  │         retried, still_missing = _download_records_batched(missing, ...)
  │         batched += retried
  │         # still_missing is discarded — no individual fallback
  │   elif still_missing:
  │     batched = _download_records_parallel(still_missing, ...)
  │   else:
  │     batched = []
  │
  └── return cached + local_records + batched
```

**Key change**: individual `_download_records_parallel` fallback after batch is removed. IDs not returned by the batch (obsolete/nonexistent accessions) are retried once via a second batch, then silently skipped. This eliminates the 1-seq/sec slowdown when hundreds of obsolete accessions would otherwise be fetched one-by-one.

`max_workers` default restored to `8` (from `3`) since the individual-request bottleneck is no longer a concern.

#### 4. `True_Members/scripts/get_structures_and_sequences.py` — pass through

```python
def get_all_structures_and_sequences(seed, bit_scores, stats, ..., local_fasta_path=None):
    ...
    records = get_records_from_accessions(
        combined_unique_ids_with_regions["id"].tolist(),
        cache_dir=shared_sequences_dir,
        local_fasta_path=local_fasta_path,
    )
```

#### 5. `True_Members/Job.py` — read config and pass through

```python
def __init__(self, ...):
    ...
    self.uniprot_fasta_path = self.true_members_config["run"].get("uniprot_fasta_path")

def run(self):
    ...
    uids_with_region_full_fasta_path, uids_with_region = get_all_structures_and_sequences(
        ...,
        local_fasta_path=self.uniprot_fasta_path,
    )
```

#### 6. `True_Members/JobScheduler.py` — propagate to per-job configs

In `create_job_config()`, add `uniprot_fasta_path` from master config to the generated per-job YAML.

### Data Flow Diagram

```
Master Config (run_true_members_config.yaml)
    │  uniprot_fasta_path: "/mpcdf/soft/.../uniprot_2026_4.fasta"
    ▼
JobScheduler.create_job_config()
    │  Writes uniprot_fasta_path into per-job config
    ▼
Per-job Config (job_config_*.yaml)
    │  
    ▼
Job.__init__()
    │  self.uniprot_fasta_path = config["run"].get("uniprot_fasta_path")
    ▼
Job.run()
    │  
    ▼
get_all_structures_and_sequences(..., local_fasta_path=path)
    │  
    ▼
get_records_from_accessions(ids, local_fasta_path=path)
    │
    ├── [1] _download_from_local_db()
    │         │
    │         ├── UniprotLocalDB.get_sequence(id) → SeqRecord | None
    │         └── Write to cache_dir/{id}.fasta (atomic .tmp → rename)
    │
    └── [2] if still missing → REST API (batch + individual fallback)
```

### Index Lifecycle

```
First call:
  check {fasta_path}.idx → not found
  → scan FASTA, build index, write .idx file (~2-5 min one-time cost)
  → load index into memory (~200MB for 2M seqs)

Subsequent calls (same job or later jobs):
  check {fasta_path}.idx → found
  → load index from text file (~5-10s for 2M entries)
  → each get_sequence() is a single seek + sequential read (~0.1ms)
```

### REST API Fallback Scenarios

The REST API is still used for:

| Case | Reason |
|---|---|
| `UPI` prefix IDs | UniParc, not in UniProtKB FASTA |
| Accession not in local index | Newer than FASTA build date |
| `local_fasta_path` is `""` (empty) | Backward compatible — pure REST API mode |

IDs not returned by the batch are retried once via a second batch, then silently skipped — no individual fallback. This avoids wasting time on obsolete accessions that will never resolve.

### Performance Estimate

| Path | Per 1000 accessions | 7408 accessions (from your run) |
|---|---|---|
| REST API (old, individual) | ~33 minutes | ~240 minutes |
| REST API (with batching) | ~30 seconds | ~3.5 minutes |
| REST API (batching, then fallback) | ~3 minutes | ~22 minutes |
| **Local FASTA** | **~0.5 seconds** | **~4 seconds** |
| Local FASTA → REST fallback (10% miss) | ~18 seconds | ~2 minutes |

The local FASTA eliminates the network bottleneck entirely and makes the fallback slowness problem irrelevant for UniProtKB accessions.

### Test Plan

1. **Unit test**: `tests/test_uniprot_local.py` — index a small mock FASTA, verify accession lookup via `extract_uid_from_recordID`, test missing accessions return `None`

### Implementation Order

1. `utils.py` — remove individual fallback from batch path; restore `max_workers=8`; change `local_fasta_path` default to `""`
2. `uniprot_local.py` — `UniprotLocalDB` class with index builder and `get_sequence()` using `extract_uid_from_recordID`
3. `utils.py` — add `_download_from_local_db()` function, wire into `get_records_from_accessions`
4. `get_structures_and_sequences.py` — pass `local_fasta_path` through
5. `Job.py` — read config field (empty string sentinel), pass through
6. `JobScheduler.py` — propagate to per-job configs
7. Config files — add `uniprot_fasta_path` field
8. `tests/test_uniprot_local.py` — unit tests with mock FASTA