# Analysis: Fallback Slowness & Local UniProt Database

## Problem: Jobs Get Progressively Slower

After the batch download completes, any IDs not returned by the batch (obsolete accessions, secondary accessions, UniParc IDs not in UniProtKB) are retried individually via `_download_records_parallel` in `utils.py:154`. The log shows this fallback running at ~1 sequence/second.

### What's Happening

```
Job 1: batch (fast) → fallback (200 IDs × ~1s each = 200s)
Job 2: batch (fast) → fallback (200 IDs × ~1s each = 200s, but slower due to cumulative UniProt pressure)
Job 3: batch (fast) → fallback (gets even slower)
...
```

Each job's fallback phase makes hundreds of individual `GET /uniprotkb/{uid}.fasta` requests. With `max_workers=3` and 5 concurrent Slurm jobs, that's 15 simultaneous connections to `rest.uniprot.org`, all making individual-accession requests. This is likely the same pattern that caused the original 19-hour problem, just at smaller scale.

### Root Cause

The batch approach (100 IDs per request) fixes the primary throughput issue, but the fallback path still uses one-request-per-ID. These fallback IDs accumulate because:

1. **Obsolete/secondary accessions** — JackHMMER returns secondary accessions that UniProt's `/stream` endpoint doesn't resolve (returns 200 but silently skips them in the output)
2. **UniParc IDs** — Separate endpoint, separate batching, may have their own resolution issues
3. **Transient server errors on chunks** — A batch that fails (5xx, timeout) dumps ALL 100 IDs into the fallback path

### Likely Cause of Throttling

**UniProt rate limiting** is the most probable cause. Evidence:
- The slowdown is progressive (later jobs suffer more) — consistent with rate-limit counters accumulating
- The slowdown is per-request, not per-batch — the batch requests themselves remain fast
- UniProt's public API has rate limits even if not publicly documented

**MPCDF Viper network throttling** is less likely because:
- Batch requests (GET with 100-query strings) remain fast even when individual requests are slow
- Network throttling would affect both paths equally
- The `/stream` endpoint likely has different rate-limit counters than the per-accession endpoint

## Opportunity: Local UniProt Database

The evcouplings config template (`evcouplings_config_template.yaml:271`) references:

```
/mpcdf/soft/RHEL_9/packages.norpm/x86_64/evcouplings/0.2.1/data/DB/uniprot/uniprot_2026_4.fasta
```

This is a **single FASTA file** containing all UniProt sequences (installed by MPCDF helpdesk for jackhmmer). It's already HMMER-indexed (`.h3f`, `.h3i`, `.h3m`, `.h3p` index files alongside).

### Using It for Sequence Fetching

Instead of hitting `rest.uniprot.org` at all, we could extract sequences directly from this local FASTA file. Benefits:

| | REST API (current) | Local FASTA (proposed) |
|---|---|---|
| Speed (per ID) | 0.1–1s (network) | <0.001s (filesystem) |
| Rate limited? | Yes | No |
| Network dependent? | Yes | No |
| Fails on slurm time limit? | Yes | No |
| Handles obsolete IDs? | Returns 404 or skips | FASTA has `sp|P12345|...` headers; we parse the primary accession |
| UniParc coverage? | Yes (separate endpoint) | No (UniProtKB only) |

### How a FASTA Lookup Would Work

1. Build an index mapping accession → byte offset (one-time, ~1–2 seconds for a 100GB FASTA)
2. For each accession: seek to offset, read until next `>` header, parse the record
3. Cache to per-accession `.fasta` files as before (backward compatible with existing cache)

The index can be built with:
- **`pyfaidx`** — Python library, reads `.fai` index file (same format as `samtools faidx`)
- **Manual index** — Read the FASTA once, build a dict of `{accession: (byte_offset, byte_length)}`
- **`samtools faidx`** — Pre-index with `samtools faidx uniprot.fasta`, then query individual accessions

### Drawbacks

1. **UniParc IDs** (`UPI` prefix) are not in UniProtKB — still need REST API for those
2. **FASTA size** — `uniprot_2026_4.fasta` is likely 80–120GB; building the index requires reading the entire file once
3. **File access** — Must be readable from all Slurm nodes (it's on `/mpcdf/soft/` shared filesystem, so this is fine)
4. **Secondary accessions** — UniProt FASTA headers have format `>sp|P00533|EGFR_HUMAN ...`. The primary accession is the second pipe-delimited field. Secondary/obsolete accessions won't have their own entries — they need resolution via the REST API or a mapping file

## Recommended Approach

### Phase 1: Reduce Fallback Pressure (immediate)

Reduce the number of IDs that hit the fallback path:

1. **Increase batch size** from 100 to 200 — UniProt's `/stream` endpoint can handle this. Fewer batches = fewer chunks that could fail and dump IDs into fallback.
2. **Retry failed batches** — Instead of sending all IDs from a failed chunk to individual fallback, split the chunk and retry the sub-chunks. Only IDs that fail in a minimal sub-chunk (e.g., batch_size=5) go to individual fallback.
3. **Pre-filter obviously invalid IDs** — JackHMMER may return IDs that aren't valid UniProt accession formats. Filter these out before batching to avoid 400s that dump the entire batch into fallback.

### Phase 2: Local FASTA Lookup (recommended, higher effort)

Add a `LocalFastaLookup` class that:

1. On first use, builds/loads a `.fai` index for the FASTA file at the path from evcouplings config
2. `get_sequence(accession)` → returns SeqRecord or None
3. Falls back to REST API for UniParc IDs and IDs not found in the local FASTA

This eliminates the network dependency entirely for the common case. With 7,408 UniProt records per job, even 1 second saved per record = 2 hours saved per job.

The FASTA path could be read from the master config (`run_true_members_config.yaml`) or discovered from the evcouplings config template.

## Questions to Decide

1. **Is the local FASTA readable from compute nodes?** The path is on `/mpcdf/soft/` which is typically a shared mount. Should be accessible, but verify on a compute node.

2. **Which approach do you want to implement first?**
   - Retry-failed-batches (low effort, helps fallback path)
   - Local FASTA lookup (higher effort, eliminates REST API for UniProtKB entirely)
   - Both

3. **Should we keep the REST API as backup?** Yes — for UniParc IDs and any accessions not found in the local FASTA (e.g., very new proteins added after the FASTA was built).