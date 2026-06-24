# Analysis: Batch Downloading vs. Connection Tuning for UniProt Downloads

## The Problem

From the logs, 32,856 UniProt records were being fetched at **~0.48 records/sec** with 8 workers, resulting in an estimated **~19 hours** for a full download — far exceeding the 7-hour Slurm time limit. Jobs predictably failed around the 12,500-record mark after 7 hours.

## Root Cause

The root cause is **per-request overhead**, not rate-limiting or parallelism. Each accession triggers:

1. A TLS handshake (reused via session pooling, but still)
2. An HTTP GET round-trip to `rest.uniprot.org`
3. Server-side lookup
4. Response download

With 32k accessions, that's 32k round-trips. Even with perfect parallelism, round-trip latency dominates.

The `pool_connections=1` bug (serializing all 8 workers behind a single TCP connection) made an already bad situation worse, but fixing it alone won't make the math work:

- 32,856 records at even 5 records/sec (a 10× improvement) = **6,571 seconds ≈ 1.8 hours**
- 5 records/sec is optimistic for individual accession lookups
- In practice, the TLS + lookup overhead per request puts a hard ceiling around 2–5 req/s even with unlimited parallelism

## Approach Comparison

|                          | Kilo's Tuning Fix                     | Claude's Batch Fix                                | Combined    |
| ------------------------ | ------------------------------------- | ------------------------------------------------- | ----------- |
| **Strategy**             | Tune HTTP pool, timeouts, concurrency | Use UniProt `/stream` bulk endpoint               | Both        |
| **Requests for 32k IDs** | 32,000 (unchanged)                    | ~220 (150 IDs/batch)                              | ~220        |
| **Estimated time**       | 2–5 hours (optimistic)                | <10 minutes                                       | <10 minutes |
| **Solves root cause?**   | No — still 1 req/ID                   | Yes — 150× fewer requests                         | Yes         |
| **Risk**                 | Low — parameter changes only          | Medium — new API surface, query syntax dependency | Medium      |

### Kilo's Changes (applied to `utils.py`)

```
pool_connections:     1 → pool_maxsize     (allows actual parallel TCP connections)
retry total:          3 → 2                 (less retry time on failures)
backoff_factor:       1.0 → 0.5             (faster retry spacing)
read timeout:         30s → 15s             (fail faster on slow/stalled requests)
max_workers UniProt:  8 → 3                 (reduce pressure: 5 jobs × 8 = 40 → 5 × 3 = 15 conns)
max_workers AF DB:    16 → 8                (same principle)
```

**Verdict:** These are correct and should be kept. They reduce wasted time, avoid serialization bugs, and prevent aggressive parallelism from triggering server-side throttling. But alone, they are **insufficient** — the per-request overhead bottleneck remains.

### Claude's Batch Changes (proposed)

```
New functions:
  _chunked(seq, size) → generator
  _build_batch_query(ids, is_uniparc) → "accession:X OR accession:Y OR ..."
  _fetch_batch_fasta_text(ids, session, is_uniparc) → GET /uniprotkb/stream?query=...&format=fasta
  download_fasta_batch_to_cache(ids, cache_dir, session, is_uniparc)
  _download_records_batched(accessions, max_workers, cache_dir, ...)

Modified:
  get_records_from_accessions → new params: batch_size=150, use_batching=True
    - Tries batched first
    - Falls back to individual requests for any IDs the batch missed
    - Backward-compatible: use_batching=False returns old behavior
```

**Verdict:** This is the **correct solution** to the throughput problem. 32k requests → 220 requests is the order-of-magnitude change needed. The fallback mechanism is well-designed — individual retries for misses preserve correctness.

## Issues with Claude's Proposed Code

The proposed code is a good starting point, but has several issues that should be addressed before merging:

### 1. Kept the old `pool_connections=1` bug "alright lets use our fix here"

```python
adapter = HTTPAdapter(pool_connections=1, pool_maxsize=pool_maxsize, max_retries=retries)
```

The batch approach only uses a few concurrent workers (on the order of 220/8 ≈ 28 chunks), so the serialization bug matters less, but it should still be fixed to `pool_connections=pool_maxsize`.

### 2. Kept old retry/timeout values "lets lower them a bit"

```python
retry_kwargs = dict(total=3, backoff_factor=1, ...)
timeout=(5, 30)   # in _fetch_fasta_text
```

With batching, these matter less (fewer requests), but reducing them is still good defensive practice.

### 3. GET request may hit URL length limits "lets use post but do research first on how the post has to look like"

```python
def _fetch_batch_fasta_text(ids, session, is_uniparc):
    query = _build_batch_query(ids, is_uniparc)
    with session.get(base_url, params={"query": query, "format": "fasta"}, ...)
```

For 150 accessions like `A0A6P2YNG8`, the query string becomes:

```
accession:A0A6P2YNG8 OR accession:A0A6P2Z3X9 OR ... (×150)
```

That's ~150 × (10 + 3 + 10 + 3) ≈ 3,900 characters. HTTP GET URLs typically have a practical limit of ~2,000–8,000 characters depending on the server. This could break for accessions with longer IDs or at scale.

**Recommendation:** Use a `POST` request with a JSON body instead, or reduce `batch_size` to 100 to stay safe. UniProt's `/stream` endpoint accepts POST with `{"ids": [...], "format": "fasta"}`.

### 4. `extract_uid_from_recordID` called on `record.id` which is already a UniProt ID "lets keep the acccession extraction"

In `download_fasta_batch_to_cache`:

```python
for record in records:
    uid = extract_uid_from_recordID(record.id) or record.id
```

`SeqIO.parse` on a FASTA from UniProt returns records whose `.id` is the full header (e.g., `sp|P12345|PROT_HUMAN`). `extract_uid_from_recordID` regex-extracts the accession from this. This is correct behavior but worth noting — if UniProt changes their FASTA header format, this regex could break.

### 5. User-Agent header "fill with waronx@gmx.de"

```python
session.headers.update({"User-Agent": "True-Members-Pipeline/1.0 (contact: REPLACE_WITH_YOUR_EMAIL)"})
```

Good practice. The contact email should be filled in before deploying.

### 6. `_download_records_batched` iterates all records to build `found_ids` "skip for now"

```python
for record in chunk_records:
    found_ids.add(extract_uid_from_recordID(record.id) or record.id)
```

With 32k records, this is fine. But the `extract_uid_from_recordID` call per-record could be avoided by pre-computing a set of expected IDs and checking membership — slightly more efficient but not critical.

### 7. Duplicate `get_records_from_accessions` return path inconsistency

```python
if missing:
    fallback_records = _download_records_parallel(missing, max_workers, cache_dir, [])
else:
    fallback_records = []
return cached_records + batched_records + fallback_records
```

When `use_batching=False`, it calls `_download_records_parallel(uncached_accessions, max_workers, cache_dir, cached_records)`. Note the different semantics: with batching, `cached_records` is NOT passed to fallback, whereas without batching it IS passed. The `initial_records` parameter in `_download_records_parallel` is only used for progress logging (calculating `target_total`, `success`), so this is a cosmetic inconsistency, not a bug.

## Recommended Implementation Plan

Combine both approaches:

### Keep from Kilo's changes:

1. `pool_connections = pool_maxsize` (non-negotiable fix)
2. `total=2, backoff_factor=0.5` (tight retries)
3. `timeout=(5, 15)` for individual fallback requests
4. Reduced `max_workers` defaults (3 for UniProt, 8 for AlphaFold) "lets do 4 for uniprot here"

### Adopt from Claude's changes (with fixes):

1. Batch downloading via `/stream` endpoint
2. Fallback to individual requests for missed IDs
3. `use_batching` toggle for backward compatibility
4. User-Agent header

### Fixes to apply to Claude's code:

1. Use `pool_connections=pool_maxsize` instead of `pool_connections=1`
2. Use `POST` instead of `GET` for batch queries to avoid URL length limits:
   ```python
   session.post(base_url, data={"ids": ids, "format": "fasta"}, timeout=(5, 60))
   ```
   This eliminates the `_build_batch_query` function entirely and is safer.
3. Keep tighter retry/timeout values from Kilo's changes

## Performance Estimate (Combined)

| Metric                | Before   | After (Kilo only)  | After (Combined)             |
| --------------------- | -------- | ------------------ | ---------------------------- |
| Requests for 32k IDs  | 32,000   | 32,000             | ~220 batches + ~50 fallbacks |
| Effective req/sec     | ~0.5     | ~2–5               | N/A (batched, server-side)   |
| Estimated time        | 19 hours | 2–5 hours          | <10 minutes                  |
| Fits in 7-hour limit? | No       | Maybe (unreliable) | Yes                          |

## Additional Consideration: AlphaFold Structure Downloads "lets not work on this here"

The AlphaFold structure downloads (`run_download_tasks`) also do one request per accession. There is no bulk endpoint for AlphaFold CIFs. However:

- Structure downloads are cached in `shared/structures/` and skipped if the file already exists
- Only one job needs to do the full download; subsequent jobs in the same pipeline run reuse cached files
- If structure downloads are also a bottleneck, consider splitting the accession list across a Slurm job array

This is outside the scope of the current fix but worth monitoring.

## Doneness Checklist

- [ ] Decide whether to implement batch downloading
- [x] ~~If yes: use `POST` to `/uniprotkb/stream` and `/uniparc/stream` instead of `GET` with query strings~~ — POST returns 500 from UniProt server; GET with `accession_id:` works reliably. URL length is ~2.5KB for 100 IDs, well under limits.
- [ ] If yes: merge Kilo's connection pool + timeout fixes into the batch-enabled code
- [ ] If yes: smoke-test with 20–50 accessions to verify query syntax and response parsing
- [ ] If no: revert `max_workers` back to 8 (since with batching it's less critical) but keep pool_connections and timeout fixes
- [ ] Fill in a real contact email for the User-Agent header
- [ ] Run `pytest tests/ -v` after changes
