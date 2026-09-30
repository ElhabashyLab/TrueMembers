"""
SequenceIdentityCalculator
 
Incrementally tracks pairwise sequence identity between search seeds and
detected hits across iterations of an iterative homolog search (e.g. jackhmmer).
 
Design
------
Per call (one finished search iteration), only the NEW cells of the
(seeds x hits) identity matrix are computed in two steps:
 
    1. current_seed        vs  (members_so_far + new_members)   -> new row
    2. previous_used_seeds vs  new_members                       -> new column

The sequence identity is computed with MMseqs2's `search` command.
Members are kept as a persistent MMseqs2 database on disk (under work_dir) 
that only ever grows by merging in the small new increment
(`createdb` on just new_members_fasta_path, then`concatdbs` to merge) 
 
MMseqs2 may not report an alignment for every pair (its k-mer prefilter can
miss very divergent pairs even at high sensitivity). Such pairs are recorded
with pident = 0.0 and flagged via `detected = False` so "below detection
threshold" stays distinguishable from "confirmed low identity" downstream.
"""

import logging
import os
import re
import subprocess
import shutil
import tempfile
import uuid
from pathlib import Path
from itertools import product
from typing import List, Optional, Tuple, Union
 
import pandas as pd

logger = logging.getLogger(__name__)
 
 
class SequenceIdentityCalculator:

    MMSEQS_OUTPUT_COLUMNS = ["query", "target", "pident"]
    MEMBERS_DB_PREFIX = "members_db_gen"

    def __init__(self, work_dir: str, mmseqs_binary: str = "mmseqs",
                 sensitivity: float = 7.5, max_evalue: str = "10000"):
        """
        work_dir:      directory that persists the growing members MMseqs2
                        database across calls/process restarts. Pass the same
                        work_dir on every call to resume state.
        mmseqs_binary:  name/path of the mmseqs executable.
        sensitivity:    -s value passed to mmseqs (higher = more sensitive, slower).
        max_evalue:     effectively-unfiltered e-value cutoff. Kept as a large
                        finite number rather than 'inf' since MMseqs2's float
                        parsing of 'inf' is not guaranteed across versions.
        """
        if shutil.which(mmseqs_binary) is None:
            raise RuntimeError(
                f"'{mmseqs_binary}' not found on PATH. Is MMseqs2 installed/loaded?"
            )
        self.mmseqs_binary = mmseqs_binary
        self.sensitivity = sensitivity
        self.max_evalue = max_evalue

        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)

        self._members_generation, self._members_db = self._discover_latest_members_db()
        logger.debug("SequenceIdentityCalculator initialized: work_dir=%s, generation=%d, db=%s",
                     self.work_dir, self._members_generation,
                     self._members_db.name if self._members_db else "None")
 
 
    def compute_sequence_identity(self,
                                   output_csv_path: str,
                                   current_seed_fasta_path: str,
                                   new_members_fasta_path: str,
                                   previous_used_seeds_fasta_path: str) -> pd.DataFrame:
        """
        Computes the identity pairs newly introduced by this iteration and
        appends them to output_csv_path (created with header if it doesn't
        exist yet). Returns the DataFrame of rows just computed.
 
        Columns written: seed_id, hit_id, pident, detected
        """
 
        with tempfile.TemporaryDirectory(prefix="seqid_scratch_") as scratch_str:
            scratch_dir = Path(scratch_str)
            new_rows = []

            has_new_members = self._fasta_has_sequences(new_members_fasta_path)
            has_seed = self._fasta_has_sequences(current_seed_fasta_path)
            has_prev_seeds = self._fasta_has_sequences(previous_used_seeds_fasta_path)
            logger.info("Computing sequence identity for collected hits.")
            logger.debug("compute_sequence_identity: has_seed=%s, has_new_members=%s, has_prev_seeds=%s, "
                         "persistent_db=%s (gen %d)",
                         has_seed, has_new_members, has_prev_seeds,
                         self._members_db.name if self._members_db else "None",
                         self._members_generation)

            new_members_db = None
            if has_new_members:
                new_members_db = self._create_db(
                    new_members_fasta_path, scratch_dir / "new_members_db"
                )

            # --- Step 1: current seed vs all members so far (persisted old + new) ---
            if has_seed and (self._members_db is not None or has_new_members):
                seed_db = self._create_db(current_seed_fasta_path, scratch_dir / "seed_db")
                seed_id = self._get_fasta_ids(current_seed_fasta_path)

                if self._members_db is not None:
                    old_ids = self._get_db_ids(self._members_db)
                    step1_old = self._search_all_pairs(
                        seed_db, self._members_db, seed_id, old_ids,
                        scratch_dir / "step1_old",
                    )
                    new_rows.append(step1_old)
                    logger.debug("Step 1 (old): %d seed x %d old members → %d pairs (%d detected)",
                                 len(seed_id), len(old_ids), len(step1_old), step1_old["detected"].sum())
                if new_members_db is not None:
                    new_ids = self._get_fasta_ids(new_members_fasta_path)
                    step1_new = self._search_all_pairs(
                        seed_db, new_members_db, seed_id, new_ids,
                        scratch_dir / "step1_new",
                    )
                    new_rows.append(step1_new)
                    logger.debug("Step 1 (new): %d seed x %d new members → %d pairs (%d detected)",
                                 len(seed_id), len(new_ids), len(step1_new), step1_new["detected"].sum())

            # --- Step 2: previously used seeds vs newly detected hits only ---
            if has_prev_seeds and new_members_db is not None:
                prev_seeds_db = self._create_db(
                    previous_used_seeds_fasta_path, scratch_dir / "prev_seeds_db"
                )
                prev_seed_ids = self._get_fasta_ids(previous_used_seeds_fasta_path)
                new_member_ids = self._get_fasta_ids(new_members_fasta_path)
                step2 = self._search_all_pairs(
                    prev_seeds_db, new_members_db,
                    prev_seed_ids, new_member_ids,
                    scratch_dir / "step2",
                )
                new_rows.append(step2)
                logger.debug("Step 2: %d prev seeds x %d new members → %d pairs (%d detected)",
                             len(prev_seed_ids), len(new_member_ids), len(step2), step2["detected"].sum())

            # --- Persist: merge new members into the running DB ---
            if new_members_db is not None:
                self._merge_new_members(new_members_db)

            result_df = (
                pd.concat(new_rows, ignore_index=True) if new_rows
                else pd.DataFrame(columns=["seed_id", "hit_id", "pident", "detected"])
            )

            self._append_to_csv(result_df, output_csv_path)
            return result_df
 
 
    def _discover_latest_members_db(self) -> Tuple[int, Optional[Path]]:
        """Finds the highest-numbered members_db_gen{N} already on disk, if any."""
        pattern = re.compile(rf"^{self.MEMBERS_DB_PREFIX}(\d+)$")
        best_gen, best_path = -1, None
        for entry in self.work_dir.iterdir():
            match = pattern.match(entry.name)
            if match:
                gen = int(match.group(1))
                if gen > best_gen:
                    best_gen, best_path = gen, entry
        if best_path is not None:
            return best_gen, best_path
        return -1, None
 
 
    def _next_members_db_path(self) -> Path:
        return self.work_dir / f"{self.MEMBERS_DB_PREFIX}{self._members_generation + 1}"


    def _merge_new_members(self, new_members_db: Path) -> None:
        """
        Merges new_members_db into the persistent members DB, writing the
        result under a fresh generation name (unique temp name first, then
        renamed into place, to avoid a half-written DB surviving a killed job).
        """
        target_path = self._next_members_db_path()
        tmp_name = self.work_dir / f".tmp_{uuid.uuid4().hex}_{target_path.name}"

        if self._members_db is None:
            logger.debug("Merging new members: first generation, copying %s → %s",
                         new_members_db.name, target_path.name)
            self._copy_db(new_members_db, tmp_name)
        else:
            logger.debug("Merging new members: concat %s + %s → %s",
                         self._members_db.name, new_members_db.name, target_path.name)
            self._concat_dbs(self._members_db, new_members_db, tmp_name)

        self._finalize_db_rename(tmp_name, target_path)

        previous_db = self._members_db
        self._members_db = target_path
        self._members_generation += 1

        if previous_db is not None:
            self._delete_db(previous_db)
        logger.debug("Members DB now at generation %d: %s",
                     self._members_generation, self._members_db.name)
 
 
    def _create_db(self, fasta_path: Union[str, Path], db_path: Path) -> Path:
        self._run_mmseqs(["createdb", str(fasta_path), str(db_path)])
        return db_path
 
 
    def _concat_dbs(self, db_a: Path, db_b: Path, out_db: Path) -> None:
        self._run_mmseqs(["concatdbs", str(db_a), str(db_b), str(out_db)])


    def _run_search(self, query_db: Path, target_db: Path, tmp_dir: Path, n_targets: int) -> pd.DataFrame:
        tmp_dir.mkdir(parents=True, exist_ok=True)
        result_db = tmp_dir / "result_db"
        result_tsv = tmp_dir / "result.tsv" # mmseqs search output in tabular format
        search_tmp = tmp_dir / "mmseqs_tmp"
        search_tmp.mkdir(parents=True, exist_ok=True)

        logger.debug("Running mmseqs search: query=%s, target=%s, sensitivity=%s, evalue=%s",
                     query_db.name, target_db.name, self.sensitivity, self.max_evalue)
        self._run_mmseqs([
            "search", str(query_db), str(target_db), str(result_db), str(search_tmp),
            "--min-seq-id", "0",
            "-e", self.max_evalue,
            "-s", str(self.sensitivity),
            "--alignment-mode", "3",
            "--max-seqs", str(n_targets)
        ])
        self._run_mmseqs([
            "convertalis", str(query_db), str(target_db), str(result_db), str(result_tsv),
            "--format-output", "query,target,pident",
        ])

        if not result_tsv.exists() or result_tsv.stat().st_size == 0:
            logger.debug("Search produced no alignments (empty result.tsv)")
            return pd.DataFrame(columns=self.MMSEQS_OUTPUT_COLUMNS)
        df = pd.read_csv(result_tsv, sep="\t", names=self.MMSEQS_OUTPUT_COLUMNS)
        logger.debug("Search returned %d alignments", len(df))
        return df
 
 
    def _run_mmseqs(self, args: List[str]) -> None:
        cmd = [self.mmseqs_binary] + args
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"mmseqs {args[0]} failed:\n{proc.stderr}")
 
 
    @staticmethod
    def _copy_db(src_db: Path, dst_db: Path) -> None:
        for src_file in src_db.parent.glob(f"{src_db.name}.*"):
            suffix = src_file.name[len(src_db.name):]
            shutil.copy2(src_file, dst_db.parent / (dst_db.name + suffix))
 
 
    @staticmethod
    def _delete_db(db_path: Path) -> None:
        for f in db_path.parent.glob(f"{db_path.name}.*"):
            f.unlink(missing_ok=True)
 
 
    @staticmethod
    def _finalize_db_rename(tmp_db: Path, final_db: Path) -> None:
        for tmp_file in tmp_db.parent.glob(f"{tmp_db.name}.*"):
            suffix = tmp_file.name[len(tmp_db.name):]
            os.replace(tmp_file, final_db.parent / (final_db.name + suffix))
 
 
    def _get_db_ids(self, db_path: Path) -> List[str]:
        """Reads the .lookup file mmseqs createdb/concatdbs produces (id, name, ...)."""
        lookup_path = db_path.parent / f"{db_path.name}.lookup"
        ids = []
        if lookup_path.exists():
            with open(lookup_path) as f:
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 2:
                        ids.append(parts[1])
        return ids
    
 
    def _search_all_pairs(self, query_db: Path, target_db: Path,
                           query_ids: List[str], target_ids: List[str],
                           tmp_dir: Path) -> pd.DataFrame:
        hits_df = self._run_search(query_db, target_db, tmp_dir, len(target_ids))
        return self._fill_missing_pairs(hits_df, query_ids, target_ids)


    def _fill_missing_pairs(self, hits_df: pd.DataFrame,
                             query_ids: List[str], target_ids: List[str]) -> pd.DataFrame:
        """
        Builds the full query x target cartesian product, left-joins MMseqs2
        hits onto it, and fills any pair with no reported alignment with
        pident = 0.0, flagged as detected = False.
        """
        all_pairs = pd.DataFrame(list(product(query_ids, target_ids)),
                                  columns=["query", "target"])

        if hits_df.empty:
            best_hits = hits_df
        else:
            best_hits = hits_df.groupby(["query", "target"], as_index=False)["pident"].max()

        merged = all_pairs.merge(best_hits, on=["query", "target"], how="left")
        merged["detected"] = merged["pident"].notna()
        merged["pident"] = merged["pident"].fillna(0.0)

        n_undetected = (~merged["detected"]).sum()
        if n_undetected > 0:
            logger.debug("fill_missing_pairs: %d / %d pairs below detection threshold",
                         n_undetected, len(merged))
        return merged.rename(columns={"query": "seed_id", "target": "hit_id"})
    
 
    @staticmethod
    def _get_fasta_ids(fasta_path: Union[str, Path]) -> List[str]:
        """First whitespace-delimited token of each header, matching MMseqs2's own ID parsing."""
        ids = []
        with open(fasta_path) as f:
            for line in f:
                if line.startswith(">"):
                    ids.append(line[1:].strip().split()[0])
        return ids
 
 
    @staticmethod
    def _fasta_has_sequences(fasta_path: Union[str, Path]) -> bool:
        fasta_path = Path(fasta_path)
        if not fasta_path.exists() or fasta_path.stat().st_size == 0:
            return False
        with open(fasta_path) as f:
            return any(line.startswith(">") for line in f)
 
 
    @staticmethod
    def _append_to_csv(df: pd.DataFrame, output_csv_path: str) -> None:
        output_path = Path(output_csv_path)
        write_header = not output_path.exists()
        df.to_csv(output_path, mode="a", header=write_header, index=False)
        