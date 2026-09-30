"""Unit tests for scripts/SequenceIdentityCalculator.py."""

import os
from pathlib import Path

import pandas as pd
import pytest

from src.scripts.SequenceIdentityCalculator import SequenceIdentityCalculator


def _make_fasta(work_dir: Path, name: str, headers_and_seqs: list[tuple[str, str]]) -> str:
    path = work_dir / name
    with open(path, "w") as f:
        for header, seq in headers_and_seqs:
            f.write(f">{header}\n{seq}\n")
    return str(path)


def _make_lookup(db_dir: Path, db_name: str, ids: list[str]) -> Path:
    db_dir.mkdir(parents=True, exist_ok=True)
    path = db_dir / f"{db_name}.lookup"
    with open(path, "w") as f:
        for i, seq_id in enumerate(ids):
            f.write(f"{i}\t{seq_id}\t{seq_id}\n")
    return path


def _make_db_dir(db_dir: Path, db_name: str) -> Path:
    """Creates a directory representing an MMseqs2 database (named after its prefix)."""
    db_path = db_dir / db_name
    db_path.mkdir(parents=True, exist_ok=True)
    return db_path


def _make_db_files(db_dir: Path, db_name: str, suffixes: list[str]) -> None:
    db_dir.mkdir(parents=True, exist_ok=True)
    for suffix in suffixes:
        (db_dir / f"{db_name}{suffix}").touch()


@pytest.fixture
def calc(tmp_path, mocker):
    mocker.patch("shutil.which", return_value="/usr/bin/mmseqs")
    work_dir = tmp_path / "seqid_work"
    return SequenceIdentityCalculator(work_dir=str(work_dir))


class TestGetFastaIds:
    def test_extracts_first_token(self, tmp_path):
        fasta = _make_fasta(tmp_path, "test.fasta", [
            ("sp|P12345|PROT_HUMAN/1-100", "MK"),
            ("tr|Q67890|PROT_MOUSE/5-50", "AIL"),
            ("UniRef90_A0A3P4APY8", "GGG"),
        ])
        result = SequenceIdentityCalculator._get_fasta_ids(fasta)
        assert result == ["sp|P12345|PROT_HUMAN/1-100", "tr|Q67890|PROT_MOUSE/5-50", "UniRef90_A0A3P4APY8"]

    def test_empty_file_returns_empty_list(self, tmp_path):
        fasta = tmp_path / "empty.fasta"
        fasta.write_text("")
        result = SequenceIdentityCalculator._get_fasta_ids(str(fasta))
        assert result == []

    def test_no_headers_returns_empty_list(self, tmp_path):
        fasta = tmp_path / "no_headers.fasta"
        fasta.write_text("MK\nAIL\n")
        result = SequenceIdentityCalculator._get_fasta_ids(str(fasta))
        assert result == []

    def test_splits_on_whitespace(self, tmp_path):
        fasta = _make_fasta(tmp_path, "test.fasta", [
            ("id_with_spaces extra_field/1-100", "MK"),
        ])
        result = SequenceIdentityCalculator._get_fasta_ids(str(fasta))
        assert result == ["id_with_spaces"]


class TestFastaHasSequences:
    def test_has_headers_returns_true(self, tmp_path):
        fasta = _make_fasta(tmp_path, "test.fasta", [("P12345", "MK")])
        assert SequenceIdentityCalculator._fasta_has_sequences(str(fasta))

    def test_empty_file_returns_false(self, tmp_path):
        fasta = tmp_path / "empty.fasta"
        fasta.write_text("")
        assert not SequenceIdentityCalculator._fasta_has_sequences(str(fasta))

    def test_no_headers_returns_false(self, tmp_path):
        fasta = tmp_path / "noseq.fasta"
        fasta.write_text("MK\nAIL\n")
        assert not SequenceIdentityCalculator._fasta_has_sequences(str(fasta))

    def test_nonexistent_file_returns_false(self, tmp_path):
        assert not SequenceIdentityCalculator._fasta_has_sequences(str(tmp_path / "nonexistent.fasta"))


class TestGetDbIds:
    def test_reads_lookup_file(self, calc):
        _make_lookup(calc.work_dir, "my_db", ["P12345", "Q67890", "A11111"])
        result = calc._get_db_ids(calc.work_dir / "my_db")
        assert result == ["P12345", "Q67890", "A11111"]

    def test_missing_lookup_returns_empty(self, calc):
        result = calc._get_db_ids(calc.work_dir / "nonexistent_db")
        assert result == []

    def test_partial_lines_skipped(self, calc):
        lookup = calc.work_dir / "my_db.lookup"
        lookup.write_text("0\tP12345\tblah\n1\t\n2\tQ67890\tblah\n")
        result = calc._get_db_ids(calc.work_dir / "my_db")
        assert result[0] == "P12345"
        assert result[-1] == "Q67890"


class TestFillMissingPairs:
    def test_all_pairs_detected(self):
        hits = pd.DataFrame({
            "query": ["A", "A", "B", "B"],
            "target": ["X", "Y", "X", "Y"],
            "pident": [0.9, 0.8, 0.7, 0.3],
        })
        result = SequenceIdentityCalculator._fill_missing_pairs(
            None, hits, ["A", "B"], ["X", "Y"]
        )
        assert len(result) == 4
        assert result["detected"].all()

    def test_missing_pairs_filled_with_zero(self):
        hits = pd.DataFrame({
            "query": ["A"],
            "target": ["X"],
            "pident": [0.9],
        })
        result = SequenceIdentityCalculator._fill_missing_pairs(
            None, hits, ["A", "B"], ["X", "Y"]
        )
        assert len(result) == 4
        detected_mask = (result["seed_id"] == "A") & (result["hit_id"] == "X")
        assert result.loc[detected_mask, "detected"].iloc[0]
        assert result.loc[detected_mask, "pident"].iloc[0] == 0.9
        undetected = result[~result["detected"]]
        assert len(undetected) == 3
        assert (undetected["pident"] == 0.0).all()

    def test_empty_hits_all_undetected(self):
        hits = pd.DataFrame(columns=["query", "target", "pident"])
        result = SequenceIdentityCalculator._fill_missing_pairs(
            None, hits, ["A"], ["X", "Y"]
        )
        assert len(result) == 2
        assert not result["detected"].any()
        assert (result["pident"] == 0.0).all()

    def test_max_pident_per_duplicate_pair(self):
        hits = pd.DataFrame({
            "query": ["A", "A"],
            "target": ["X", "X"],
            "pident": [0.5, 0.9],
        })
        result = SequenceIdentityCalculator._fill_missing_pairs(
            None, hits, ["A"], ["X"]
        )
        assert len(result) == 1
        assert result["pident"].iloc[0] == 0.9


class TestDiscoverLatestMembersDb:
    def test_empty_dir_returns_none(self, calc):
        gen, db = calc._discover_latest_members_db()
        assert gen == -1
        assert db is None

    def test_finds_highest_generation(self, calc):
        _make_db_dir(calc.work_dir, "members_db_gen1")
        _make_db_dir(calc.work_dir, "members_db_gen3")
        _make_db_dir(calc.work_dir, "members_db_gen2")
        gen, db = calc._discover_latest_members_db()
        assert gen == 3
        assert db == calc.work_dir / "members_db_gen3"

    def test_ignores_non_matching_dirs(self, calc):
        _make_db_dir(calc.work_dir, "members_db_gen1")
        _make_db_dir(calc.work_dir, "other_dir")
        gen, db = calc._discover_latest_members_db()
        assert gen == 1

    def test_prefix_no_false_match(self, calc):
        _make_db_dir(calc.work_dir, "members_db_gen5")
        _make_db_dir(calc.work_dir, "members_db_gen50")
        gen, db = calc._discover_latest_members_db()
        assert gen == 50
        assert db == calc.work_dir / "members_db_gen50"


class TestNextMembersDbPath:
    def test_increments_generation(self, calc):
        calc._members_generation = 3
        path = calc._next_members_db_path()
        assert path == calc.work_dir / "members_db_gen4"

    def test_first_generation(self, calc):
        calc._members_generation = -1
        path = calc._next_members_db_path()
        assert path == calc.work_dir / "members_db_gen0"


class TestCopyDb:
    def test_copies_all_suffix_files(self, tmp_path):
        _make_db_files(tmp_path, "src_db", [".index", ".lookup", ".dbtype", ".source", ".h"])
        SequenceIdentityCalculator._copy_db(
            tmp_path / "src_db", tmp_path / "dst_db"
        )
        for suffix in [".index", ".lookup", ".dbtype", ".source", ".h"]:
            assert (tmp_path / f"dst_db{suffix}").exists()

    def test_does_not_match_other_prefixes(self, tmp_path):
        _make_db_files(tmp_path, "src_db", [".index"])
        _make_db_files(tmp_path, "src_db_other", [".txt"])
        SequenceIdentityCalculator._copy_db(
            tmp_path / "src_db", tmp_path / "dst_db"
        )
        assert (tmp_path / "dst_db.index").exists()
        assert not (tmp_path / "dst_db_other.txt").exists()


class TestDeleteDb:
    def test_deletes_all_suffix_files(self, tmp_path):
        _make_db_files(tmp_path, "old_db", [".index", ".lookup", ".dbtype"])
        SequenceIdentityCalculator._delete_db(tmp_path / "old_db")
        for suffix in [".index", ".lookup", ".dbtype"]:
            assert not (tmp_path / f"old_db{suffix}").exists()

    def test_no_files_is_noop(self, tmp_path):
        SequenceIdentityCalculator._delete_db(tmp_path / "nonexistent")

    def test_does_not_match_other_prefixes(self, tmp_path):
        _make_db_files(tmp_path, "old_db", [".index"])
        _make_db_files(tmp_path, "old_db_keep", [".index"])
        SequenceIdentityCalculator._delete_db(tmp_path / "old_db")
        assert not (tmp_path / "old_db.index").exists()
        assert (tmp_path / "old_db_keep.index").exists()


class TestFinalizeDbRename:
    def test_renames_all_suffix_files(self, tmp_path):
        _make_db_files(tmp_path, ".tmp_abc_members_db_gen5", [".index", ".lookup"])
        SequenceIdentityCalculator._finalize_db_rename(
            tmp_path / ".tmp_abc_members_db_gen5",
            tmp_path / "members_db_gen5",
        )
        assert not (tmp_path / ".tmp_abc_members_db_gen5.index").exists()
        assert (tmp_path / "members_db_gen5.index").exists()
        assert (tmp_path / "members_db_gen5.lookup").exists()


class TestMergeNewMembers:
    def test_first_generation_copies(self, calc, mocker):
        mocker.patch.object(calc, "_copy_db")
        mocker.patch.object(calc, "_concat_dbs")
        mocker.patch.object(calc, "_finalize_db_rename")
        mocker.patch.object(calc, "_delete_db")

        new_db = calc.work_dir / "new_db"
        calc._members_db = None
        calc._members_generation = -1
        calc._merge_new_members(new_db)

        calc._copy_db.assert_called_once()
        calc._concat_dbs.assert_not_called()
        calc._finalize_db_rename.assert_called_once()
        calc._delete_db.assert_not_called()
        assert calc._members_generation == 0

    def test_subsequent_generation_concats(self, calc, mocker):
        mocker.patch.object(calc, "_copy_db")
        mocker.patch.object(calc, "_concat_dbs")
        mocker.patch.object(calc, "_finalize_db_rename")
        mocker.patch.object(calc, "_delete_db")

        existing_db = calc.work_dir / "members_db_gen2"
        calc._members_db = existing_db
        calc._members_generation = 2
        new_db = calc.work_dir / "new_db"
        calc._merge_new_members(new_db)

        calc._copy_db.assert_not_called()
        calc._concat_dbs.assert_called_once_with(existing_db, new_db, mocker.ANY)
        calc._finalize_db_rename.assert_called_once()
        calc._delete_db.assert_called_once_with(existing_db)
        assert calc._members_generation == 3


class TestComputeSequenceIdentity:
    def test_no_inputs_returns_empty_df(self, calc, mocker, tmp_path):
        mocker.patch.object(calc, "_run_mmseqs")
        mocker.patch.object(calc, "_merge_new_members")

        empty_fasta = tmp_path / "empty.fasta"
        empty_fasta.write_text("")
        csv_path = str(tmp_path / "out.csv")

        result = calc.compute_sequence_identity(
            output_csv_path=csv_path,
            current_seed_fasta_path=str(empty_fasta),
            new_members_fasta_path=str(empty_fasta),
            previous_used_seeds_fasta_path=str(empty_fasta),
        )
        assert len(result) == 0
        assert list(result.columns) == ["seed_id", "hit_id", "pident", "detected"]
        calc._merge_new_members.assert_not_called()

    def test_seed_vs_persisted_db(self, calc, mocker, tmp_path):
        mocker.patch.object(calc, "_run_mmseqs")

        _make_db_dir(calc.work_dir, "members_db_gen0")
        _make_lookup(calc.work_dir, "members_db_gen0", ["P12345", "Q67890"])
        calc._members_db = calc.work_dir / "members_db_gen0"
        calc._members_generation = 0

        seed_fasta = _make_fasta(tmp_path, "seed.fasta", [("sp|SEED|blah", "MK")])
        empty_fasta = tmp_path / "empty.fasta"
        empty_fasta.write_text("")

        mocker.patch.object(calc, "_merge_new_members")

        fake_results = pd.DataFrame({
            "query": ["sp|SEED|blah", "sp|SEED|blah"],
            "target": ["P12345", "Q67890"],
            "pident": [0.8, 0.6],
        })
        mocker.patch.object(calc, "_run_search", return_value=fake_results)

        csv_path = str(tmp_path / "out.csv")
        result = calc.compute_sequence_identity(
            output_csv_path=csv_path,
            current_seed_fasta_path=seed_fasta,
            new_members_fasta_path=str(empty_fasta),
            previous_used_seeds_fasta_path=str(empty_fasta),
        )

        assert len(result) == 2
        assert result["detected"].all()
        assert os.path.exists(csv_path)

    def test_seed_and_prev_seeds_vs_new_members(self, calc, mocker, tmp_path):
        mocker.patch.object(calc, "_run_mmseqs")
        mocker.patch.object(calc, "_copy_db")
        mocker.patch.object(calc, "_concat_dbs")
        mocker.patch.object(calc, "_finalize_db_rename")
        mocker.patch.object(calc, "_delete_db")

        seed_fasta = _make_fasta(tmp_path, "seed.fasta", [("SEED_ID", "MK")])
        prev_fasta = _make_fasta(tmp_path, "prev.fasta", [
            ("PREV1", "MK"), ("PREV2", "AIL")
        ])
        new_fasta = _make_fasta(tmp_path, "new.fasta", [
            ("NEW1", "GGG"), ("NEW2", "PPP")
        ])

        calc._members_db = None
        calc._members_generation = -1

        call_count = [0]

        def fake_search(self, query_db, target_db, tmp_dir):
            call_count[0] += 1
            call_num = call_count[0]
            if call_num == 1:
                return pd.DataFrame({
                    "query": ["SEED_ID", "SEED_ID"],
                    "target": ["NEW1", "NEW2"],
                    "pident": [0.5, 0.4],
                })
            if call_num == 2:
                return pd.DataFrame({
                    "query": ["PREV1", "PREV1", "PREV2"],
                    "target": ["NEW1", "NEW2", "NEW2"],
                    "pident": [0.3, 0.2, 0.1],
                })
            return pd.DataFrame(columns=["query", "target", "pident"])

        mocker.patch.object(calc, "_run_search", side_effect=fake_search)

        csv_path = str(tmp_path / "out.csv")
        result = calc.compute_sequence_identity(
            output_csv_path=csv_path,
            current_seed_fasta_path=seed_fasta,
            new_members_fasta_path=new_fasta,
            previous_used_seeds_fasta_path=prev_fasta,
        )

        assert len(result) == 6
        assert result["detected"].sum() == 5

    def test_appends_to_existing_csv(self, calc, mocker, tmp_path):
        mocker.patch.object(calc, "_run_mmseqs")
        mocker.patch.object(calc, "_merge_new_members")

        empty_fasta = tmp_path / "empty.fasta"
        empty_fasta.write_text("")

        csv_path = str(tmp_path / "out.csv")
        pd.DataFrame({
            "seed_id": ["OLD_SEED"], "hit_id": ["OLD_HIT"],
            "pident": [0.9], "detected": [True],
        }).to_csv(csv_path, index=False)

        result = calc.compute_sequence_identity(
            output_csv_path=csv_path,
            current_seed_fasta_path=str(empty_fasta),
            new_members_fasta_path=str(empty_fasta),
            previous_used_seeds_fasta_path=str(empty_fasta),
        )

        df = pd.read_csv(csv_path)
        assert len(df) == 1
        assert df["seed_id"].iloc[0] == "OLD_SEED"


class TestWithRealMmmseqsData:
    """Tests that exercise real MMseqs2 output files from test_data/output_mmseq2."""

    def test_get_fasta_ids_real_data(self, mmseqs_fasta_1):
        ids = SequenceIdentityCalculator._get_fasta_ids(mmseqs_fasta_1)
        assert len(ids) == 5
        assert ids[0] == "tr|A0AAU9EQU2|A0AAU9EQU2_9BACT"
        assert "tr|A0A369A485|A0A369A485_9FLAO" in ids

    def test_fasta_has_sequences_real_data(self, mmseqs_fasta_1, mmseqs_fasta_2):
        assert SequenceIdentityCalculator._fasta_has_sequences(mmseqs_fasta_1)
        assert SequenceIdentityCalculator._fasta_has_sequences(mmseqs_fasta_2)

    def test_get_db_ids_real_lookup_db1(self, calc, mmseqs_db1):
        ids = calc._get_db_ids(Path(mmseqs_db1))
        assert len(ids) == 5
        assert "A0AAU9EQU2" in ids

    def test_get_db_ids_real_lookup_db2(self, calc, mmseqs_db2):
        ids = calc._get_db_ids(Path(mmseqs_db2))
        assert len(ids) == 1
        assert ids[0] == "A0AAU9EQU2"

    def test_get_db_ids_real_lookup_db3(self, calc, mmseqs_db3):
        ids = calc._get_db_ids(Path(mmseqs_db3))
        assert len(ids) == 6

    def test_parse_real_mmseqs_result_tsv(self, mmseqs_result_tsv, calc):
        """Parse a real MMseqs2 convertalis TSV output and run through _fill_missing_pairs."""
        from src.scripts.SequenceIdentityCalculator import SequenceIdentityCalculator as SIC

        df = pd.read_csv(
            mmseqs_result_tsv, sep="\t",
            names=SIC.MMSEQS_OUTPUT_COLUMNS,
        )
        assert len(df) == 3
        assert list(df.columns) == ["query", "target", "pident"]
        assert df.loc[df["query"] == "A0AAU9EQU2", "pident"].iloc[0] == 100.0

        result = calc._fill_missing_pairs(
            df, ["A0AAU9EQU2", "A0A401XKQ8", "A0A369A485"],
            ["A0AAU9EQU2"],
        )
        assert len(result) == 3
        assert result["detected"].all()
        assert result.loc[result["seed_id"] == "A0AAU9EQU2", "pident"].iloc[0] == 100.0