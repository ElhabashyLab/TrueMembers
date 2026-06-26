"""Unit tests for scripts/hmm_domain_filter.py — HMM parsing and domtblout parsing."""

import pandas as pd
import pytest

from True_Members.scripts.hmm_domain_filter import (
    get_hmms_with_ga_cutoffs,
    parse_ga,
    parse_hmmsearch_results,
    parse_name,
    run_hmm_domain_filter,
)
from True_Members.scripts.types import JobStats


class TestParseName:
    def test_extracts_name(self):
        assert parse_name("NAME  HisKA") == "HisKA"

    def test_extracts_name_with_numbers(self):
        assert parse_name("NAME  3dh4A01-i2") == "3dh4A01-i2"

    def test_single_word(self):
        assert parse_name("NAME") is None

    def test_empty_string(self):
        assert parse_name("") is None


class TestParseGA:
    def test_extracts_ga_score(self):
        assert parse_ga("GA  22.8  22.8;") == 22.8

    def test_larger_score(self):
        assert parse_ga("GA  105.0  105.0;") == 105.0

    def test_no_semicolon(self):
        assert parse_ga("GA  22.8  22.8") == 22.8

    def test_insufficient_parts(self):
        assert parse_ga("GA") is None
        assert parse_ga("GA  22.8") is None


class TestGetHmmsWithGaCutoffs:
    def test_extracts_all_hmms(self, hmms_dir):
        cutoffs = get_hmms_with_ga_cutoffs(hmms_dir)
        assert len(cutoffs) >= 2
        assert "HisKA" in cutoffs
        assert "HATPase_c" in cutoffs

    def test_empty_directory(self, tmp_path):
        cutoffs = get_hmms_with_ga_cutoffs(str(tmp_path))
        assert cutoffs == {}

    def test_nonexistent_directory(self):
        cutoffs = get_hmms_with_ga_cutoffs("/nonexistent/path/xyz")
        assert cutoffs == {}


def _make_hit(seq_id, binscore_val):
    hit = type("FakeHit", (), {
        "id": seq_id,
        "hsps": [type("FakeHSP", (), {"bitscore": binscore_val, "evalue": 1e-30})()],
    })()
    return hit


def _make_qresult(hmm_name, hit_id, binscore_vals):
    qresult = type("FakeQResult", (), {
        "id": hmm_name,
        "hits": [_make_hit(hid, b) for hid, b in zip(hit_id, binscore_vals)],
    })()
    return qresult


class TestParseHmmsearchResults:
    def test_builds_pivot_table(self, mocker):
        fake_results = [
            _make_qresult("HisKA", ["Q5SK37/507-569", "A0A3P4APY8/507-569", "P12345/508-565"], [100.0, 85.0, 70.0]),
            _make_qresult("HATPase_c", ["Q5SK37/507-569", "A0A3P4APY8/507-569", "P12345/508-565"], [85.0, 67.0, 54.0]),
        ]
        mocker.patch("True_Members.scripts.hmm_domain_filter.SearchIO.parse", return_value=fake_results)

        hmms_with_cutoffs = {"HisKA": 25.0, "HATPase_c": 25.0}
        result = parse_hmmsearch_results("dummy_path", hmms_with_cutoffs)

        assert isinstance(result, pd.DataFrame)
        assert "id" in result.columns
        assert "HisKA" in result.columns
        assert "HATPase_c" in result.columns
        assert len(result) == 3

    def test_all_domains_present_for_first_id(self, mocker):
        fake_results = [
            _make_qresult("HisKA", ["Q5SK37/507-569", "A0A3P4APY8/507-569"], [100.0, 85.0]),
            _make_qresult("HATPase_c", ["Q5SK37/507-569"], [85.0]),
        ]
        mocker.patch("True_Members.scripts.hmm_domain_filter.SearchIO.parse", return_value=fake_results)

        hmms_with_cutoffs = {"HisKA": 25.0, "HATPase_c": 25.0}
        result = parse_hmmsearch_results("dummy_path", hmms_with_cutoffs)

        q5sk37 = result[result["id"] == "Q5SK37"].iloc[0]
        assert q5sk37["HisKA"] == True
        assert q5sk37["HATPase_c"] == True

    def test_high_cutoff_rejects_all(self, mocker):
        fake_results = [
            _make_qresult("HisKA", ["P12345/508-565"], [70.0]),
            _make_qresult("HATPase_c", ["P12345/508-565"], [54.0]),
        ]
        mocker.patch("True_Members.scripts.hmm_domain_filter.SearchIO.parse", return_value=fake_results)

        hmms_with_cutoffs = {"HisKA": 500.0, "HATPase_c": 500.0}
        result = parse_hmmsearch_results("dummy_path", hmms_with_cutoffs)

        upi = result[result["id"] == "P12345"].iloc[0]
        assert upi["HisKA"] == False
        assert upi["HATPase_c"] == False

    def test_no_cutoff_falls_back_to_evalue(self, mocker):
        fake_results = [
            _make_qresult("HisKA", ["Q5SK37/507-569"], [100.0]),
        ]
        mocker.patch("True_Members.scripts.hmm_domain_filter.SearchIO.parse", return_value=fake_results)

        hmms_with_cutoffs = {"HisKA": None, "HATPase_c": None}
        result = parse_hmmsearch_results("dummy_path", hmms_with_cutoffs)

        q5sk37 = result[result["id"] == "Q5SK37"].iloc[0]
        assert q5sk37["HisKA"] == True

    def test_empty_searchio_result_returns_empty_df(self, mocker):
        mocker.patch("True_Members.scripts.hmm_domain_filter.SearchIO.parse", return_value=[])

        hmms_with_cutoffs = {"HisKA": 25.0}
        result = parse_hmmsearch_results("dummy_path", hmms_with_cutoffs)

        assert isinstance(result, pd.DataFrame)
        assert result.empty

    def test_extra_hmm_in_cutoffs_filled_with_false(self, mocker):
        fake_results = [
            _make_qresult("HisKA", ["Q5SK37/507-569", "A0A3P4APY8/507-569", "P12345/508-565"], [100.0, 85.0, 70.0]),
        ]
        mocker.patch("True_Members.scripts.hmm_domain_filter.SearchIO.parse", return_value=fake_results)

        hmms_with_cutoffs = {"HisKA": 25.0, "ExtraDomain": 30.0}
        result = parse_hmmsearch_results("dummy_path", hmms_with_cutoffs)

        assert "ExtraDomain" in result.columns
        assert result["ExtraDomain"].tolist() == [False, False, False]


class TestRunHmmDomainFilter:
    def test_empty_seed_containing_sequences_returns_empty_df(self, mocker, hmms_dir, tmp_path):
        mock_records = []
        mocker.patch(
            "True_Members.scripts.hmm_domain_filter.get_records_from_fasta",
            return_value=mock_records,
        )
        mock_press = mocker.patch(
            "True_Members.scripts.hmm_domain_filter.press_hmms",
            return_value="/fake/domain_validator.hmm",
        )
        mock_subprocess = mocker.patch("subprocess.run")

        stats = JobStats()
        empty_seed_containing = pd.DataFrame(columns=["id", "qtmscore"])

        result = run_hmm_domain_filter(
            stats=stats,
            cpus=1,
            uids_with_region_fasta_path="/fake/uids.fasta",
            seed_containing_sequences=empty_seed_containing,
            HMMs_DIR=hmms_dir,
            hmm_domain_filtering_dir=str(tmp_path),
            root_dir=str(tmp_path),
            TM_SCORE_CUTOFF=0.6,
            hmms_with_ga_cutoffs={"HisKA": 25.0},
        )

        assert isinstance(result, pd.DataFrame)
        assert result.empty or list(result.columns) == ["id"]
        assert stats.hmm_search_results == 0
        mock_subprocess.assert_not_called()

    def test_empty_seed_containing_no_ids_match_records(self, mocker, hmms_dir, tmp_path):
        mock_record = type("FakeRecord", (), {"id": "P12345/10-20"})()
        mock_records = [mock_record]
        mocker.patch(
            "True_Members.scripts.hmm_domain_filter.get_records_from_fasta",
            return_value=mock_records,
        )
        mocker.patch(
            "True_Members.scripts.hmm_domain_filter.press_hmms",
            return_value="/fake/domain_validator.hmm",
        )
        mock_subprocess = mocker.patch("subprocess.run")

        stats = JobStats()
        seed_containing = pd.DataFrame({"id": ["Q99999"], "qtmscore": [0.85]})

        result = run_hmm_domain_filter(
            stats=stats,
            cpus=1,
            uids_with_region_fasta_path="/fake/uids.fasta",
            seed_containing_sequences=seed_containing,
            HMMs_DIR=hmms_dir,
            hmm_domain_filtering_dir=str(tmp_path),
            root_dir=str(tmp_path),
            TM_SCORE_CUTOFF=0.6,
            hmms_with_ga_cutoffs={"HisKA": 25.0},
        )

        assert isinstance(result, pd.DataFrame)
        assert stats.hmm_search_results == 0
        mock_subprocess.assert_not_called()