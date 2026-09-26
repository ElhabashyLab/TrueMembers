"""Unit tests for scripts/utils.py — ID parsing, region parsing, FASTA operations."""

import os

import yaml

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from True_Members.scripts.utils import (
    _build_batch_query,
    _build_uniprot_url,
    _chunked,
    _collect_cached_records,
    _log_interval,
    copy_file,
    extract_uid_from_recordID,
    extract_region_from_recordID,
    filter_fasta_file,
    get_records_from_fasta,
    read_yaml_config,
    write_fasta_from_records,
)


class TestExtractUidFromRecordID:
    def test_standard_uniprot_accession(self):
        result = extract_uid_from_recordID("tr|A0A1A9EVS2|A0A1A9EVS2_9GAMM/553-621")
        assert result == "A0A1A9EVS2"

    def test_bare_accession_with_region(self):
        result = extract_uid_from_recordID("Q5SK37/507-569")
        assert result == "Q5SK37"

    def test_uniref90_cluster_id(self):
        result = extract_uid_from_recordID("UniRef90_A0A3P4APY8/507-569")
        assert result == "A0A3P4APY8"

    def test_upi_ids_return_none(self):
        result = extract_uid_from_recordID("UniRef90_UPI000542879D/506-563")
        assert result is None

    def test_underscore_suffix(self):
        result = extract_uid_from_recordID("tr|A0A6S5RMM9|A0A6S5RMM9_9GAMM/544-610")
        assert result == "A0A6S5RMM9"

    def test_record_with_gn_field(self):
        result = extract_uid_from_recordID("sp|P12345|PROT_HUMAN/1-100")
        assert result == "P12345"

    def test_no_region(self):
        result = extract_uid_from_recordID("tr|A0A1A9EVS2|A0A1A9EVS2_9GAMM")
        assert result == "A0A1A9EVS2"

    def test_no_match_returns_none(self):
        result = extract_uid_from_recordID("not_a_valid_id")
        assert result is None

    def test_empty_string(self):
        result = extract_uid_from_recordID("")
        assert result is None


class TestExtractRegionFromRecordID:
    def test_standard_region(self):
        start, end = extract_region_from_recordID("tr|A0A1A9EVS2|A0A1A9EVS2_9GAMM/553-621")
        assert start == 553
        assert end == 621

    def test_single_digit_region(self):
        start, end = extract_region_from_recordID("Q5SK37/5-9")
        assert start == 5
        assert end == 9

    def test_large_region(self):
        start, end = extract_region_from_recordID("A0A3P4APY8/1-9999")
        assert start == 1
        assert end == 9999

    def test_no_region_returns_none(self):
        start, end = extract_region_from_recordID("A0A3P4APY8")
        assert start is None
        assert end is None

    def test_empty_string(self):
        start, end = extract_region_from_recordID("")
        assert start is None
        assert end is None


class TestFilterFastaFile:
    def test_filters_correct_ids(self, tmp_path):
        records = [
            SeqRecord(Seq("MKAIL"), id="tr|P12345|PROT_HUMAN/1-100"),
            SeqRecord(Seq("MKAIL"), id="sp|Q67890|PROT_MOUSE/1-50"),
            SeqRecord(Seq("MKAIL"), id="tr|A11111|PROT_FISH/10-60"),
        ]
        fasta_path = tmp_path / "test.fasta"
        write_fasta_from_records(records, str(fasta_path))

        out_path = str(tmp_path / "filtered.fasta")
        result = filter_fasta_file(str(fasta_path), ["P12345", "Q67890"], out_path)
        filtered = get_records_from_fasta(result)
        assert len(filtered) == 2

    def test_empty_id_list(self, tmp_path):
        records = [
            SeqRecord(Seq("MKAIL"), id="P12345/1-100"),
        ]
        fasta_path = tmp_path / "test.fasta"
        write_fasta_from_records(records, str(fasta_path))
        out_path = str(tmp_path / "filtered.fasta")
        result = filter_fasta_file(str(fasta_path), [], out_path)
        filtered = get_records_from_fasta(result)
        assert len(filtered) == 0

    def test_missing_id_warns(self, tmp_path):
        records = [
            SeqRecord(Seq("MKAIL"), id="P12345/1-100"),
        ]
        fasta_path = tmp_path / "test.fasta"
        write_fasta_from_records(records, str(fasta_path))
        out_path = str(tmp_path / "filtered.fasta")
        result = filter_fasta_file(str(fasta_path), ["NONEXISTENT"], out_path)
        filtered = get_records_from_fasta(result)
        assert len(filtered) == 0


class TestGetRecordsFromFasta:
    def test_parses_single_record(self, tmp_path):
        records = [SeqRecord(Seq("MK"), id="P12345")]
        fasta_path = tmp_path / "test.fasta"
        write_fasta_from_records(records, str(fasta_path))
        result = get_records_from_fasta(str(fasta_path))
        assert len(result) == 1
        assert result[0].id == "P12345"

    def test_parses_multiple_records(self, tmp_path):
        records = [
            SeqRecord(Seq("MK"), id="P12345"),
            SeqRecord(Seq("AIL"), id="Q67890"),
            SeqRecord(Seq("GGG"), id="A11111"),
        ]
        fasta_path = tmp_path / "test.fasta"
        write_fasta_from_records(records, str(fasta_path))
        result = get_records_from_fasta(str(fasta_path))
        assert len(result) == 3

    def test_empty_file_returns_empty_list(self, tmp_path):
        fasta_path = tmp_path / "empty.fasta"
        fasta_path.write_text("")
        result = get_records_from_fasta(str(fasta_path))
        assert result == []


class TestWriteFastaFromRecords:
    def test_writes_and_reads_back(self, tmp_path):
        records = [SeqRecord(Seq("MK"), id="P12345")]
        out_path = str(tmp_path / "out.fasta")
        result = write_fasta_from_records(records, out_path)
        assert result == out_path
        parsed = get_records_from_fasta(out_path)
        assert len(parsed) == 1

    def test_filters_none_records(self, tmp_path):
        records = [
            SeqRecord(Seq("MK"), id="P12345"),
            None,
            SeqRecord(Seq("AIL"), id="Q67890"),
        ]
        out_path = str(tmp_path / "out.fasta")
        write_fasta_from_records(records, out_path)
        parsed = get_records_from_fasta(out_path)
        assert len(parsed) == 2

    def test_empty_records_creates_empty_file(self, tmp_path):
        out_path = str(tmp_path / "out.fasta")
        write_fasta_from_records([], out_path)
        parsed = get_records_from_fasta(out_path)
        assert parsed == []


class TestReadYamlConfig:
    def test_reads_simple_config(self, tmp_path):
        config_path = tmp_path / "config.yaml"
        config_path.write_text("key: value\nnumber: 42\n")
        result = read_yaml_config(str(config_path))
        assert result == {"key": "value", "number": 42}

    def test_reads_nested_config(self, tmp_path):
        config_path = tmp_path / "config.yaml"
        config_path.write_text("global:\n  iterations: 5\nrun:\n  seed_uid: P12345\n")
        result = read_yaml_config(str(config_path))
        assert result["global"]["iterations"] == 5
        assert result["run"]["seed_uid"] == "P12345"


class TestCopyFile:
    def test_copies_file_successfully(self, tmp_path):
        src = tmp_path / "src.txt"
        src.write_text("hello")
        dst = str(tmp_path / "subdir" / "dst.txt")
        assert copy_file(str(src), dst)
        assert os.path.exists(dst)
        with open(dst) as f:
            assert f.read() == "hello"

    def test_returns_false_on_missing_source(self, tmp_path):
        assert not copy_file(str(tmp_path / "nonexistent.txt"), str(tmp_path / "dst.txt"))


class TestLogInterval:
    def test_small_total(self):
        assert _log_interval(3) == 1

    def test_medium_total(self):
        assert _log_interval(50) == 10

    def test_large_total(self):
        assert _log_interval(1000) == 10

    def test_zero_total(self):
        assert _log_interval(0) == 1


class TestBuildUniprotUrl:
    def test_uniprotkb_id(self):
        url = _build_uniprot_url("P12345")
        assert url == "https://rest.uniprot.org/uniprotkb/P12345.fasta"

    def test_uniparc_id(self):
        url = _build_uniprot_url("UPI0000123456")
        assert url == "https://rest.uniprot.org/uniparc/UPI0000123456.fasta"


class TestChunked:
    def test_even_split(self):
        assert list(_chunked([1, 2, 3, 4], 2)) == [[1, 2], [3, 4]]

    def test_uneven_split(self):
        assert list(_chunked([1, 2, 3, 4, 5], 2)) == [[1, 2], [3, 4], [5]]

    def test_chunk_size_larger_than_list(self):
        assert list(_chunked([1, 2], 10)) == [[1, 2]]

    def test_empty_list(self):
        assert list(_chunked([], 5)) == []

    def test_single_item(self):
        assert list(_chunked([42], 1)) == [[42]]


class TestBuildBatchQuery:
    def test_uniprotkb_ids(self):
        result = _build_batch_query(["P12345", "Q67890"], is_uniparc=False)
        assert result == "accession_id:P12345 OR accession_id:Q67890"

    def test_uniparc_ids(self):
        result = _build_batch_query(["UPI0001", "UPI0002"], is_uniparc=True)
        assert result == "upi:UPI0001 OR upi:UPI0002"

    def test_single_id(self):
        result = _build_batch_query(["P12345"], is_uniparc=False)
        assert result == "accession_id:P12345"

    def test_empty_list(self):
        result = _build_batch_query([], is_uniparc=False)
        assert result == ""


class TestCollectCachedRecords:
    def test_all_cached(self, tmp_path):
        cache_dir = str(tmp_path / "cache")
        os.makedirs(cache_dir)
        records = [SeqRecord(Seq("MK"), id="P12345")]
        write_fasta_from_records(records, os.path.join(cache_dir, "P12345.fasta"))

        found, uncached = _collect_cached_records(["P12345"], cache_dir)
        assert len(found) == 1
        assert found[0].id == "P12345"
        assert uncached == []

    def test_mixed_cached_and_uncached(self, tmp_path):
        cache_dir = str(tmp_path / "cache")
        os.makedirs(cache_dir)
        records = [SeqRecord(Seq("MK"), id="P12345")]
        write_fasta_from_records(records, os.path.join(cache_dir, "P12345.fasta"))

        found, uncached = _collect_cached_records(["P12345", "Q67890"], cache_dir)
        assert len(found) == 1
        assert uncached == ["Q67890"]

    def test_none_cached(self, tmp_path):
        cache_dir = str(tmp_path / "cache")
        os.makedirs(cache_dir)
        found, uncached = _collect_cached_records(["P12345", "Q67890"], cache_dir)
        assert found == []
        assert uncached == ["P12345", "Q67890"]

    def test_empty_accessions(self, tmp_path):
        cache_dir = str(tmp_path / "cache")
        os.makedirs(cache_dir)
        found, uncached = _collect_cached_records([], cache_dir)
        assert found == []
        assert uncached == []