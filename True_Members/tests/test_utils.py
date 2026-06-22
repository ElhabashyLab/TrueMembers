"""Unit tests for scripts/utils.py — ID parsing, region parsing, FASTA operations."""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from True_Members.scripts.utils import (
    extract_uid_from_recordID,
    extract_region_from_recordID,
    filter_fasta_file,
    write_fasta_from_records,
    get_records_from_fasta,
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