"""Unit tests for scripts/uniprot_local.py — UniprotLocalDB indexing and lookup."""

import os
import tempfile

from True_Members.scripts.uniprot_local import UniprotLocalDB

MOCK_FASTA = """>sp|P00533|EGFR_HUMAN Epidermal growth factor receptor
MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEV
VLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLENLQIIRGNMYYENSYALA
>tr|A0A6P2YNG8|A0A6P2YNG8_HUMAN Uncharacterized protein
MSTGDSFETRFEIPIWVDANDNFPGATVVEVTLGGPSGTPVTVEVGDTLKLTVTGLKEDR
GLSLDRLLEVAGQVGERVTVRLDGDRVSSRLVRGGQDTIEVTAANLGQSSSRSVTV
>sp|P04637|P53_HUMAN Cellular tumor antigen p53
MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK
>tr|Q9Y6K9|Q9Y6K9_HUMAN IKBKG
MNRHLWKSQLCEMVKPSLKRHPTLF
>sp|P68871|HBB_HUMAN Hemoglobin subunit beta
MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLS
"""


class TestUniprotLocalDB:
    def setup_method(self):
        self.tmpdir = tempfile.mkdtemp(prefix="uniprot_test_")
        self.fasta_path = os.path.join(self.tmpdir, "test.fasta")
        with open(self.fasta_path, "w") as f:
            f.write(MOCK_FASTA)

    def teardown_method(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_build_index_creates_idx_file(self):
        db = UniprotLocalDB(self.fasta_path)
        db.ensure_index()
        assert os.path.exists(self.fasta_path + ".idx")

    def test_index_get_sequence_returns_record(self):
        db = UniprotLocalDB(self.fasta_path)
        record = db.get_sequence("P00533")
        assert record is not None
        assert record.id == "sp|P00533|EGFR_HUMAN"
        assert record.seq.startswith("MRPSGTAG")

    def test_get_sequence_tr_entry(self):
        db = UniprotLocalDB(self.fasta_path)
        record = db.get_sequence("A0A6P2YNG8")
        assert record is not None
        assert record.id == "tr|A0A6P2YNG8|A0A6P2YNG8_HUMAN"

    def test_get_sequence_missing_accession(self):
        db = UniprotLocalDB(self.fasta_path)
        record = db.get_sequence("P99999")
        assert record is None

    def test_get_sequence_short_protein(self):
        db = UniprotLocalDB(self.fasta_path)
        record = db.get_sequence("Q9Y6K9")
        assert record is not None
        assert len(record.seq) == 25

    def test_get_sequence_multi_line(self):
        db = UniprotLocalDB(self.fasta_path)
        record = db.get_sequence("P04637")
        assert record is not None
        assert len(record.seq) > 40

    def test_index_reuse(self):
        db1 = UniprotLocalDB(self.fasta_path)
        db1.ensure_index()
        idx_mtime = os.path.getmtime(self.fasta_path + ".idx")

        db2 = UniprotLocalDB(self.fasta_path)
        db2.ensure_index()
        assert os.path.getmtime(self.fasta_path + ".idx") == idx_mtime

    def test_all_accessions_resolve(self):
        db = UniprotLocalDB(self.fasta_path)
        for acc in ["P00533", "A0A6P2YNG8", "P04637", "Q9Y6K9", "P68871"]:
            assert db.get_sequence(acc) is not None, f"{acc} not found"

    def test_binary_search_correctness(self):
        db = UniprotLocalDB(self.fasta_path)
        db.ensure_index()
        for acc in ["P00533", "A0A6P2YNG8", "P04637", "Q9Y6K9", "P68871"]:
            offset = db._binary_search_offset(acc)
            assert offset is not None, f"{acc} not found by binary search"

    def test_binary_search_missing(self):
        db = UniprotLocalDB(self.fasta_path)
        db.ensure_index()
        assert db._binary_search_offset("P99999") is None
        assert db._binary_search_offset("ZZZZZZ") is None
