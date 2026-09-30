from pathlib import Path

import pytest

TEST_DATA_DIR = Path(__file__).resolve().parent.parent / "test_data"


@pytest.fixture
def annotation_csv_small() -> str:
    return str(TEST_DATA_DIR / "output_evcouplings" / "Q5SK37_507-569_b0.8_annotation.csv")


@pytest.fixture
def evcouplings_test_dir() -> str:
    return str(TEST_DATA_DIR / "output_evcouplings")


@pytest.fixture
def hmms_dir() -> str:
    return str(TEST_DATA_DIR / "test_hmm_domain_filters")


@pytest.fixture
def mmseqs_test_dir() -> str:
    return str(TEST_DATA_DIR / "output_mmseq2")


@pytest.fixture
def mmseqs_fasta_1(mmseqs_test_dir) -> str:
    return str(Path(mmseqs_test_dir) / "test_file_1.fasta")


@pytest.fixture
def mmseqs_fasta_2(mmseqs_test_dir) -> str:
    return str(Path(mmseqs_test_dir) / "test_file_2.fasta")


@pytest.fixture
def mmseqs_db1(mmseqs_test_dir) -> str:
    return str(Path(mmseqs_test_dir) / "db1.db")


@pytest.fixture
def mmseqs_db2(mmseqs_test_dir) -> str:
    return str(Path(mmseqs_test_dir) / "db2.db")


@pytest.fixture
def mmseqs_db3(mmseqs_test_dir) -> str:
    return str(Path(mmseqs_test_dir) / "db3.db")


@pytest.fixture
def mmseqs_result_tsv(mmseqs_test_dir) -> str:
    return str(Path(mmseqs_test_dir) / "resultstsv")