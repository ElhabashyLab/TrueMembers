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