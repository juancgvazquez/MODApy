from pathlib import Path
import pickle
import pytest


TEST_DATA_PATH = Path("tests/test_data/")


@pytest.fixture(scope="session")
def parsed_vcf():
    with open(TEST_DATA_PATH / "parsed_vcf.pkl", "rb") as f:
        return pickle.load(f)