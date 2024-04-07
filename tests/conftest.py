import pickle
from pathlib import Path

import pytest


TEST_DATA_PATH = Path("tests/test_data/")


@pytest.fixture(scope="session")
def test_data_path():
    return TEST_DATA_PATH


@pytest.fixture(scope="session")
def parsed_vcf():
    with open(TEST_DATA_PATH / "parsed_vcf_pat1.pkl", "rb") as f:
        return pickle.load(f)


@pytest.fixture(scope="session")
def single_result():
    with open(TEST_DATA_PATH / "single_pat1.pkl", "rb") as f:
        return pickle.load(f)


@pytest.fixture(scope="session")
def duos_result():
    with open(TEST_DATA_PATH / "duos_pat1_pat2.pkl", "rb") as f:
        return pickle.load(f)


@pytest.fixture(scope="session")
def trios_result():
    with open(TEST_DATA_PATH / "trios_pat1_pat2_pat3.pkl", "rb") as f:
        return pickle.load(f)
