from pathlib import Path

from MODApy.vcfmgr import ParsedVCF

TEST_DATA_PATH = Path("tests/test_data")


def test_from_vcf(parsed_vcf):
    # Create a temporary VCF file for testing
    df = ParsedVCF.from_vcf(str(TEST_DATA_PATH / "TEST.vcf"))
    print(df.head())
    print(parsed_vcf.head())
    # Check that the dataframe has the correct shape and columns
    assert isinstance(df, ParsedVCF)
    assert set(df.columns) == set(parsed_vcf.columns)
    assert df.shape == parsed_vcf.shape
    assert df.equals(parsed_vcf)
