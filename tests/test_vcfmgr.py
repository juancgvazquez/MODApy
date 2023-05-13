from pathlib import Path
import tempfile
import pandas as pd
from MODApy.vcfmgr import ParsedVCF

TEST_DATA_PATH = Path("tests/test_data")

def test_from_vcf():
    # Create a temporary VCF file for testing
    vcf = tempfile.NamedTemporaryFile(suffix=".vcf")
    # with open(vcf.name, "w") as f:
    #     f.write(
    #         "##fileformat=VCFv4.2\n"
    #         "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
    #         "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">\n"
    #         "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    #         "1\t12345\t.\tA\tC\t50\tPASS\tAC=1;AF=0.5;AN=2\n"
    #         "1\t23456\trs123\tT\tG\t100\tPASS\tAC=1,2;AF=0.5,0.75;AN=4\n"
    #     )

    # Test the from_vcf method
    df = ParsedVCF.from_vcf(str(TEST_DATA_PATH / "TEST.vcf"))

    # Check that the dataframe has the correct shape and columns
    assert isinstance(df, pd.DataFrame)
    assert df.shape == (2, 11)
    assert set(df.columns) == {
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "ID",
        "QUAL",
        "FILTER",
        "AC",
        "AF",
        "AN",
        "VARTYPE",
    }

    # Check that the values in the dataframe are correct
    assert df.loc[0, "CHROM"] == "1"
    assert df.loc[0, "POS"] == 12345
    assert df.loc[0, "REF"] == "A"
    assert df.loc[0, "ALT"] == "C"
    assert df.loc[0, "ID"] == "."
    assert df.loc[0, "QUAL"] == 50.0
    assert df.loc[0, "FILTER"] == "PASS"
    assert df.loc[0, "AC"] == 1
    assert df.loc[0, "AF"] == 0.5
    assert df.loc[0, "AN"] == 2
    assert df.loc[0, "VARTYPE"] == "SNV"

    assert df.loc[1, "CHROM"] == "1"
    assert df.loc[1, "POS"] == 23456
    assert df.loc[1, "REF"] == "T"
    assert df.loc[1, "ALT"] == "G"
    assert df.loc[1, "ID"] == "rs123"
    assert df.loc[1, "QUAL"] == 100.0
    assert df.loc[1, "FILTER"] == "PASS"
    assert df.loc[1, "AC"] == 1
    assert df.loc[1, "AF"] == 0.5
    assert df.loc[1, "AN"] == 4
    assert df.loc[1, "VARTYPE"] == "SNV"
