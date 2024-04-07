import pickle

from MODApy import vcfanalysis, vcfmgr

import pandas as pd


def generate_parsed_vcf_data():
    """Generate a pickled ParsedVCF object from a VCF file."""
    df = vcfmgr.ParsedVCF.from_vcf("../test_data/test_pat1.vcf")
    with open("../test_data/test_pat1.pkl", "wb") as f:
        pickle.dump(df, f)


def generate_single_data():
    """Generate a pickled single result."""
    result = vcfanalysis.single(
        "../test_data/test_pat1.vcf", "../test_data/test_panel.xlsx"
    )
    df = pd.read_excel(result)
    with open("../test_data/single_pat1.pkl", "wb") as f:
        pickle.dump(df, f)


def generate_duos_data():
    """Generate a pickled duos result."""
    result_path = vcfanalysis.duos(
        "../test_data/test_pat1.vcf",
        "../test_data/test_pat2.vcf",
        VennPlace="A",
        Panel="../test_data/test_panel.xlsx",
        Filter=["ZIGOSITY HOM"],
    )
    df = pd.read_excel(result_path)
    with open("../test_data/duos_pat1_pat2.pkl", "wb") as f:
        pickle.dump(df, f)


def generate_trios_data():
    """Generate a pickled trios result."""
    result_path = vcfanalysis.trios(
        "../test_data/test_pat1.vcf",
        "../test_data/test_pat2.vcf",
        "../test_data/test_pat3.vcf",
        VennPlace="A:B:C",
        Panel="../test_data/test_panel.xlsx",
        Filter=["VARTYPE SNP"],
    )
    df = pd.read_excel(result_path)
    with open("../test_data/trios_pat1_pat2_pat3.pkl", "wb") as f:
        pickle.dump(df, f)


if __name__ == "__main__":
    generate_parsed_vcf_data()
    generate_single_data()
    generate_duos_data()
    generate_trios_data()
