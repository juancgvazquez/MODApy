from MODApy import vcfanalysis

import pandas as pd


def test_single_succesful(single_result, test_data_path):
    # Run the single function
    result = vcfanalysis.single(
        str(test_data_path / "test_pat1.vcf"), str(test_data_path / "test_panel.xlsx")
    )
    df = pd.read_excel(result)
    assert df.equals(single_result)


def test_duos_with_panel_succesful(duos_result, test_data_path):
    # Run the duos function
    result_path = vcfanalysis.duos(
        str(test_data_path / "test_pat1.vcf"),
        str(test_data_path / "test_pat2.vcf"),
        VennPlace="A",
        Panel=str(test_data_path / "test_panel.xlsx"),
        Filter=["ZIGOSITY HOM"],
    )
    df = pd.read_excel(result_path)
    assert df.equals(duos_result)


def test_trios_with_panel_succesful(trios_result, test_data_path):
    # Run the trios function
    result_path = vcfanalysis.trios(
        str(test_data_path / "test_pat1.vcf"),
        str(test_data_path / "test_pat2.vcf"),
        str(test_data_path / "test_pat3.vcf"),
        VennPlace="A:B:C",
        Panel=str(test_data_path / "test_panel.xlsx"),
        Filter=["VARTYPE SNP"],
    )
    df = pd.read_excel(result_path)
    assert df.equals(trios_result)
