import pandas as pd

from MODApy.vcfmgr import ParsedVCF


def panelrun(panel, vcffile):
    pldf = pd.ExcelFile(panel).parse('GeneList')
    gsymbollist = list(pldf.GeneSymbol.unique())
    if type(vcffile) == ParsedVCF:
        vcfdf = vcffile
    elif type(vcffile) == str:
        vcfdf = ParsedVCF.from_vcf(vcffile)

    df_final = check_panel(gsymbollist, vcfdf)
    df_final.name = vcfdf.name

    return df_final


def check_panel(genelist, vcf: pd.DataFrame):
    result = pd.DataFrame()
    for gene in genelist:
        result = result.append(vcf.loc[vcf['GENE_NAME'] == gene])
    return result
