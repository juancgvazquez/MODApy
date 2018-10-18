import pandas as pd
from MoDAPy import vcfmgr


def panelrun(panel, vcffile):
	pldf = pd.ExcelFile(panel).parse('GeneList')
	gsymbollist = list(pldf.GeneSymbol.unique())
	if type(vcffile) == pd.DataFrame:
		vcfdf = vcffile
	elif type(vcffile) == str:
		vcfdf = vcfmgr.vcf_to_df(vcffile)

	df_final = check_panel(gsymbollist,vcfdf)
	df_final.name = vcfdf.name

	return df_final


def check_panel(genelist, vcf: pd.DataFrame):
	result = pd.DataFrame()
	for gene in genelist:
		result = result.append(vcf.loc[vcf['Gene_Name'].str.contains(gene)])
	return result
