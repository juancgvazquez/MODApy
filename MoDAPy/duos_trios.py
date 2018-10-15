import pandas as pd
from MoDAPy import vcfmgr


def duos(vcf1, vcf2):
	vcf1df = vcfmgr.vcf_to_df(vcf1) # type:pd.DataFrame
	vcf2df = vcfmgr.vcf_to_df(vcf2) # type:pd.DataFrame
	AyB = vcf1df.loc[vcf1df.index.isin(vcf2df.index)].copy()
	A = vcf1df.loc[~vcf1df.index.isin(vcf2df.index)].copy()
	B = vcf2df.loc[~vcf2df.index.isin(vcf1df.index)].copy()
	AyB['Duos'] = ':'.join([vcf1df.name, vcf2df.name])
	A['Duos'] = vcf1df.name
	B['Duos'] = vcf2df.name
	df_final = pd.concat([A,B,AyB], sort=False)
	droplist = ['AC', 'SAMPLES_AF', 'AN', 'DP', 'FS', 'MLEAC', 'MLEAF', 'MQ', 'QD', 'SOR', 'dbSNPBuildID']
	df_final = df_final.drop(columns=droplist)
	df_final.name = ':'.join([vcf1df.name, vcf2df.name])
	return df_final


def trios(vcf1, vcf2, vcf3):
	vcf1df = vcfmgr.vcf_to_df(vcf1) # type:pd.DataFrame
	vcf2df = vcfmgr.vcf_to_df(vcf2) # type:pd.DataFrame
	vcf3df = vcfmgr.vcf_to_df(vcf3) # type:pd.DataFrame
	A = vcf1df.loc[(~vcf1df.index.isin(vcf2df.index)) & (~vcf1df.index.isin(vcf3df.index))].copy()
	A['Trios'] = vcf1df.name
	B = vcf2df.loc[(~vcf2df.index.isin(vcf1df.index)) & (~vcf2df.index.isin(vcf3df.index))].copy()
	B['Trios'] = vcf2df.name
	C = vcf3df.loc[(~vcf3df.index.isin(vcf1df.index)) & (~vcf3df.index.isin(vcf2df.index))].copy()
	C['Trios'] = vcf3df.name
	ABC = vcf1df.loc[(vcf1df.index.isin(vcf2df.index)) & (vcf1df.index.isin(vcf3df.index))].copy()
	ABC['Trios'] = ':'.join([vcf1df.name,vcf2df.name,vcf3df.name])
	AB = vcf1df.loc[(vcf1df.index.isin(vcf2df.index)) & (~vcf1df.index.isin(vcf3df.index))].copy()
	AB['Trios'] = ':'.join([vcf1df.name,vcf2df.name])
	AC = vcf1df.loc[(vcf1df.index.isin(vcf3df.index)) & (~vcf1df.index.isin(vcf2df.index))].copy()
	AC['Trios'] = ':'.join([vcf1df.name,vcf3df.name])
	BC = vcf2df.loc[(vcf2df.index.isin(vcf3df.index)) & (~vcf2df.index.isin(vcf1df.index))].copy()
	BC['Trios'] = ':'.join([vcf2df.name,vcf3df.name])
	dfs = [A,B,C,AB,AC,BC,ABC]
	df_final = pd.concat(dfs, sort=False)
	droplist = ['AC','SAMPLES_AF','AN','DP','FS','MLEAC','MLEAF','MQ','QD','SOR','dbSNPBuildID']
	df_final = df_final.drop(columns=droplist)
	df_final.name = ':'.join([vcf1df.name,vcf2df.name,vcf3df.name])
	return df_final