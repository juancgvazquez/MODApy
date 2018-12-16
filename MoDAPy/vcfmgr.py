from collections import OrderedDict

import cyvcf2
import pandas as pd


class ParsedVCF(pd.DataFrame):

    @staticmethod
    def _divide(x, y):
        try:
            return float(x) / y
        except:
            return x

    @property
    def _constructor(self):
        return ParsedVCF

    @classmethod
    def from_vcf(cls, vcf):
        try:
            pVCF = cyvcf2.Reader(vcf)
        except:
            print('error loading vcf')
            exit(1)

        variantsDict = OrderedDict()
        for variant in pVCF:
            variantsDict[variant.CHROM + '_' + str(variant.POS) + '_' + variant.REF + '_' + ','.join(variant.ALT)] = {
                'ID': variant.ID, 'QUAL': variant.QUAL, 'FILTER': variant.FILTER}
            variantsDict[
                variant.CHROM + '_' + str(variant.POS) + '_' + variant.REF + '_' + ','.join(variant.ALT)].update(
                {k: v for (k, v) in variant.INFO})

        df1 = pd.DataFrame.from_dict(variantsDict, orient='index')
        if 'ANN' in df1.columns:
            anndf = df1['ANN']
            annhead = pVCF.get_header_type('ANN')['Description'].strip('"Functional annotations: \'"')
            annheaderlist = [x.strip() for x in annhead.split('|')]
            anndf = anndf.str.split(',', expand=True).stack()
            anndf = anndf.str.split('|', expand=True)
            anndf.columns = annheaderlist
            df1.drop(columns=['ANN'], inplace=True)
            anndf.index = anndf.index.droplevel(1)
            vcfdf = df1.join(anndf, how='inner')
        else:
            vcfdf = df1

        vcfdf.index = vcfdf.index.str.split('_', expand=True)
        vcfdf.index.names = ['CHROM', 'POS', 'REF', 'ALT']
        vcfdf.columns = vcfdf.columns.str.upper()
        if 'HOM' in vcfdf.columns:
            vcfdf['HOM'] = vcfdf['HOM'].map({True: 'HOM'})
            vcfdf.HOM.fillna('HET', inplace=True)
            vcfdf.rename(columns={'HOM': 'ZIGOSITY'}, inplace=True)
        if 'ESP6500_MAF' in vcfdf.columns:
            vcfdf[['ESP6500_MAF_EA', 'ESP6500_MAF_AA', 'ESP6500_MAF_ALL']] = vcfdf['ESP6500_MAF'].str.split(',',
                                                                                                            expand=True)
            vcfdf['ESP6500_MAF_EA'] = vcfdf['ESP6500_MAF_EA'].apply(cls._divide, args=(100,))
            vcfdf['ESP6500_MAF_AA'] = vcfdf['ESP6500_MAF_AA'].apply(cls._divide, args=(100,))
            vcfdf['ESP6500_MAF_ALL'] = vcfdf['ESP6500_MAF_ALL'].apply(cls._divide, args=(100,))
            vcfdf.drop(columns=['ESP6500_MAF'], inplace=True)
        if 'ESP6500_PH' in vcfdf.columns:
            vcfdf[['PolyPhen_Pred', 'PolyPhen_Score']] = vcfdf['ESP6500_PH'].str.split(':', 1, expand=True)
            vcfdf['PolyPhen_Pred'] = vcfdf['PolyPhen_Pred'].str.strip('.').str.strip('.,')
            vcfdf['PolyPhen_Score'] = vcfdf['PolyPhen_Score'].str.split(',').str[0]
            vcfdf.drop(columns=['ESP6500_PH'], inplace=True)
        vcfdf.rename(columns={'ANNOTATION': 'EFFECT', 'ANNOTATION_IMPACT': 'IMPACT', 'ID': 'RSID'},
                     inplace=True)
        vcfdf.fillna(inplace=True)
        result = vcfdf.pipe(ParsedVCF)
        try:
            result.name = pVCF.samples[0]
        except:
            result.name = vcf.split('/')[-1]

        return result

    def duos(self, vcf2):
        vcf2df = ParsedVCF.from_vcf(vcf2)
        AyB = self.loc[self.index.isin(vcf2df.index)].copy()
        A = self.loc[~self.index.isin(vcf2df.index)].copy()
        B = vcf2df.loc[~vcf2df.index.isin(self.index)].copy()
        AyB['DUOS'] = ':'.join([self.name, vcf2df.name])
        A['DUOS'] = self.name
        B['DUOS'] = vcf2df.name
        df_final = pd.concat([A, B, AyB], sort=False)
        droplist = ['AC', 'SAMPLES_AF', 'AN', 'DP', 'FS', 'MLEAC', 'MLEAF', 'MQ', 'QD', 'SOR', 'DBSNPBUILDID']
        df_final = df_final.drop(columns=droplist)
        df_final.name = ':'.join([self.name, vcf2df.name])
        return df_final

    def trios(self, vcf2, vcf3):

        vcf2df = ParsedVCF.from_vcf(vcf2)
        print(vcf2df.name)
        vcf3df = ParsedVCF.from_vcf(vcf3)
        print(vcf3df.name)
        A = self.loc[(~self.index.isin(vcf2df.index)) & (~self.index.isin(vcf3df.index))].copy()
        A['TRIOS'] = self.name
        B = vcf2df.loc[(~vcf2df.index.isin(self.index)) & (~vcf2df.index.isin(vcf3df.index))].copy()
        B['TRIOS'] = vcf2df.name
        C = vcf3df.loc[(~vcf3df.index.isin(self.index)) & (~vcf3df.index.isin(vcf2df.index))].copy()
        C['TRIOS'] = vcf3df.name
        ABC = self.loc[(self.index.isin(vcf2df.index)) & (self.index.isin(vcf3df.index))].copy()
        ABC['TRIOS'] = ':'.join([self.name, vcf2df.name, vcf3df.name])
        AB = self.loc[(self.index.isin(vcf2df.index)) & (~self.index.isin(vcf3df.index))].copy()
        AB['TRIOS'] = ':'.join([self.name, vcf2df.name])
        AC = self.loc[(self.index.isin(vcf3df.index)) & (~self.index.isin(vcf2df.index))].copy()
        AC['TRIOS'] = ':'.join([self.name, vcf3df.name])
        BC = vcf2df.loc[(vcf2df.index.isin(vcf3df.index)) & (~vcf2df.index.isin(self.index))].copy()
        BC['TRIOS'] = ':'.join([vcf2df.name, vcf3df.name])
        dfs = [A, B, AB, C, AC, BC, ABC]
        df_final = pd.concat(dfs, sort=False)
        droplist = ['AC', 'SAMPLES_AF', 'AN', 'DP', 'FS', 'MLEAC', 'MLEAF', 'MQ', 'QD', 'SOR', 'DBSNPBUILDID']
        df_final = df_final.drop(columns=droplist)
        df_final.name = ':'.join([self.name, vcf2df.name, vcf3df.name])
        return df_final
