import logging
from collections import OrderedDict

import cyvcf2
import pandas as pd

log = logging.getLogger(__name__)


class ParsedVCF(pd.DataFrame):
    """ParsedVCF class is an extension of a Pandas DataFrame class"""

    @staticmethod
    def _divide(x, y):
        """
        Method to divide x on y, needed for dividing freqs.
        Parameters
        ----------
        x
            The dividend
        y
            The divisor
        Returns result or x.
        """
        try:
            return float(x) / y
        except:
            return x

    @property
    def _constructor(self):
        """
        Method needed to keep the resulting objects as ParsedVCF type.
        """
        return ParsedVCF

    @classmethod
    def from_vcf(cls, vcf):
        """
        Method that creates a ParsedVCF (a DataFrame) from a vcf file
        Parameters
        ----------
        vcf
            Path to the vcf to parse.
        """
        pVCF = cyvcf2.Reader(vcf)
        variantsDict = OrderedDict()
        for variant in pVCF:
            variantsDict[variant.CHROM + '+' + str(variant.POS) + '+' + variant.REF + '+' + ','.join(variant.ALT)] = {
                'ID': variant.ID, 'QUAL': variant.QUAL, 'FILTER': variant.FILTER}
            variantsDict[
                variant.CHROM + '+' + str(variant.POS) + '+' + variant.REF + '+' + ','.join(variant.ALT)].update(
                {k: v for (k, v) in variant.INFO})

        df1 = pd.DataFrame.from_dict(variantsDict, orient='index')
        df1.index = df1.index.str.split('+', expand=True)
        df1.index.names = ['CHROM', 'POS', 'REF', 'ALT']
        df1.reset_index(inplace=True)
        splitdf = df1.loc[df1['ALT'].str.contains(',') == True].copy()
        ALT = splitdf['ALT'].astype(str).str.split(',', n=1, expand=True).stack().rename('ALT')
        ALT.index = ALT.index.droplevel(-1)
        ALT = ALT.to_frame()
        splitdf = splitdf.join(ALT, lsuffix='_x', rsuffix='_y')
        splitdf['ALT'] = splitdf['ALT_y'].combine_first(splitdf['ALT_x'])
        splitdf.drop(columns=['ALT_y', 'ALT_x'], inplace=True)
        splitdf.reset_index(inplace=True)
        splitdf.drop(columns='index', inplace=True)
        odd = splitdf.iloc[::2].copy()
        even = splitdf.iloc[1::2].copy()
        splitlist = ['ID', 'AC', 'SAMPLES_AF', 'MLEAC', 'MLEAF', 'VARTYPE', 'dbSNPBuildID']
        if 'ANN' in splitlist:
            splitlist.remove('ANN')
        for col in splitlist:
            odd[col] = odd[col].astype(str).str.split(',', n=1).str[0]
            even[col] = even[col].apply(lambda x:
                                        x if len(str(x).split(',')) <= 1 else str(x).split(',', maxsplit=1)[1])
        splitdf = pd.concat([odd, even]).sort_index().replace(to_replace=['\(', '\)'], value='', regex=True)
        splitdf = splitdf[['CHROM', 'POS', 'REF', 'ALT'] + splitlist]
        df1 = df1.merge(splitdf, on=['CHROM', 'POS', 'REF'], how='left')
        splitlist.append('ALT')
        xlist = [x + '_x' for x in splitlist]
        ylist = [y + '_y' for y in splitlist]

        for col in splitlist:
            df1[col] = df1[col + '_y'].combine_first(df1[col + '_x'])
        df1.drop(columns=xlist + ylist, inplace=True)
        df1['POS'] = df1['POS'].astype(int)
        df1.sort_values(by=['CHROM', 'POS'], inplace=True)
        if 'ANN' in df1.columns:
            anndf = df1['ANN']
            annhead = pVCF.get_header_type('ANN')['Description'].strip('"Functional annotations: \'"')
            annheaderlist = [x.strip() for x in annhead.split('|')]
            anndf = anndf.str.split(',', expand=True).stack()
            anndf = anndf.str.split('|', expand=True)
            anndf.columns = annheaderlist
        vcfdf = df1.copy()
        vcfdf.drop(columns=['ANN'], inplace=True)
        anndf.index = anndf.index.droplevel(1)
        vcfdf = vcfdf.join(anndf, how='inner')
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
            vcfdf = vcfdf.where(pd.notnull(vcfdf), None)
            vcfdf.sort_values(by=['CHROM', 'POS']).reset_index(inplace=True)
        IMPACT_SEVERITY = {
            'exon_loss_variant': 1,
            'frameshift_variant': 2,
            'stop_gained': 3,
            'stop_lost': 4,
            'start_lost': 5,
            'splice_acceptor_variant': 6,
            'splice_donor_variant': 7,
            'disruptive_inframe_deletion': 8,
            'inframe_insertion': 9,
            'disruptive_inframe_insertion': 10,
            'inframe_deletion': 11,
            'missense_variant': 12,
            'splice_region_variant': 13,
            'stop_retained_variant': 14,
            'initiator_codon_variant': 15,
            'synonymous_variant': 16,
            'start_retained': 17,
            'coding_sequence_variant': 18,
            '5_prime_UTR_variant': 19,
            '3_prime_UTR_variant': 20,
            '5_prime_UTR_premature_start_codon_gain_variant': 21,
            'intron_variant': 22,
            'non_coding_exon_variant': 23,
            'upstream_gene_variant': 24,
            'downstream_gene_variant': 25,
            'TF_binding_site_variant': 26,
            'regulatory_region_variant': 27,
            'intergenic_region': 28,
            'transcript': 29
        }
        vcfdf2 = vcfdf.copy()
        vcfdf2['sorter'] = vcfdf2['EFFECT'].str.split('&').str[0].replace(IMPACT_SEVERITY)
        vcfdf2.loc[vcfdf2['HGVS.C'].str.contains('null'), 'HGVS.C'] = None
        vcfdf2['sorter2'] = [x[0] == x[1] for x in zip(vcfdf2['ALT'], vcfdf2['ALLELE'])]
        vcfdf2['chrsort'] = vcfdf2['CHROM'].str.strip('chr').replace({'X': 30, 'Y': 40}).astype(int)
        vcfdf2 = vcfdf2.sort_values(by=['CHROM', 'POS', 'sorter2', 'sorter'],
                                    ascending=[True, True, False, True]).drop_duplicates(['CHROM', 'POS', 'REF', 'ALT'])
        result = vcfdf2.pipe(ParsedVCF)
        try:
            result.name = pVCF.samples[0]
        except:
            result.name = vcf.split('/')[-1]

        return result

    def to_macrogen_xls(self, outpath):
        macrogen_cols = ['CHROM', 'POS', 'REF', 'ALT', 'DP', 'AD', 'QUAL', 'MQ', 'Zygosity', 'FILTER', 'Effect',
                         'Putative_Impact', 'Gene_Name', 'Feature_Type', 'Feature_ID', 'Transcript_BioType',
                         'Rank/Total', 'HGVS.c', 'HGVS.p', 'REF_AA', 'ALT_AA', 'cDNA_pos', 'cDNA_length', 'CDS_pos',
                         'CDS_length', 'AA_pos', 'AA_length', 'Distance', 'dbSNP142_ID', '1000Gp3_AF', '1000Gp3_AFR_AF',
                         '1000Gp3_AMR_AF', '1000Gp3_EAS_AF', '1000Gp3_EUR_AF', '1000Gp3_SAS_AF', 'ESP6500_MAF_EA',
                         'ESP6500_MAF_AA', 'ESP6500_MAF_ALL', 'SIFT_score', 'SIFT_pred', 'Polyphen2_HDIV_score',
                         'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred', 'CLINVAR_CLNSIG',
                         'CLINVAR_CLNDSDB', 'CLINVAR_CLNDSDBID', 'CLINVAR_CLNDBN', 'CLINVAR_CLNREVSTAT',
                         'CLINVAR_CLNACC']
        df1 = self.sort_values(['chrsort', 'POS'])[
            [x.upper() for x in macrogen_cols if x.upper() in self.columns]].copy()
        df1.to_excel(outpath)

    # TODO:USE THIS PANEL!
    def panel(self, panel):
        pldf = pd.ExcelFile(panel).parse('GeneList')
        geneSymbolList = list(pldf.GeneSymbol.unique())
        panel_df = pd.DataFrame()
        for gene in geneSymbolList:
            panel_df = panel_df.append(self.loc[self['GENE_NAME'] == gene])
        panel_df.name = self.name
        panel_df = panel_df.pipe(ParsedVCF)
        return panel_df

    def duos(self, vcf2):
        """
        Method to compare two vcfs, using CHROM, POS, REF and ALT columns as index.
        Parameters
        ----------
        vcf2
            VCF file to compare to.

        Returns a Dataframe containing a new column 'DUOS', that indicates in which file is the variant.
        """
        vcf2df = ParsedVCF.from_vcf(vcf2)
        self.set_index(['CHROM', 'POS', 'REF', 'ALT'], inplace=True)
        vcf2df.set_index(['CHROM', 'POS', 'REF', 'ALT'], inplace=True)
        AyB = self.loc[self.index.isin(vcf2df.index)].copy()
        A = self.loc[~self.index.isin(vcf2df.index)].copy()
        B = vcf2df.loc[~vcf2df.index.isin(self.index)].copy()
        AyB['DUOS'] = ':'.join([self.name, vcf2df.name])
        A['DUOS'] = self.name
        B['DUOS'] = vcf2df.name
        duos_df = pd.concat([A, B, AyB], sort=False)
        droplist = ['AC', 'SAMPLES_AF', 'AN', 'DP', 'FS', 'MLEAC', 'MLEAF', 'MQ', 'QD', 'SOR', 'DBSNPBUILDID']
        duos_df.drop(columns=droplist, inplace=True, errors='ignore')
        duos_df.name = ':'.join([self.name, vcf2df.name])
        return duos_df

    def trios(self, vcf2, vcf3):
        """
        Method to compare three vcfs, using CHROM, POS, REF and ALT columns as index.
        Parameters
        ----------
        vcf2
            VCF file to compare to.

        vcf3
            VCF file to compare to.

        Returns a Dataframe containing a new column 'DUOS', that indicates in which file is the variant.
        """
        vcf2df = ParsedVCF.from_vcf(vcf2)
        vcf3df = ParsedVCF.from_vcf(vcf3)
        self.set_index(['CHROM', 'POS', 'REF', 'ALT'], inplace=True)
        vcf2df.set_index(['CHROM', 'POS', 'REF', 'ALT'], inplace=True)
        vcf3df.set_index(['CHROM', 'POS', 'REF', 'ALT'], inplace=True)
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
        trios_df = pd.concat(dfs, sort=False)
        droplist = ['AC', 'SAMPLES_AF', 'AN', 'DP', 'FS', 'MLEAC', 'MLEAF', 'MQ', 'QD', 'SOR', 'DBSNPBUILDID']
        trios_df.drop(columns=droplist, inplace=True, errors='ignore')
        trios_df.name = ':'.join([self.name, vcf2df.name, vcf3df.name])
        return trios_df
