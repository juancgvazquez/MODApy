import logging
import multiprocessing as mp
from collections import OrderedDict

import cyvcf2
import numpy as np
import pandas as pd

from MODApy.cfg import cfg

logger = logging.getLogger(__name__)


class ParsedVCF(pd.DataFrame):
    _metadata = ['name']

    @property
    def _constructor(self):
        return ParsedVCF

    @classmethod
    def from_vcf(cls, vcf):
        """
        Method that creates a ParsedVCF1 (a DataFrame) from a vcf file
        Parameters
        ----------
        vcf
            Path to the vcf to parse.
        """

        def aminoChange(value: str):
            try:
                value = value.replace('p.', '')
                if value[:3] != value[-3:]:
                    return 'CHANGE'
                else:
                    return '.'
            except:
                return '.'

        def divide(x, y):
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

        logger.info('Parsing VCF File. %s' % vcf)
        pVCF = cyvcf2.Reader(vcf)
        try:
            name = pVCF.samples[0]
        except:
            name = vcf.split('/')[-1]
        variantsDict = OrderedDict()
        for variant in pVCF:
            variantsDict[variant.CHROM + '+' + str(variant.POS) + '+' + variant.REF + '+' + ','.join(variant.ALT)] = {
                'ID': variant.ID, 'QUAL': variant.QUAL, 'FILTER': variant.FILTER}
            variantsDict[
                variant.CHROM + '+' + str(variant.POS) + '+' + variant.REF + '+' + ','.join(variant.ALT)].update(
                {k: v for (k, v) in variant.INFO})

        df1 = pd.DataFrame.from_dict(variantsDict, orient='index')
        del variantsDict
        df1.index = df1.index.str.split('+', expand=True)
        df1.index.names = ['CHROM', 'POS', 'REF', 'ALT']
        df1.reset_index(inplace=True)
        splitdf = df1.loc[df1['ALT'].str.contains(',') == True].copy()
        ALT = splitdf['ALT'].astype(str).str.split(',', n=1, expand=True).stack().rename('ALT')
        ALT.index = ALT.index.droplevel(-1)
        ALT = ALT.to_frame()
        splitdf = splitdf.join(ALT, lsuffix='_x', rsuffix='_y')
        del ALT
        splitdf['ALT'] = splitdf['ALT_y'].combine_first(splitdf['ALT_x'])
        splitdf.drop(columns=['ALT_y', 'ALT_x'], inplace=True)
        splitdf.reset_index(inplace=True)
        splitdf.drop(columns='index', inplace=True)
        odd = splitdf.iloc[::2].copy()
        even = splitdf.iloc[1::2].copy()
        splitlist = ['ID', 'AC', 'SAMPLES_AF', 'MLEAC', 'MLEAF', 'VARTYPE', 'dbSNPBuildID']
        splitlist += [x for x in df1.columns if x.startswith(('1000', 'CLINVAR'))]
        for col in splitlist:
            odd[col] = odd[col].astype(str).str.split(',', n=1).str[0]
            even[col] = even[col].apply(lambda x:
                                        x if len(str(x).split(',')) <= 1 else str(x).split(',', maxsplit=1)[1])
        splitdf = pd.concat([odd, even]).sort_index().replace(to_replace=['\(', '\)'], value='', regex=True)
        del odd, even
        splitdf = splitdf[['CHROM', 'POS', 'REF', 'ALT'] + splitlist]
        df1 = df1.merge(splitdf, on=['CHROM', 'POS', 'REF'], how='left')
        splitlist.append('ALT')
        xlist = [x + '_x' for x in splitlist]
        ylist = [y + '_y' for y in splitlist]
        del splitdf  # ya no uso más splitdf así que lo borro
        for col in splitlist:
            df1[col] = df1[col + '_y'].combine_first(df1[col + '_x'])
        del splitlist  # ya no uso más splitlist
        df1.drop(columns=xlist + ylist, inplace=True)
        del xlist, ylist  # ya no uso más esto.
        df1['POS'] = df1['POS'].astype(int)
        if 'ANN' in df1.columns:
            anndf = df1['ANN']
            annhead = pVCF.get_header_type('ANN')['Description'].strip('"Functional annotations: \'"')
            annheaderlist = [x.strip() for x in annhead.split('|')]
            anndf = anndf.str.split(',', expand=True).stack()
            anndf = anndf.str.split('|', expand=True)
            anndf.columns = annheaderlist
            df1.drop(columns='ANN', inplace=True)
            anndf.index = anndf.index.droplevel(1)
            df1 = df1.join(anndf, how='inner')
            del anndf
            del annhead
            del annheaderlist
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
            df1['sorter'] = df1['Annotation'].str.split('&').str[0].replace(IMPACT_SEVERITY)
            df1.loc[df1['HGVS.c'].str.contains('null'), 'HGVS.c'] = None
            df1['sorter2'] = [x[0] == x[1] for x in zip(df1['ALT'], df1['Allele'])]
            df1 = df1.sort_values(by=['CHROM', 'POS', 'sorter2', 'sorter'],
                                  ascending=[True, True, False, True]).drop_duplicates(['CHROM', 'POS', 'REF', 'ALT'])
            df1.drop(columns=['sorter', 'sorter2'], inplace=True)
            del IMPACT_SEVERITY
        df1.columns = df1.columns.str.upper()
        if 'HOM' in df1.columns:
            df1['HOM'] = df1['HOM'].map({True: 'HOM', False: 'HET'})
            df1.drop(columns='HET', inplace=True)
            df1.rename(columns={'HOM': 'ZIGOSITY'}, inplace=True)

        if 'ESP6500_MAF' in df1.columns:
            df1[['ESP6500_MAF_EA', 'ESP6500_MAF_AA', 'ESP6500_MAF_ALL']] = df1['ESP6500_MAF'].str.split(',',
                                                                                                        expand=True)
            df1['ESP6500_MAF_EA'] = df1['ESP6500_MAF_EA'].apply(divide, args=(100,))
            df1['ESP6500_MAF_AA'] = df1['ESP6500_MAF_AA'].apply(divide, args=(100,))
            df1['ESP6500_MAF_ALL'] = df1['ESP6500_MAF_ALL'].apply(divide, args=(100,))
            df1.drop(columns=['ESP6500_MAF'], inplace=True)
        if 'ESP6500_PH' in df1.columns:
            df1[['POLYPHEN_PRED', 'POLYPHEN_SCORE']] = df1['ESP6500_PH'].str.split(':', 1, expand=True)
            df1['POLYPHEN_PRED'] = df1['POLYPHEN_PRED'].str.strip('.').str.strip('.,')
            df1['POLYPHEN_SCORE'] = df1['POLYPHEN_SCORE'].str.split(',').str[0]
            df1.drop(columns=['ESP6500_PH'], inplace=True)
            df1.rename(columns={'ANNOTATION': 'EFFECT', 'ANNOTATION_IMPACT': 'IMPACT', 'ID': 'RSID'},
                       inplace=True)
        numcols = list()
        for x in pVCF.header_iter():
            if x.type == 'INFO':
                if x['Type'] in ['Float', 'Integer']:
                    numcols.append(x['ID'])
        numcols += ['ESP6500_MAF_EA', 'ESP6500_MAF_AA', 'ESP6500_MAF_ALL']
        numcols = list(set([x.upper() for x in numcols for y in df1.columns if x.upper() == y]))
        df1[numcols] = df1[numcols].apply(pd.to_numeric, errors='coerce', axis=1)
        df1 = df1.round(6)

        if ('CLINVAR_CLNSIG' in df1.columns):
            clinvartranslation = {'255': 'other', '0': 'Uncertain significance', '1': 'not provided', '2': 'Benign',
                                  '3': 'Likely Benign', '4': 'Likely pathogenic', '5': 'Pathogenic',
                                  '6': 'drug response',
                                  '7': 'histocompatibility'}
            for k, v in clinvartranslation.items():
                df1['CLINVAR_CLNSIG'] = df1['CLINVAR_CLNSIG'].str.replace(k, v)

        df1['AMINOCHANGE'] = df1['HGVS.P'].apply(aminoChange)
        df1.replace(['nan', '', np.nan], '.', inplace=True)
        df1.replace([[None], '.'], inplace=True, regex=True)
        df1 = df1.astype('str')
        df1['POS'] = df1['POS'].astype(int)
        df1 = df1.pipe(ParsedVCF)
        df1.name = name
        return df1

    @classmethod
    def mp_parser(cls, *vcfs, cores=int(cfg['GENERAL']['cores'])):
        if len(vcfs) < 1:
            logger.error('No vcfs provided!')
            exit(1)
        elif len(vcfs) == 1:
            pvcfs = list()
            pvcfs.append(ParsedVCF.from_vcf(vcfs[0]))
        else:
            try:
                [x + '' for x in vcfs]
            except:
                logger.error('All mp_parser args must be strings')
            else:
                logger.info('Starting Multi-Parser')
                if cores is None:
                    cores = mp.cpu_count()
                if cores > 1:
                    if len(vcfs) <= cores - 1:
                        pool = mp.Pool(processes=len(vcfs))
                    else:
                        pool = mp.Pool(processes=cores - 1)
                else:
                    pool = mp.Pool(processes=cores)
                pvcfs = pool.map(cls.from_vcf, (x for x in vcfs))
                pool.close()
                pool.join()
        return pvcfs

    def to_macrogen_xls(self, outpath):
        macrogen_cols = ['CHROM', 'POS', 'REF', 'ALT', 'DP', 'AD', 'QUAL', 'MQ', 'Zygosity', 'FILTER', 'Effect',
                         'IMPACT', 'Gene_Name', 'Feature_Type', 'Feature_ID', 'Transcript_BioType',
                         'Rank/Total', 'HGVS.c', 'HGVS.p', 'REF_AA', 'ALT_AA', 'cDNA_pos/cDNA_length',
                         'CDS_pos/CDS_length', 'AA_pos/AA_length', 'Distance', 'dbSNP142_ID', '1000Gp3_AF',
                         '1000Gp3_AFR_AF', '1000Gp3_AMR_AF', '1000Gp3_EAS_AF', '1000Gp3_EUR_AF', '1000Gp3_SAS_AF',
                         'ESP6500_MAF_EA', 'ESP6500_MAF_AA', 'ESP6500_MAF_ALL', 'SIFT_score', 'SIFT_pred',
                         'Polyphen2_HDIV_score', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred',
                         'CLINVAR_CLNSIG', 'CLINVAR_CLNDSDB', 'CLINVAR_CLNDSDBID', 'CLINVAR_CLNDBN',
                         'CLINVAR_CLNREVSTAT', 'CLINVAR_CLNACC']
        self['chrsort'] = self['CHROM'].replace({'X': 30, 'Y': 40})
        df1 = self.sort_values(['chrsort', 'POS'])[
            [x.upper() for x in macrogen_cols if x.upper() in self.columns]].copy()
        df1.to_excel(outpath)

    def panel(self, panel):
        logger.info('Analyzing Panel')
        try:
            pldf = pd.ExcelFile(panel).parse('GeneList')
            geneSymbolList = list(pldf.GeneSymbol.unique())
        except:
            logger.error('There was an error parsing GeneList')
            logger.debug('', exc_info=True)
            exit(1)
        panel_df = pd.DataFrame()
        for gene in geneSymbolList:
            panel_df = panel_df.append(self.loc[self['GENE_NAME'] == gene])
        panel_df = panel_df.pipe(ParsedVCF)
        panel_df.name = self.name
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
        # chequeo si el segundo es un vcf parseado o una ruta a un vcf
        if isinstance(vcf2, str):
            pvcf2 = ParsedVCF.from_vcf(vcf2)
        elif isinstance(vcf2, ParsedVCF):
            pvcf2 = vcf2

        indcols = ['QUAL', 'FILTER', 'DP', 'FS', 'MQ', 'SOR', 'QD', 'SET', 'BASEQRANKSUM', 'CLIPPINGRANKSUM',
                   'MQRANKSUM', 'READPOSRANKSUM', 'AC', 'SAMPLES_AF', 'MLEAC', 'MLEAF', 'DBSNPBUILDID']

        # chequeo si alguno es un duos y dropeo columnas individuales
        if ('VENN' in self.columns) & ('VENN' in pvcf2.columns):
            logger.error(
                'Only one argument can be a duos. Function cannot currently combine more than one duos (resulting'
                'in a trios). Both arguments provided where duos or trios')
            exit(1)
        elif 'VENN' in self.columns:
            logger.info('Running TRIOS analysis on %s' % ':'.join([self.name, pvcf2.name]))
            indicator = 'TRIOS'
            left = pvcf2
            right = self
            left.drop(columns=indcols, inplace=True)
        elif 'VENN' in pvcf2.columns:
            logger.info('Running TRIOS analysis on %s' % ':'.join([self.name, pvcf2.name]))
            indicator = 'TRIOS'
            left = self
            right = pvcf2
            left.drop(columns=indcols, inplace=True)
        else:
            logger.info('Running DUOS analysis on %s' % ':'.join([self.name, pvcf2.name]))
            indicator = 'DUOS'
            left = self
            right = pvcf2
            left.drop(columns=indcols, inplace=True)
            right.drop(columns=indcols, inplace=True)

        # Hago el merge
        mergedVCF = left.merge(right, on=['CHROM', 'POS', 'REF', 'ALT'], how='outer',
                               suffixes=('_' + self.name, '_' + pvcf2.name),
                               indicator=indicator)

        # columnas que deberían ser iguales y columnas que podrían ser distintas
        difcols = [x.replace('_' + self.name, '') for x in mergedVCF.columns if '_' + self.name in x]
        eqcols = [x for x in difcols if any(y in x for y in ['1000GP3', 'CLINVAR', 'ESP6500', 'RSID', 'POLYPHEN'])]
        difcols = set(difcols) - set(eqcols)
        # columnas a eliminar
        dropcols = [x + '_' + self.name for x in difcols]
        dropcols += [x + '_' + pvcf2.name for x in difcols]
        dropcols += [x + '_' + pvcf2.name for x in eqcols]
        dropcols += [x + '_' + self.name for x in eqcols]
        # armo un dataframe para combinar las columnas
        tmp = mergedVCF.loc[mergedVCF[indicator] == 'both'][dropcols]
        # combino las que deberían ser iguales
        for col in eqcols:
            mergedVCF[col] = mergedVCF[col + '_' + self.name].combine_first(mergedVCF[col + '_' + pvcf2.name])
        # combino las que podrían ser diferentes, si son diferentes (para variantes en ambos archivos, no las combino.
        for col in difcols:
            if all(tmp[col + '_' + self.name] == tmp[col + '_' + pvcf2.name]):
                mergedVCF[col] = mergedVCF[col + '_' + self.name].combine_first(mergedVCF[col + '_' + pvcf2.name])
            else:
                dropcols.remove(col + '_' + self.name)
                dropcols.remove(col + '_' + pvcf2.name)
        # dropeo las columans que combine
        mergedVCF.drop(columns=dropcols, inplace=True)
        # armo la columna indicadora para Duos y Trios
        if indicator == 'DUOS':
            mergedVCF['DUOS'].replace(
                {'left_only': self.name, 'right_only': pvcf2.name, 'both': self.name + ':' + pvcf2.name}, inplace=True)
            mergedVCF.rename(columns={'DUOS': 'VENN'}, inplace=True)
        if indicator == 'TRIOS':
            mergedVCF['PATIENT'] = None
            mergedVCF['PATIENT'] = np.where(mergedVCF['TRIOS'] == 'left_only', left.name, mergedVCF['PATIENT'])
            mergedVCF['PATIENT'] = np.where((mergedVCF['TRIOS'] == 'both'), mergedVCF['VENN'] + ':' + left.name,
                                            mergedVCF['PATIENT'])
            mergedVCF['PATIENT'] = np.where((mergedVCF['TRIOS'] == 'right_only'), mergedVCF['VENN'],
                                            mergedVCF['PATIENT'])
            mergedVCF.drop(columns=['VENN', 'TRIOS'], inplace=True)
            mergedVCF.rename(columns={'PATIENT': 'VENN', 'ZIGOSITY': 'ZIGOSITY_' + left.name}, inplace=True)
        mergedVCF.fillna('.', inplace=True)
        mergedVCF = mergedVCF.pipe(ParsedVCF)
        mergedVCF.name = ':'.join([self.name, pvcf2.name])
        return mergedVCF
