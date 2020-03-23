from MODApy.cfg import cfg
import pandas as pd
import numpy as np
import cyvcf2
from collections import OrderedDict
import logging
import multiprocessing as mp
import os

import matplotlib

matplotlib.use('agg')
import matplotlib_venn as venn
import matplotlib.pyplot as plt


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
        if(len(splitdf)>0):
            ALT = splitdf['ALT'].astype(str).str.split(
                ',', n=1, expand=True).stack().rename('ALT')
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
        splitlist = ['ID', 'AC', 'AF', 'SAMPLES_AF',
                     'MLEAC', 'MLEAF', 'VARTYPE', 'dbSNPBuildID']
        splitlist = [x for x in splitlist if x in df1.columns]
        splitlist += [x for x in df1.columns if x.startswith(
            ('1000', 'CLINVAR'))]
        for col in splitlist:
            odd[col] = odd[col].astype(str).str.split(',', n=1).str[0]
            even[col] = even[col].apply(lambda x:
                                        x if len(str(x).split(',')) <= 1 else str(x).split(',', maxsplit=1)[1])
        splitdf = pd.concat([odd, even]).sort_index().replace(
            to_replace=['\(', '\)'], value='', regex=True)
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
            annhead = pVCF.get_header_type('ANN')['Description'].strip(
                '"Functional annotations: \'"')
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
            df1['sorter'] = df1['Annotation'].str.split(
                '&').str[0].replace(IMPACT_SEVERITY)
            df1.loc[df1['HGVS.c'].str.contains('null'), 'HGVS.c'] = None
            df1['sorter2'] = [x[0] == x[1]
                              for x in zip(df1['ALT'], df1['Allele'])]
            df1 = df1.sort_values(by=['CHROM', 'POS', 'sorter2', 'sorter'],
                                  ascending=[True, True, False, True]).drop_duplicates(['CHROM', 'POS', 'REF', 'ALT'])
            df1.drop(columns=['sorter', 'sorter2'], inplace=True)
            del IMPACT_SEVERITY
        df1.columns = df1.columns.str.upper()
        if 'HGVS.P' in df1.columns:
            df1['AMINOCHANGE'] = df1['HGVS.P'].apply(aminoChange)
        if 'HOM' in df1.columns:
            df1['HOM'] = df1['HOM'].replace(
                {True: 'HOM', np.nan: 'HET', None: 'HET'})
            df1.drop(columns='HET', inplace=True)
            df1.rename(columns={'HOM': 'ZIGOSITY'}, inplace=True)
        if 'ESP6500_MAF' in df1.columns:
            df1[['ESP6500_MAF_EA', 'ESP6500_MAF_AA', 'ESP6500_MAF_ALL']] = df1['ESP6500_MAF'].str.split(',',
                                                                                                        expand=True)
            df1['ESP6500_MAF_EA'] = df1['ESP6500_MAF_EA'].apply(
                divide, args=(100,))
            df1['ESP6500_MAF_AA'] = df1['ESP6500_MAF_AA'].apply(
                divide, args=(100,))
            df1['ESP6500_MAF_ALL'] = df1['ESP6500_MAF_ALL'].apply(
                divide, args=(100,))
            df1.drop(columns=['ESP6500_MAF'], inplace=True)
        if 'ESP6500_PH' in df1.columns:
            df1[['POLYPHEN_PRED', 'POLYPHEN_SCORE']
                ] = df1['ESP6500_PH'].str.split(':', 1, expand=True)
            df1['POLYPHEN_PRED'] = df1['POLYPHEN_PRED'].str.strip(
                '.').str.strip('.,')
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
        numcols = list(
            set([x.upper() for x in numcols for y in df1.columns if x.upper() == y]))
        df1[numcols] = df1[numcols].apply(
            pd.to_numeric, errors='coerce', axis=1)
        df1 = df1.round(6)

        if ('CLINVAR_CLNSIG' in df1.columns):
            clinvartranslation = {'255': 'other', '0': 'Uncertain significance', '1': 'not provided', '2': 'Benign',
                                  '3': 'Likely Benign', '4': 'Likely pathogenic', '5': 'Pathogenic',
                                  '6': 'drug response',
                                  '7': 'histocompatibility'}
            for k, v in clinvartranslation.items():
                df1['CLINVAR_CLNSIG'] = df1['CLINVAR_CLNSIG'].str.replace(k, v)

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
        if len(panel_df) < 1:
            logger.error(
                'After running Panel, resulting Dataframe holds no variants')
        return panel_df

    def duos(self, vcf2, VENNPLACE=None):
        """
        Method to compare two vcfs, using CHROM, POS, REF and ALT columns as index.
        Parameters
        ----------
        vcf2
            VCF file to compare to.

        Returns a Dataframe containing a new column 'DUOS', that indicates in which file is the variant.
        """

        # chequeo si el segundo es un vcf parseado o una ruta a un vcf
        def _trios_stats(self, names):
            logger.info('Calculating Trios statistics')
            trios = self.groupby('VENN', sort=False).size()
            A, B, C = names.split(':')
            names = [A, B]
            names.append(A + ':' + B)
            names.append(C)
            names.append(A + ':' + C)
            names.append(B + ':' + C)
            names.append(':'.join([A, B, C]))
            trios = trios.reindex(names).fillna(0)
            trios = trios.astype(int)
            triosgraph = plt.figure()
            venn.venn3(trios, set_labels=[A, B, C], set_colors=['b', 'r', 'g'])
            triosgraph.savefig('./venn.png', dpi=triosgraph.dpi)
            triosgraph.clf()

        def _duos_stats(self, names):
            logger.info('Calculating Duos statistics')
            duos = self.groupby('VENN', sort=False).size()
            A, B = names.split(':')
            names = [A, B, names]
            duos = duos.reindex(names).fillna(0)
            duosgraph = plt.figure()
            venn.venn2(duos, set_labels=[A, B], set_colors=['b', 'r'])
            duosgraph.savefig('./venn.png', dpi=duosgraph.dpi)
            duosgraph.clf()

        if isinstance(vcf2, str):
            pvcf2 = ParsedVCF.from_vcf(vcf2)
        elif isinstance(vcf2, ParsedVCF):
            pvcf2 = vcf2

        indcols = ['QUAL', 'FILTER', 'DP', 'FS', 'MQ', 'SOR', 'QD', 'SET', 'BASEQRANKSUM', 'CLIPPINGRANKSUM',
                   'MQRANKSUM', 'READPOSRANKSUM', 'AC', 'SAMPLES_AF', 'MLEAC', 'MLEAF', 'DBSNPBUILDID']
        indself = [x for x in indcols if x in self.columns]
        indpvcf2 = [x for x in indcols if x in pvcf2.columns]
        self.drop(columns=indself, inplace=True)
        pvcf2.drop(columns=indpvcf2, inplace=True)
        del indcols, indself, indpvcf2
        # chequeo si alguno es un duos y dropeo columnas individuales
        if ('VENN' in self.columns) & ('VENN' in pvcf2.columns):
            logger.error(
                'Only one argument can be a duos. Function cannot currently combine more than one duos (resulting'
                'in a trios). Both arguments provided where duos or trios')
            exit(1)
        elif 'VENN' in self.columns:
            logger.info('Running TRIOS analysis on %s' %
                        ':'.join([self.name, pvcf2.name]))
            indicator = 'TRIOS'
            left = pvcf2
            right = self
        elif 'VENN' in pvcf2.columns:
            logger.info('Running TRIOS analysis on %s' %
                        ':'.join([self.name, pvcf2.name]))
            indicator = 'TRIOS'
            left = self
            right = pvcf2
        else:
            logger.info('Running DUOS analysis on %s' %
                        ':'.join([self.name, pvcf2.name]))
            indicator = 'DUOS'
            left = self
            right = pvcf2

        # Hago el merge
        mergedVCF = left.merge(right, on=['CHROM', 'POS', 'REF', 'ALT'], how='outer',
                               suffixes=('_' + self.name, '_' + pvcf2.name),
                               indicator=indicator)

        # columnas que deberían ser iguales y columnas que podrían ser distintas
        difcols = [x.replace('_' + self.name, '')
                   for x in mergedVCF.columns if '_' + self.name in x]
        eqcols = [x for x in difcols if any(
            y in x for y in ['1000GP3', 'CLINVAR', 'ESP6500', 'RSID', 'POLYPHEN'])]
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
            mergedVCF[col] = mergedVCF[col + '_' +
                                       self.name].combine_first(mergedVCF[col + '_' + pvcf2.name])
        # combino las que podrían ser diferentes, si son diferentes (para variantes en ambos archivos, no las combino.
        for col in difcols:
            if all(tmp[col + '_' + self.name] == tmp[col + '_' + pvcf2.name]):
                mergedVCF[col] = mergedVCF[col + '_' +
                                           self.name].combine_first(mergedVCF[col + '_' + pvcf2.name])
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
            names = ':'.join([self.name, pvcf2.name])
            _duos_stats(mergedVCF, names)
            if VENNPLACE is not None:
                print(VENNPLACE)
                if VENNPLACE == 'A':
                    mergedVCF = mergedVCF[mergedVCF['VENN'] == self.name]
                elif VENNPLACE == 'B':
                    mergedVCF = mergedVCF[mergedVCF['VENN'] == pvcf2.name]
                elif VENNPLACE == 'A:B':
                    mergedVCF = mergedVCF[mergedVCF['VENN']
                                          == self.name + ':' + pvcf2.name]
                else:
                    logger.error(
                        'VENNPLACE can only be A, B or A:B for a trios analysis')
                    logger.debug('', exc_info=True)
                    exit(1)
        if indicator == 'TRIOS':
            mergedVCF['PATIENT'] = None
            mergedVCF['PATIENT'] = np.where(
                mergedVCF['TRIOS'] == 'left_only', left.name, mergedVCF['PATIENT'])
            mergedVCF['PATIENT'] = np.where((mergedVCF['TRIOS'] == 'both'), mergedVCF['VENN'] + ':' + left.name,
                                            mergedVCF['PATIENT'])
            mergedVCF['PATIENT'] = np.where((mergedVCF['TRIOS'] == 'right_only'), mergedVCF['VENN'],
                                            mergedVCF['PATIENT'])
            mergedVCF.drop(columns=['VENN', 'TRIOS'], inplace=True)
            mergedVCF.rename(
                columns={'PATIENT': 'VENN', 'ZIGOSITY': 'ZIGOSITY_' + left.name}, inplace=True)
            names = self.name + ':' + pvcf2.name
            _trios_stats(mergedVCF, names)

            if VENNPLACE is not None:
                names = names.split(':')
                if VENNPLACE == 'A':
                    mergedVCF = mergedVCF[
                        (mergedVCF['VENN'].str.contains(names[0])) & ~(mergedVCF['VENN'].str.contains(names[1])) & ~(
                            mergedVCF['VENN'].str.contains(names[2]))]
                elif VENNPLACE == 'B':
                    mergedVCF = mergedVCF[
                        ~(mergedVCF['VENN'].str.contains(names[0])) & (mergedVCF['VENN'].str.contains(names[1])) & ~(
                            mergedVCF['VENN'].str.contains(names[2]))]
                elif VENNPLACE == 'C':
                    mergedVCF = mergedVCF[
                        ~(mergedVCF['VENN'].str.contains(names[0])) & ~(mergedVCF['VENN'].str.contains(names[1])) & (
                            mergedVCF['VENN'].str.contains(names[2]))]
                elif VENNPLACE == 'A:B':
                    mergedVCF = mergedVCF[
                        (mergedVCF['VENN'].str.contains(names[0])) & (mergedVCF['VENN'].str.contains(names[1])) & ~(
                            mergedVCF['VENN'].str.contains(names[2]))]
                elif VENNPLACE == 'A:C':
                    mergedVCF = mergedVCF[
                        (mergedVCF['VENN'].str.contains(names[0])) & ~(mergedVCF['VENN'].str.contains(names[1])) & (
                            mergedVCF['VENN'].str.contains(names[2]))]
                elif VENNPLACE == 'B:C':
                    mergedVCF = mergedVCF[
                        ~(mergedVCF['VENN'].str.contains(names[0])) & (mergedVCF['VENN'].str.contains(names[1])) & (
                            mergedVCF['VENN'].str.contains(names[2]))]
                elif VENNPLACE == 'A:B:C':
                    mergedVCF = mergedVCF[
                        (mergedVCF['VENN'].str.contains(names[0])) & (mergedVCF['VENN'].str.contains(names[1])) & (
                            mergedVCF['VENN'].str.contains(names[2]))]
                else:
                    logger.error(
                        'VENNPLACE can only be A, B, C, A:B, A:C, B:C or A:B:C in a trios analysis')
                    logger.debug('', exc_info=True)
                    exit(1)
        if len(mergedVCF) < 1:
            logger.error(
                'After running Duos/Trios, resulting Dataframe does not hold any variants')
            exit(1)
        mergedVCF.fillna('.', inplace=True)
        mergedVCF = mergedVCF.pipe(ParsedVCF)
        mergedVCF.name = ':'.join([self.name, pvcf2.name])
        return mergedVCF

    def general_stats(self):
        if 'VENN' in self.columns:
            colstats = ['CHROM', 'VARTYPE', 'IMPACT', 'EFFECT']
        else:
            colstats = ['CHROM', 'ZIGOSITY', 'VARTYPE', 'IMPACT', 'EFFECT']
        if set(colstats).issubset(self.columns):
            logger.info('Calculating General Statistics')
            vcfstats = self.groupby(colstats).size().to_frame(name='count')
            vcfstats.name = 'stats'
            plt.pie([item for sublist in vcfstats.groupby('CHROM').count().values for item in sublist],
                    labels=vcfstats.groupby('CHROM').size().index.values)
            my_circle = plt.Circle((0, 0), 0.7, color='white')
            chromVars = plt.gcf()
            chromVars.gca().add_artist(my_circle)
            chromVars.savefig('./general.png')
            return vcfstats

    def vcf_to_excel(self, outpath):
        # logger.info('Getting Frequencies from VariantsDB')
        # if os.path.exists(cfg['PATHS']['dbpath']):
        #     if cfg['PATHS']['dbpath'].rsplit('.')[-1].lower() == 'xlsx':
        #         variantsDB = pd.read_excel(cfg['PATHS']['dbpath'])
        #     elif cfg['PATHS']['dbpath'].rsplit('.')[-1].lower() == 'csv':
        #         variantsDB = pd.read_csv(cfg['PATHS']['dbpath'])
        #     else:
        #         logger.error('VariantsDBPath must be a xlsx or csv file')
        #         exit(1)
        # #self = self.merge(variantsDB[['CHROM', 'POS', 'REF', 'ALT', 'FREQ', 'ALLELE_FREQ']],
        #                   on=['CHROM', 'POS', 'REF', 'ALT'], how='left')
        # #self.rename(columns={'FREQ': 'VARDB_FREQ'}, inplace=True)
        # #self['VARDB_FREQ'] = pd.to_numeric(self['VARDB_FREQ'], errors='coerce')
        # #self['VARDB_FREQ'].round(6)
        # logger.info('Formating Excel File')
        os.makedirs(outpath.rsplit('/', maxsplit=1)[0], exist_ok=True)
        output = pd.ExcelWriter(outpath)
        cols_selected = cfg["OUTPUT"]["columnsorder"].replace(',', ' ').split()
        if 'VENN' in self.columns:
            if 'ZIGOSITY' in cols_selected:
                cols_selected += [x for x in self.columns if 'ZIGOSITY' in x]
        finalcols = [x for x in cols_selected if x in self.columns]
        self = self[finalcols]
        self = self.sort_values(by=finalcols[0])
        workbook = output.book
        datasheet = workbook.add_worksheet('DATA')
        statsheet = workbook.add_worksheet('STATISTICS')

        output.sheets['DATA'] = datasheet

        formatnum = workbook.add_format({'num_format': '0.00000'})
        # for i, col in enumerate(self.columns):
        datasheet.set_column(0, len(self.columns), 15, formatnum)

        formatpos = workbook.add_format({'num_format': '###,###,###'})
        datasheet.set_column(finalcols.index(
            'POS'), finalcols.index('POS'), 15, formatpos)
        datasheet.set_column(finalcols.index(
            'RSID'), finalcols.index('RSID'), 15)
        # Light red fill with dark red text.
        highformat = workbook.add_format(
            {'bg_color': '#FFC7CE', 'font_color': '#9C0006', 'bold': True})
        # Light yellow fill with dark yellow text.
        modformat = workbook.add_format(
            {'bg_color': '#FFFF99', 'font_color': '#9C6500', 'bold': True})
        # Light orange fill with dark orange text.
        moderformat = workbook.add_format(
            {'bg_color': '#FFCC99', 'font_color': '#FF6600', 'bold': True})
        # Green fill with dark green text.
        lowformat = workbook.add_format(
            {'bg_color': '#C6EFCE', 'font_color': '#006100', 'bold': True})
        datasheet.conditional_format(0, cols_selected.index('IMPACT'),
                                     len(self), cols_selected.index('IMPACT'),
                                     {'type': 'text', 'criteria': 'containing', 'value': 'HIGH',
                                      'format': highformat, })
        datasheet.conditional_format(0, cols_selected.index('IMPACT'),
                                     len(self), cols_selected.index('IMPACT'),
                                     {'type': 'text', 'criteria': 'containing', 'value': 'MODIFIER',
                                      'format': modformat})
        datasheet.conditional_format(0, cols_selected.index('IMPACT'),
                                     len(self), cols_selected.index('IMPACT'),
                                     {'type': 'text', 'criteria': 'containing', 'value': 'MODERATE',
                                      'format': moderformat})
        datasheet.conditional_format(0, cols_selected.index('IMPACT'),
                                     len(self), cols_selected.index('IMPACT'),
                                     {'type': 'text', 'criteria': 'containing', 'value': 'LOW', 'format': lowformat})
        logger.info('Writing Excel File')
        self.to_excel(output, sheet_name='DATA',
                      merge_cells=False, index=False, header=True)

        if (self.reset_index().index.max() < 32150):
            logger.info('Redirecting IDs and GENEs to URLs')
            try:
                colid = cols_selected.index('RSID')
                colgen = cols_selected.index('GENE_NAME')
                row = 2
                for x in zip(self['RSID'], self['GENE_NAME']):
                    if type(x[0]) == str:
                        urlrs = "https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=%s"
                        rsvalue = (x[0].replace(';', ',').split(','))[0]
                        datasheet.write_url('%s%i' % (chr(colid + 65), (row)),
                                            urlrs % rsvalue, string=rsvalue)
                    if type(x[1]) == str:
                        urlgen = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s"
                        datasheet.write_url('%s%i' % (chr(colgen + 65), (row)),
                                            urlgen % x[1], string=x[1])
                    row += 1
            except Exception as e:
                logger.error(e, exc_info=True)
        datasheet.autofilter(0, 0, len(self), len(finalcols))
        stats = self.general_stats()
        output.sheets['STATISTICS'] = statsheet
        try:
            stats.to_excel(output, sheet_name='STATISTICS')
        except Exception as e:
            logger.error(
                'Could not print statistics. Error was {}'.format(e), exc_info=True)
        try:
            statsheet.insert_image('H2', './general.png')
        except Exception as e:
            logger.error(
                'Could not print stats graphs. Error was {}'.format(e), exc_info=True)
        if os.path.isfile('./venn.png'):
            statsheet.insert_image('H25', './venn.png')
        output.save()
        try:
            os.remove('./general.png')
        except:
            logger.debug('Could not remove general.png')
        try:
            os.remove('./venn.png')
        except:
            logger.debug('Could not remove venn.png')
