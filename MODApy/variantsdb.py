import glob
import logging
import os

import cyvcf2
import numpy as np
import pandas as pd

from MODApy.cfg import variantsDBPath, patientPath, cfg, resultsPath
from MODApy.vcfmgr import ParsedVCF

logger = logging.getLogger(__name__)


# TODO: FIND A MORE EFFICIENT WAY TO SUM EMPTY

class VariantsDB(pd.DataFrame):
    @property
    def _constructor(self):
        return VariantsDB

    @classmethod
    def from_exceldb(cls, excelpath):
        if os.path.exists(excelpath):
            try:
                db = pd.read_excel(excelpath)

            except:
                logger.error('There was an error parsing excel File')
                logger.debug('', exc_info=True)
                exit(1)
        else:
            logger.error('Path to excel file incorrect.')
            exit(1)
        db.set_index(['CHROM', 'POS', 'REF', 'ALT', 'GENE_NAME',
                      'HGVS.C', 'HGVS.P'], inplace=True)
        db = db.pipe(VariantsDB)
        return db

    @classmethod
    def from_csvdb(cls, csvpath):
        if os.path.exists(csvpath):
            try:
                files = [f for f in glob.glob(csvpath + '/*.csv')]
                if len(files) == 0:
                    return None
                dfs = [pd.read_csv(x) for x in files]
                db = pd.concat(dfs, sort=True)
                del files
                del dfs
            except:
                logger.error('There was an error parsing CSV File')
                logger.debug('', exc_info=True)
                exit(1)
        else:
            logger.error('Path to CSV file incorrect.')
            exit(1)
        db.set_index(['CHROM', 'POS', 'REF', 'ALT', 'GENE_NAME',
                      'HGVS.C', 'HGVS.P'], inplace=True)
        db = db.pipe(VariantsDB)
        return db

    @classmethod
    def buildDB(cls):
        def patientLister(db=None):
            vcfspath = []
            for dirpath, dirnames, filenames in os.walk(patientPath):
                for filename in [f for f in filenames if f.lower().endswith('final.vcf')]:
                    vcfspath.append(os.path.join(dirpath, filename))
            try:
                vcfsnames = [cyvcf2.Reader(x).samples[0] for x in vcfspath]
            except:
                logger.info(
                    'No Sample name in one of the vcfs files. Using File Names Instead')
                vcfsnames = [x.rsplit('/', maxsplit=1)
                             [-1].strip('.final.vcf') for x in vcfspath]

            if db is not None:
                addpatnames = [x for x in vcfsnames if x not in db.columns]
                if len(addpatnames) >= 1:
                    logger.info('Adding Patients: {}'.format(
                        [x for x in addpatnames]))
                else:
                    logger.error('No Patients to Add')
                    exit(1)
                patientslist = [
                    x for x in vcfspath for y in addpatnames if y in x]
            else:
                patientslist = vcfspath

            return patientslist

        def dbbuilder(patientslist, db=None):
            logger.info('Parsing Patients')
            pvcfs = ParsedVCF.mp_parser(*patientslist)
            pvcfs = [x[['CHROM', 'POS', 'REF', 'ALT', 'ZIGOSITY',
                        'GENE_NAME', 'HGVS.C', 'HGVS.P']] for x in pvcfs]
            for df in pvcfs:
                if 'ZIGOSITY' not in df.columns:
                    df['ZIGOSITY'] = 'UNKWN'
            pvcfs = [x.rename(columns={'ZIGOSITY': x.name})
                     for x in pvcfs if 'ZIGOSITY' in x.columns]
            logger.info('Merging parsed patients toDB')
            if db is not None:
                db = db.reset_index()
                pvcfs.insert(0, db)
            pvcfs = [x.set_index(
                ['CHROM', 'POS', 'REF', 'ALT', 'GENE_NAME', 'HGVS.C', 'HGVS.P']) for x in pvcfs]
            tempdb1 = pd.concat(pvcfs, axis=1, join='outer')
            tempdb2 = tempdb1.reset_index().groupby(['CHROM', 'POS', 'REF', 'ALT']).agg(
                {'GENE_NAME': ' | '.join, 'HGVS.P': ' | '.join, 'HGVS.C': ' | '.join}).reset_index()
            pvcfs2 = [x.reset_index().drop(
                columns=['GENE_NAME', 'HGVS.C', 'HGVS.P']) for x in pvcfs]
            pvcfs2.insert(0, tempdb2)
            pvcfs2 = [x.set_index(['CHROM', 'POS', 'REF', 'ALT'])
                      for x in pvcfs2]
            db = pd.concat(pvcfs2, axis=1, join='outer')
            colslist = ['GENE_NAME', 'HGVS.C', 'HGVS.P']
            for col in colslist:
                db[col] = db[col].apply(
                    lambda x: ' | '.join(set(x.split(' | '))))
            db = db.reset_index().set_index(
                ['CHROM', 'POS', 'REF', 'ALT', 'GENE_NAME', 'HGVS.C', 'HGVS.P'])
            db.replace({'.': np.nan}, inplace=True)
            db = db.pipe(VariantsDB)
            db = db.calcfreqs()
            return db

        try:
            logger.info('Checking DB File')
            if variantsDBPath.rsplit('.')[-1].lower() == 'xlsx':
                db = VariantsDB.from_exceldb(
                    variantsDBPath.rsplit('/', maxsplit=1)[0])
                patientslist = patientLister(db)
            elif variantsDBPath.rsplit('.')[-1].lower() == 'csv':
                db = VariantsDB.from_csvdb(
                    variantsDBPath.rsplit('/', maxsplit=1)[0])
                patientslist = patientLister(db)
            else:
                logger.error('VariantsDBPath must be a xlsx or csv file')
                exit(1)
        except:
            exit()
            logger.info('No DB Found, Building new Variants DB')
            patientslist = patientLister()
            db = None
        sublists = [patientslist[i:i + int(cfg['GENERAL']['cores'])] for i in
                    range(0, len(patientslist), int(cfg['GENERAL']['cores']))]
        for l in sublists:
            db = dbbuilder(l, db)
        return db

    def addPatientToDB(self, patient):
        if patient.rsplit('/')[-1].strip('.final.vcf') in self.columns:
            logger.error('Patient already is in DB')
            exit(1)
        if isinstance(patient, str):
            pvcf = ParsedVCF.from_vcf(patient)
        elif isinstance(patient, ParsedVCF):
            pvcf = patient
        else:
            logger.error(
                'Patient must be either a path to vcf or a ParsedVCF object')
            logger.debug('', exc_info=True)
            exit(1)
        pvcf = pvcf[['CHROM', 'POS', 'REF', 'ALT',
                     'ZIGOSITY', 'GENE_NAME', 'HGVS.C', 'HGVS.P']]
        if 'ZIGOSITY' not in pvcf.columns:
            pvcf['ZIGOSITY'] = 'UNKWN'
        pvcf.rename(columns={'ZIGOSITY': pvcf.name}, inplace=True)
        pvcf.set_index(['CHROM', 'POS', 'REF', 'ALT',
                        'GENE_NAME', 'HGVS.C', 'HGVS.P'], inplace=True)
        db = pd.concat([self, pvcf], axis=1, join='outer')
        db = db.pipe(VariantsDB)
        db = db.calcfreqs()
        return db

    def to_VarDBXLS(self):
        logger.info('Writing DB to Excel')
        self.reset_index(inplace=True)
        self['POS'] = self['POS'].astype(int)
        self.sort_values(['CHROM', 'POS'], inplace=True)
        os.makedirs(variantsDBPath.rsplit('/', maxsplit=1)[0], exist_ok=True)
        vdbpath = variantsDBPath.rsplit('.', maxsplit=1)[0]
        output = pd.ExcelWriter(variantsDBPath)
        workbook = output.book
        datasheet = workbook.add_worksheet('VariantSDB')
        output.sheets['VariantsDB'] = datasheet
        formatpos = workbook.add_format({'num_format': '###,###,###'})
        self['POS'] = self['POS'].astype(int)
        datasheet.set_column('B:B', 15, formatpos)
        for chrom in self['CHROM'].unique():
            self[self['CHROM'] == chrom].to_excel(vdbpath + str(chrom) + '.csv', index=False, float_format='%.5f',
                                                  merge_cells=False)
        output.save()
        logger.info('Xlsx DB construction complete')

    def to_VarDBCSV(self):
        logger.info('Writing DB to CSV')
        self.reset_index(inplace=True)
        self['POS'] = self['POS'].astype(int)
        self.sort_values(['CHROM', 'POS'], inplace=True)
        vdbpath = variantsDBPath.rsplit('.', maxsplit=1)[0]
        os.makedirs(variantsDBPath.rsplit('/', maxsplit=1)[0], exist_ok=True)
        for chrom in self['CHROM'].unique():
            self[self['CHROM'] == chrom].to_csv(
                vdbpath + str(chrom) + '.csv', index=False, float_format='%.5f')
        logger.info('DB construction complete')

    def calcfreqs(self):
        logger.info('Calculating Variant Frequencies')
        patients = self.columns.tolist()
        if 'FREQ' in patients:
            patients.remove('FREQ')
        if 'ALLELE_FREQ' in patients:
            patients.remove('ALLELE_FREQ')
        self.replace({'.': np.nan}, inplace=True)
        self['FREQ'] = (self[patients].notnull().sum(axis=1) / len(patients))
        self['ALLELE_FREQ'] = self[patients].apply(
            lambda x: ((x.str.contains('HOM') * 2 + x.str.contains('HET') * 1).sum()) / len(patients * 2), axis=1)
        cols = self.columns.tolist()
        cols.remove('FREQ')
        cols.remove('ALLELE_FREQ')
        self = self[['ALLELE_FREQ', 'FREQ'] + cols]
        self.replace({np.nan: '.'}, inplace=True)
        self.pipe(VariantsDB)
        return self

    def annotate_excel(self, df, fileName):
        logger.info('Annotating Excel file')
        cols_to_drop = ['FREQ', 'ALLELE_FREQ', 'VARDB_FREQ']
        df.drop(
            columns=[x for x in cols_to_drop if x in df.columns], inplace=True)
        df = df.merge(self.reset_index()[['CHROM', 'POS', 'REF', 'ALT', 'FREQ', 'ALLELE_FREQ']],
                      on=['CHROM', 'POS', 'REF', 'ALT'], how='left')
        df.rename(columns={'FREQ': 'VARDB_FREQ'}, inplace=True)
        df['VARDB_FREQ'] = pd.to_numeric(df['VARDB_FREQ'], errors='coerce')
        df['VARDB_FREQ'].round(6)
        outpath = resultsPath + 'vDBannotated' + '/' + \
                  fileName.rsplit(
                      '.', maxsplit=1)[0].replace('.annotated', '') + '.annotated.xlsx'
        logger.info(outpath)
        firstcols = ['GENE_NAME', 'AMINOCHANGE', 'HGVS.P', 'HGVS.C', 'RSID', 'IMPACT', 'EFFECT', 'VARDB_FREQ',
                     'ALLELE_FREQ']
        lastcols = [x for x in df.columns if x not in firstcols]
        output = pd.ExcelWriter(outpath)
        workbook = output.book
        datasheet = workbook.add_worksheet('DATA')
        statsheet = workbook.add_worksheet('STATISTICS')
        output.sheets['STATISTICS'] = statsheet

        output.sheets['DATA'] = datasheet
        formatnum = workbook.add_format({'num_format': '0.00000'})
        # for i, col in enumerate(self.columns):
        datasheet.set_column(0, len(df.columns), 15, formatnum)

        formatpos = workbook.add_format({'num_format': '###,###,###'})
        datasheet.set_column(df.columns.to_list().index(
            'POS'), df.columns.to_list().index('POS'), 15, formatpos)
        datasheet.set_column(df.columns.to_list().index(
            'RSID'), df.columns.to_list().index('RSID'), 15)
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
        datasheet.conditional_format(0, df.columns.to_list().index('IMPACT'),
                                     len(df), df.columns.to_list().index(
                'IMPACT'),
                                     {'type': 'text', 'criteria': 'containing', 'value': 'HIGH',
                                      'format': highformat, })
        datasheet.conditional_format(0, df.columns.to_list().index('IMPACT'),
                                     len(df), df.columns.to_list().index(
                'IMPACT'),
                                     {'type': 'text', 'criteria': 'containing', 'value': 'MODIFIER',
                                      'format': modformat})
        datasheet.conditional_format(0, df.columns.to_list().index('IMPACT'),
                                     len(df), df.columns.to_list().index(
                'IMPACT'),
                                     {'type': 'text', 'criteria': 'containing', 'value': 'MODERATE',
                                      'format': moderformat})
        datasheet.conditional_format(0, df.columns.to_list().index('IMPACT'),
                                     len(df), df.columns.to_list().index(
                'IMPACT'),
                                     {'type': 'text', 'criteria': 'containing', 'value': 'LOW', 'format': lowformat})
        logger.info('Writing Excel File')
        df[firstcols + lastcols].to_excel(output, sheet_name='DATA',
                                          merge_cells=False, index=False, header=True)

        if (df.reset_index().index.max() < 32150):
            logger.info('Redirecting IDs and GENEs to URLs')
            try:
                colid = df.columns.to_list().index('RSID')
                colgen = df.columns.to_list().index('GENE_NAME')
                row = 2
                for x in zip(df['RSID'], df['GENE_NAME']):
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
        stats = ParsedVCF.general_stats(df)
        stats.to_excel(output, sheet_name='STATISTICS')
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
        datasheet.autofilter(0, 0, len(self), len(df.columns))
        output.save()
        logger.info('File saved to %s' % outpath)
