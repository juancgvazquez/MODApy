import logging
import os

import numpy as np
import pandas as pd

from MODApy.cfg import variantsDBPath, patientPath, cfg
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
        db.set_index(['CHROM', 'POS', 'REF', 'ALT','GENE_NAME'], inplace=True)
        db = db.pipe(VariantsDB)
        return db

    @classmethod
    def from_csvdb(cls, csvpath):
        if os.path.exists(csvpath):
            try:
                db = pd.read_csv(csvpath)
            except:
                logger.error('There was an error parsing CSV File')
                logger.debug('', exc_info=True)
                exit(1)
        else:
            logger.error('Path to CSV file incorrect.')
            exit(1)
        db.set_index(['CHROM', 'POS', 'REF', 'ALT','GENE_NAME'], inplace=True)
        db = db.pipe(VariantsDB)
        return db

    @classmethod
    def buildDB(cls):
        def patientLister(db=None):
            vcfspath = []
            for dirpath, dirnames, filenames in os.walk(patientPath):
                for filename in [f for f in filenames if f.lower().endswith('final.vcf')]:
                    vcfspath.append(os.path.join(dirpath, filename))
            vcfsnames = [x.rsplit('/', maxsplit=1)[-1].strip('.final.vcf') for x in vcfspath]

            if db is not None:
                addpatnames = [x for x in vcfsnames if x not in db.columns]
                if len(addpatnames) >= 1:
                    logger.info('Adding Patients: {}'.format([x for x in addpatnames]))
                else:
                    logger.error('No Patients to Add')
                    exit(1)
                patientslist = [x for x in vcfspath for y in addpatnames if y in x]
            else:
                patientslist = vcfspath

            return patientslist

        def dbbuilder(patientslist, db=None):
            logger.info('Parsing Patients')
            pvcfs = ParsedVCF.mp_parser(*patientslist)
            pvcfs = [x[['CHROM', 'POS', 'REF', 'ALT', 'ZIGOSITY', 'GENE_NAME']] for x in pvcfs]
            for df in pvcfs:
                if 'ZIGOSITY' not in df.columns:
                    df['ZIGOSITY'] = 'UNKWN'
            pvcfs = [x.rename(columns={'ZIGOSITY': x.name}) for x in pvcfs if 'ZIGOSITY' in x.columns]
            pvcfs = [x.set_index(['CHROM', 'POS', 'REF', 'ALT','GENE_NAME']) for x in pvcfs]
            if db is not None:
                pvcfs.insert(0, db)
            logger.info('Merging parsed patients toDB')
            db = pd.concat(pvcfs, axis=1, join='outer')
            db.replace({'.': np.nan}, inplace=True)
            db = db.pipe(VariantsDB)
            db = db.calcfreqs()
            return db

        if os.path.exists(variantsDBPath):
            if variantsDBPath.rsplit('.')[-1].lower() == 'xlsx':
                logger.info('Parsing XLSX DB File')
                db = VariantsDB.from_exceldb(variantsDBPath)
                patientslist = patientLister(db)
            elif variantsDBPath.rsplit('.')[-1].lower() == 'csv':
                logger.info('Parsing CSV DB File')
                db = VariantsDB.from_csvdb(variantsDBPath)
                patientslist = patientLister(db)
            else:
                logger.error('VariantsDBPath must be a xlsx or csv file')
                exit(1)

        else:
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
            logger.error('Patient must be either a path to vcf or a ParsedVCF object')
            logger.debug('', exc_info=True)
            exit(1)
        pvcf = pvcf[['CHROM', 'POS', 'REF', 'ALT', 'ZIGOSITY', 'GENE_NAME']]
        if 'ZIGOSITY' not in pvcf.columns:
            pvcf['ZIGOSITY'] = 'UNKWN'
        pvcf.rename(columns={'ZIGOSITY': pvcf.name}, inplace=True)
        pvcf.set_index(['CHROM', 'POS', 'REF', 'ALT','GENE_NAME'], inplace=True)
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
        output = pd.ExcelWriter(variantsDBPath)
        workbook = output.book
        datasheet = workbook.add_worksheet('VariantSDB')
        output.sheets['VariantsDB'] = datasheet
        formatpos = workbook.add_format({'num_format': '###,###,###'})
        self['POS'] = self['POS'].astype(int)
        datasheet.set_column('B:B', 15, formatpos)
        self.to_excel(output, sheet_name='VariantsDB', index=False, merge_cells=False)
        output.save()
        logger.info('Xlsx DB construction complete')

    def to_VarDBCSV(self):
        logger.info('Writing DB to CSV')
        self.reset_index(inplace=True)
        self['POS'] = self['POS'].astype(int)
        self.sort_values(['CHROM', 'POS'], inplace=True)
        os.makedirs(variantsDBPath.rsplit('/', maxsplit=1)[0], exist_ok=True)
        self.to_csv(variantsDBPath, index=False)
        logger.info('DB construction complete')

    def calcfreqs(self):
        logger.info('Calculating Variant Frequencies')
        patients = self.columns.tolist()
        if 'FREQ' in patients:
            patients.remove('FREQ')
        self.replace({'.': np.nan}, inplace=True)
        self['FREQ'] = (self[patients].notnull().sum(axis=1) / len(patients))
        cols = self.columns.tolist()
        cols.remove('FREQ')
        self = self[['FREQ'] + cols]
        self.replace({np.nan: '.'}, inplace=True)
        self.pipe(VariantsDB)
        return self
