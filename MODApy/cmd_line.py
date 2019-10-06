#!/usr/bin/env python
import argparse
import logging
import os
import shlex
import subprocess
from sys import argv

import pandas as pd

from MODApy import cfg, pipeline, vcfmgr, downloader, variantsdb, coverage
from MODApy.version import __version__

logger = logging.getLogger(__name__)


class Parser(object):
    '''
    Main parser for command line arguments. Uses ArgParse
    '''

    def __init__(self):
        parser = argparse.ArgumentParser(description="Multi-Omics Data Analisis for Python", usage='''MODApy <command> [<args>]

        Commands:
        launcher        Run MoDAPy Web Interface
        variantsDB      Work with Variants Database
        addPatient      Download Patient Data to Patients folder. Receives both url or xls/xlsx
        pipeline        Run pipeline on FastQ file/s
        single          Run study on a single patient
        duos            Run Duos analysis on two selected patients
        trios           Run Trios analysis on three selected patients
        coverageStats   Generate coverages stats for bam file or list of files

        For more info on any of these commands, use "cmd_line.py <command> -h
        
        You can check the package version using -v or --version"''')

        parser.add_argument("command", help="Select command to run")
        parser.add_argument("-v", "--version", action='version', version='MODApy ' + __version__)
        # exclude all arguments but the first one
        args = parser.parse_args(argv[1:2])
        if not hasattr(self, args.command):
            logger.error('Unrecognized commands')
            parser.print_help()
            exit(1)

        # goes to each command
        getattr(self, args.command)()

    def launcher(self):
        logger.info('Launching Web Interface')
        cmd = 'R -e shiny::runApp(\\"' + cfg.rootDir + '/MODApy-Shiny.R\\"\\,port=3838\\,host=\\"0.0.0.0\\")'
        print(cmd)
        webapp = subprocess.Popen(shlex.split(cmd), stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
        webapp.wait()
        output, error = webapp.communicate()
        logger.debug(output)
        logger.debug(error)
        logger.info('Web Interface Closed')

    def coverageStats(self):
        parser = argparse.ArgumentParser(
            description="Downloads Patient Data to Patients folder. Receives both url or xls/xlsx")
        parser.add_argument('Gene_Exon_Bed_File',
                            help='A Gene and Exon Bed file, in the format CHROM START END GENE_EXON STRAND')
        parser.add_argument("-Panel",
                            help="Filepath of a bed file containing a group of genes or exons of interest")
        parser.add_argument('Bam_files',
                            help='Bam files or list of files to calculate coverage stats (can use wildcard) Example: /home/bams/*.bam',
                            nargs='*')
        args = parser.parse_args(argv[2:])
        Bam_files = list(args.Bam_files)
        bed_file = args.Gene_Exon_Bed_File
        panel_file = args.Panel
        coverage.main(Bam_files, bed_file, panel_file)

    def addPatient(self):
        parser = argparse.ArgumentParser(
            description="Downloads Patient Data to Patients folder. Receives both url or xls/xlsx")
        parser.add_argument('FileorURL', help="URL to download or filepath to xls or xlsx files containing urls")
        args = parser.parse_args(argv[2:])
        fileorurl = args.FileorURL
        downloader.get_links(fileorurl)

    def variantsDB(self):
        parser = argparse.ArgumentParser(description="Work with Variants DB")
        parser.add_argument('-buildDB', action='store_true',
                            help='Build Variants DataBase with all patients that are not currently in it.')
        parser.add_argument('-addPatientToDB',
                            help='Adds a single patient to DB. Must supply path to vcf')
        parser.add_argument('-annotate',
                            help='Adds data from variantsdb to analysis done in modapy. Must supply path to excel output from modapy')
        args = parser.parse_args(argv[2:])
        if args.buildDB:
            db = variantsdb.VariantsDB.buildDB()
            db.to_VarDBCSV()
        if args.addPatientToDB:
            patient = cfg.patientPath + args.addPatientToDB
            db = variantsdb.VariantsDB.from_exceldb(variantsdb.variantsDBPath)
            db = db.addPatientToDB(patient)
            db.to_VarDBCSV()
        if args.annotate:
            fileName = args.annotate.rsplit('/', maxsplit=1)[1]
            patient = pd.read_excel(args.annotate)
            db = variantsdb.VariantsDB.from_csvdb(variantsdb.variantsDBPath.rsplit('/', maxsplit=1)[0])
            patient = db.annotate_excel(patient, fileName)

    def pipeline(self):
        # Description for pipeline usage
        parser = argparse.ArgumentParser(description="Run a Pipeline from FASTQ to VCF")
        parser.add_argument("-Pipeline", required=True, help="File name of the Pipeline inside Pipelines folder")
        parser.add_argument("-FQ", required=True,
                            help="Patient FastQ1 File Path - It needs to match exactly "
                                 "the filename found inside Patients folder. Only this one is needed for Single End."
                                 "Two FastQs will be needed for Paired End (usage: -FQ Fastq1 -FQ Fastq2",
                            action='append')
        parser.add_argument("-keeptmp", action="store_true", default=False,
                            help="Keep Temp files, otherwise just creates annotated vcf file.")
        parser.add_argument("-startStep", default=0, type=int, help="Defines step to start running pipeline.")
        parser.add_argument("-endStep", default=0, type=int, help="Defines step to start running pipeline.")

        # ignore first argument
        args = parser.parse_args(argv[2:])
        pipe = cfg.pipelinesPath + args.Pipeline

        checkFile(pipe, args.Pipeline.split('.')[-1])

        newpipe = pipeline.Pipeline.from_json(pipe)

        if len(args.FQ) > 2:
            logger.error('Only Two FASTQ files allowed. The Input for FastQ Files was: ', str(args.FQ))
            return exit(1)

        elif len(args.FQ) == 2:
            fq1 = cfg.patientPath + args.FQ[0]
            fq2 = cfg.patientPath + args.FQ[1]
            checkFile(fq1, '.' + fq1.split('.')[-1])
            checkFile(fq2, '.' + fq2.split('.')[-1])
            if args.keeptmp:
                newpipe.runpipeline(fq1, fq2, keeptmp=True, startStep=args.startStep, endStep=args.endStep)
            else:
                newpipe.runpipeline(fq1, fq2, startStep=args.startStep, endStep=args.endStep)
            return 0
        else:
            fq1 = cfg.patientPath + args.FQ[0]
            fq2 = ''
            checkFile(fq1, '.' + fq1.split('.')[-1])
            if args.keeptmp:
                newpipe.runpipeline(fq1, keeptmp=True, startStep=args.startStep, endStep=args.endStep)
            else:
                newpipe.runpipeline(fq1, startStep=args.startStep, endStep=args.endStep)
            return 0

    def single(self):
        # Description for panel usage
        parser = argparse.ArgumentParser(description="Run study on a single patient")
        parser.add_argument("-Panel", required=True, help="File name of Panel inside Panels folder")
        parser.add_argument("-Patient", required=True,
                            help="Patient File Path - It needs to match exactly to the one found inside Patients folder")
        # ignore first argument
        try:
            args = parser.parse_args(argv[2:])
            panel = cfg.panelsPath + args.Panel + '.xlsx'
            patient = cfg.patientPath + args.Patient
            ptCheck = checkFile(patient, '.vcf')
            pnCheck = checkFile(panel, '.xlsx')
            logger.info("Running %s on patient %s" % (str(args.Panel), str(args.Patient)))
            result = vcfmgr.ParsedVCF.from_vcf(patient).panel(panel)
            outpath = cfg.resultsPath + 'Panels/' + result.name + '/' + result.name + '_' + args.Panel + '.xlsx'
            os.makedirs(os.path.dirname(outpath), exist_ok=True)
            result.vcf_to_excel(outpath)
            logger.info('Single Analisis Complete')
            logger.info('File available at:%s' % outpath)
        except:
            logger.info('Single Analisis Failed')
        return 0

    def duos(self):
        # Description for duos usage
        parser = argparse.ArgumentParser(description="Run Duos Study on two patients")
        parser.add_argument("-Patient1", required=True,
                            help="Patient 1 File Path - It needs to match exactly to the one found inside Patients folder")
        parser.add_argument("-Patient2", required=True,
                            help="Patient 2 File Path - It needs to match exactly to the one found inside Patients folder")
        parser.add_argument("--VennPlace", const=None, choices=['A', 'B', 'A:B'], nargs='?',
                            help="Place in a Venn Diagram to obtain variants from")
        parser.add_argument("--Panel", nargs='?', const=None, help="Panel to run on Duos study")
        parser.add_argument("--Filter", nargs='?', const=None,
                            help="Filter to apply. This function will filter out every row that includes the given text"
                                 " in the given column. For filtering Empty data, TEXT keyword is 'Empty'",
                            metavar=("COLUMN TEXT"), action='append')
        # ignore first argument
        try:
            args = parser.parse_args(argv[2:])
            patient1 = cfg.patientPath + args.Patient1
            patient2 = cfg.patientPath + args.Patient2
            # Checks file existence and type for patients
            pt1Check = checkFile(patient1, '.vcf')
            pt2Check = checkFile(patient2, '.vcf')
            logger.info("Running Duos Study on %s and %s" % (str(args.Patient1), str(args.Patient2)))
            pvcfs = vcfmgr.ParsedVCF.mp_parser(patient1, patient2)
            result = pvcfs[0].duos(pvcfs[1], VENNPLACE=args.VennPlace)
            resultname = result.name
            outpath = cfg.resultsPath + 'Duos/' + result.name.replace(':', '_') + '/' + result.name.replace(':', '_')
            result.name = resultname
            if args.VennPlace is not None:
                outpath = outpath + '_Venn' + args.VennPlace.replace(':', '_')
            if args.Panel is not None:
                logger.info('Running panel {}'.format(args.Panel))
                panel = cfg.panelsPath + args.Panel + '.xlsx'
                checkFile(panel, '.xlsx')
                result = result.panel(panel)
                result.name = resultname
                outpath = outpath + '_P' + args.Panel
            if args.Filter[0] is not None:
                for x in args.Filter:
                    if (len(x.split())) != 2:
                        logger.error('--Filter accepts only two arguments. Usage: --Filter COLUMN_NAME TEXT_TO_FILTER')
                        exit(1)
                    else:
                        x = x.split()
                        if x[1] == 'Empty':
                            result = result[result[x[0]] != '']
                        else:
                            result = result[~result[x[0]].str.contains(x[1])]
                    result.name = resultname
                    outpath = outpath + '_F' + str(x[0]) + str(x[1])
            outpath = outpath + '.xlsx'
            logger.info('Writing Result File')
            result.vcf_to_excel(outpath)
            logger.info('Duos Analisis Complete')
            logger.info('File available at:%s' % outpath)
        except:
            logger.info('Duos Analisis Failed')
        return 0

    def trios(self):
        parser = argparse.ArgumentParser(description="Run Trios Study on two patients")
        parser.add_argument("-Patient1", required=True,
                            help="Patient 1 File Path - It needs to match exactly to the one found inside Patients folder")
        parser.add_argument("-Patient2", required=True,
                            help="Patient 2 File Path - It needs to match exactly to the one found inside Patients folder")
        parser.add_argument("-Patient3", required=True,
                            help="Patient 3 File Path - It needs to match exactly to the one found inside Patients folder")
        parser.add_argument("--Panel", nargs='?', const=None, help="Panel to run on Trios study")
        parser.add_argument("--Filter", nargs='?', const=None,
                            help="Filter to apply. This function will filter out every row that includes the given text"
                                 " in the given column. For filtering Empty data, TEXT keyword is 'Empty'",
                            metavar=("COLUMN TEXT"), action='append')
        parser.add_argument("--VennPlace", default=None, const=None, nargs='?',
                            choices=['A', 'B', 'C', 'A:B', 'A:C', 'B:C', 'A:B:C', 'ALL'],
                            help="Place in a Venn Diagram to obtain variants from")
        try:
            # ignore first argument
            args = parser.parse_args(argv[2:])
            patient1 = cfg.patientPath + args.Patient1
            patient2 = cfg.patientPath + args.Patient2
            patient3 = cfg.patientPath + args.Patient3
            # Checks file existence and type for patients
            pt1Check = checkFile(patient1, '.vcf')
            pt2Check = checkFile(patient2, '.vcf')
            pt3Check = checkFile(patient3, '.vcf')
            logger.info(
                "Running Trios Study on %s, %s and %s" % (str(args.Patient1), str(args.Patient2), str(args.Patient3)))
            pvcfs = vcfmgr.ParsedVCF.mp_parser(patient1, patient2, patient3)
            result = pvcfs[0].duos(pvcfs[1]).duos(pvcfs[2], VENNPLACE=args.VennPlace)
            resultname = result.name
            outpath = cfg.resultsPath + 'Trios/' + result.name.replace(':', '_') + '/' + result.name.replace(':', '_')
            result.name = resultname
            if args.VennPlace is not None:
                outpath = outpath + '_Venn' + args.VennPlace.replace(':', '_')
            # check if there is a Panel Requested
            if args.Panel:
                logger.info('Running panel {}'.format(args.Panel))
                panel = cfg.panelsPath + args.Panel + '.xlsx'
                checkFile(panel, '.xlsx')
                result = result.panel(panel)
                result.name = resultname
                outpath = outpath + '_Panel' + args.Panel
            # check if there is a Filter Requested
            if args.Filter[0] is not None:
                for x in args.Filter:
                    if (len(x.split())) != 2:
                        logger.error('--Filter accepts only two arguments. Usage: --Filter COLUMN_NAME TEXT_TO_FILTER')
                        exit(1)
                    else:
                        x = x.split()
                        if x[1] == 'Empty':
                            result = result[result[x[0]] != '']
                        else:
                            result = result[~result[x[0]].str.contains(x[1])]
                        result.name = resultname
                        outpath = outpath + '_Filter' + str(x[0]) + str(x[1])
            outpath = outpath + '.xlsx'
            result.vcf_to_excel(outpath)
            logger.info('Trios Analisis Complete')
            logger.info('File available at:%s' % outpath)
        except:
            logger.info('Trios Analisis Failed')
        return 0


def checkFile(filePath, extension):
    if os.path.isfile(filePath):
        fileName, fileExtension = os.path.splitext(filePath)
        if extension == fileExtension:
            return True
    else:
        logger.error("%s couldn't be found. Please check if file exists and that it's extension is %s" % (filePath,
                                                                                                          extension))
        exit(1)


def main():
    cfg.setup_logging()
    Parser()


if __name__ == "__main__":
    main()
