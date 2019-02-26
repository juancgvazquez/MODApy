#!/usr/bin/env python
import argparse
import logging
import os
import shlex
import subprocess
from sys import argv

from MODApy import filemgr, cfg, panelmdl, pipeline, vcfmgr
from MODApy.version import __version__

logger = logging.getLogger(__name__)


class Parser(object):
    '''
    Main parser for command line arguments. Uses ArgParse
    '''

    def __init__(self):
        parser = argparse.ArgumentParser(description="Multi-Omics Data Analisis for Python", usage='''MODApy <command> [<args>]

        Commands:
        launcher Run MoDAPy Web Interface
        pipeline Run pipeline on FastQ file/s
        single   Run study on a single patient
        duos     Run Duos analysis on two selected patients
        trios    Run Trios analysis on three selected patients

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
        # TODO:Set shiny app inside MODAPY
        logger.info('Launching Web Interface')
        cmd = 'R --vanilla -e shiny::runApp(\\"' + cfg.rootDir + '/MODApy-Shiny.R\\"\\,port=8081\\,launch\\.browser=TRUE)'
        webapp = subprocess.Popen(shlex.split(cmd), stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        webapp.wait()
        logger.info('Web Interface Closed')

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

        # ignore first argument
        args = parser.parse_args(argv[2:])
        pipe = cfg.pipelinesPath + args.Pipeline

        filemgr.checkFile(pipe, args.Pipeline.split('.')[-1])

        newpipe = pipeline.Pipeline.from_json(pipe)

        if len(args.FQ) > 2:
            logger.error('Only Two FASTQ files allowed. The Input for FastQ Files was: ', str(args.FQ))
            return exit(1)

        elif len(args.FQ) == 2:
            fq1 = cfg.patientPath + args.FQ[0]
            fq2 = cfg.patientPath + args.FQ[1]
            filemgr.checkFile(fq1, '.' + fq1.split('.')[-1])
            filemgr.checkFile(fq2, '.' + fq2.split('.')[-1])
            if args.keeptmp:
                newpipe.runpipeline(fq1, fq2, keeptmp=True, startStep=args.startStep)
            else:
                newpipe.runpipeline(fq1, fq2, startStep=args.startStep)
            return 0
        else:
            fq1 = cfg.patientPath + args.FQ[0]
            fq2 = ''
            filemgr.checkFile(fq1, '.' + fq1.split('.')[-1])
            if args.keeptmp:
                newpipe.runpipeline(fq1, keeptmp=True, startStep=args.startStep)
            else:
                newpipe.runpipeline(fq1, startStep=args.startStep)
            return 0

    def single(self):
        # Description for panel usage
        parser = argparse.ArgumentParser(description="Run study on a single patient")
        parser.add_argument("-Panel", required=True, help="File name of Panel inside Panels folder")
        parser.add_argument("-Patient", required=True,
                            help="Patient File Path - It needs to match exactly to the one found inside Patients folder")
        # ignore first argument
        args = parser.parse_args(argv[2:])
        panel = cfg.panelsPath + args.Panel + '.xlsx'
        patient = cfg.patientPath + args.Patient
        ptCheck = filemgr.checkFile(patient, '.vcf')
        pnCheck = filemgr.checkFile(panel, '.xlsx')

        logger.info("Running %s on patient %s" % (str(args.Panel), str(args.Patient)))
        result = panelmdl.panelrun(panel, patient)
        outpath = cfg.resultsPath + 'Panels/' + result.name + '/' + result.name + '_' + args.Panel + '.xlsx'
        os.makedirs(os.path.dirname(outpath), exist_ok=True)
        filemgr.df_to_excel(result, outpath)
        logger.info('Single Analisis Complete')
        return 0

    def duos(self):
        # Description for duos usage
        parser = argparse.ArgumentParser(description="Run Duos Study on two patients")
        parser.add_argument("-Patient1", required=True,
                            help="Patient 1 File Path - It needs to match exactly to the one found inside Patients folder")
        parser.add_argument("-Patient2", required=True,
                            help="Patient 2 File Path - It needs to match exactly to the one found inside Patients folder")
        parser.add_argument("--VennPlace", default='ALL', const='ALL', choices=['A', 'B', 'A:B', 'ALL'], nargs='?',
                            help="Place in a Venn Diagram to obtain variants from")
        parser.add_argument("--Panel", nargs='?', const=None, help="Panel to run on Duos study")
        parser.add_argument("--Filter", nargs='?', const=None,
                            help="Filter to apply. This function will filter out every row that includes the given text"
                                 " in the given column. For filtering Empty data, TEXT keyword is 'Empty'",
                            metavar=("COLUMN TEXT"), action='append')
        # ignore first argument
        args = parser.parse_args(argv[2:])
        patient1 = cfg.patientPath + args.Patient1
        patient2 = cfg.patientPath + args.Patient2
        # Checks file existence and type for patients
        pt1Check = filemgr.checkFile(patient1, '.vcf')
        pt2Check = filemgr.checkFile(patient2, '.vcf')
        logger.info("Running Duos Study on %s and %s" % (str(args.Patient1), str(args.Patient2)))
        result = vcfmgr.ParsedVCF.from_vcf(patient1).duos(patient2)
        resultname = result.name
        outpath = cfg.resultsPath + 'Duos/' + result.name.replace(':', '_')
        filemgr.getstats(result, 1)
        # check if there is a special Venn Place Requested
        if args.VennPlace == 'A':
            if len(result[result['DUOS'] == result.name.split(':')[0]]) > 0:
                result = result[result['DUOS'] == result.name.split(':')[0]]
                outpath += '_Venn' + resultname.split(':')[0]
            else:
                logger.info('No variants present in selected Venn Place!')
                exit(1)
        if args.VennPlace == 'B':
            if len(result[result['DUOS'] == result.name.split(':')[1]]) > 0:
                result = result[result['DUOS'] == result.name.split(':')[1]]
                outpath += '_Venn' + resultname.split(':')[1]
            else:
                logger.info('No variants present in selected Venn Place!')
                logger.error('Exiting due to no variants present in selected Venn Place!')
                exit(1)
        if args.VennPlace == 'A:B':
            if len(result[result['DUOS'] == ':'.join([result.name.split(':')[0], result.name.split(':')[1]])]) > 0:
                result = result[result['DUOS'] == ':'.join([result.name.split(':')[0], result.name.split(':')[1]])]
                outpath += '_Venn' + resultname.replace(':', '_')
            else:
                logger.info('No variants present in selected Venn Place!')
                logger.error('Exiting due to no variants present in selected Venn Place!')
                exit(1)
        result.name = resultname
        if args.Panel is not None:
            panel = cfg.panelsPath + args.Panel + '.xlsx'
            filemgr.checkFile(panel, '.xlsx')
            result = panelmdl.panelrun(panel, result)
            if len(result == 0):
                logger.error('No variants present after running selected Panel on selected Data')
                exit(1)
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
        filemgr.df_to_excel(result, outpath)
        logger.info('Duos Analisis Complete')
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
        parser.add_argument("--VennPlace", default='ALL', const='ALL', nargs='?',
                            choices=['A', 'B', 'C', 'A:B', 'A:C', 'B:C', 'A:B:C', 'ALL'],
                            help="Place in a Venn Diagram to obtain variants from")
        # ignore first argument
        args = parser.parse_args(argv[2:])
        patient1 = cfg.patientPath + args.Patient1
        patient2 = cfg.patientPath + args.Patient2
        patient3 = cfg.patientPath + args.Patient3
        # Checks file existence and type for patients
        pt1Check = filemgr.checkFile(patient1, '.vcf')
        pt2Check = filemgr.checkFile(patient2, '.vcf')
        pt3Check = filemgr.checkFile(patient3, '.vcf')
        logger.info(
            "Running Trios Study on %s, %s and %s" % (str(args.Patient1), str(args.Patient2), str(args.Patient3)))
        result = vcfmgr.ParsedVCF.from_vcf(patient1).trios(patient2, patient3)
        resultname = result.name
        outpath = cfg.resultsPath + 'Trios/' + result.name.replace(':', '_')
        filemgr.getstats(result, 1)
        # check if there is a special Venn Place Requested
        if args.VennPlace == 'A':
            if (len(result[result['TRIOS'] == result.name.split(':')[0]])) > 0:
                result = result[result['TRIOS'] == result.name.split(':')[0]]
                outpath += '_Venn' + resultname.split(':')[0]
            else:
                logger.info('No variants present in selected Venn Place!')
                logger.error('Exiting due to no variants present in selected Venn Place!')
                exit(1)
        if args.VennPlace == 'B':
            if (len(result[result['TRIOS'] == result.name.split(':')[1]])) > 0:
                result = result[result['TRIOS'] == result.name.split(':')[1]]
                outpath += '_Venn' + resultname.split(':')[1]
            else:
                logger.info('No variants present in selected Venn Place!')
                logger.error('Exiting due to no variants present in selected Venn Place!')
                exit(1)
        if args.VennPlace == 'C':
            if (len(result[result['TRIOS'] == result.name.split(':')[2]])) > 0:
                result = result[result['TRIOS'] == result.name.split(':')[2]]
                outpath += '_Venn' + resultname.split(':')[2]
            else:
                logger.info('No variants present in selected Venn Place!')
                logger.error('Exiting due to no variants present in selected Venn Place!')
                exit(1)
        if args.VennPlace == 'A:B':
            if (len(result[result['TRIOS'] == ':'.join([result.name.split(':')[0], result.name.split(':')[1]])])) > 0:
                result = result[result['TRIOS'] == ':'.join([result.name.split(':')[0], result.name.split(':')[1]])]
                outpath += '_Venn' + '_'.join([resultname.split(':')[0], resultname.split(':')[1]])
            else:
                logger.info('No variants present in selected Venn Place!')
                logger.error('Exiting due to no variants present in selected Venn Place!')
                exit(1)
        if args.VennPlace == 'A:C':
            if (len(result[result['TRIOS'] == ':'.join([result.name.split(':')[0], result.name.split(':')[2]])])) > 0:
                result = result[result['TRIOS'] == ':'.join([result.name.split(':')[0], result.name.split(':')[2]])]
                outpath += '_Venn' + '_'.join([resultname.split(':')[0], resultname.split(':')[2]])
            else:
                logger.info('No variants present in selected Venn Place!')
                logger.error('Exiting due to no variants present in selected Venn Place!')
                exit(1)
        if args.VennPlace == 'B:C':
            if (len(result[result['TRIOS'] == ':'.join([result.name.split(':')[1], result.name.split(':')[2]])])) > 0:
                result = result[result['TRIOS'] == ':'.join([result.name.split(':')[1], result.name.split(':')[2]])]
                outpath += '_Venn' + '_'.join([resultname.split(':')[1], resultname.split(':')[2]])
            else:
                logger.info('No variants present in selected Venn Place!')
                logger.error('Exiting due to no variants present in selected Venn Place!')
                exit(1)
        if args.VennPlace == 'A:B:C':
            if (len(result[result['TRIOS'] == ':'.join(
                    [result.name.split(':')[0], result.name.split(':')[1], result.name.split(':')[2]])])) > 0:
                result = result[result['TRIOS'] == ':'.join([result.name.split(':')[0], result.name.split(':')[1]])]
                outpath += '_Venn' + '_'.join(
                    [resultname.split(':')[0], resultname.split(':')[1], resultname.split(':')[2]])
            else:
                logger.info('No variants present in selected Venn Place!')
                logger.error('Exiting due to no variants present in selected Venn Place!')
                exit(1)
        result.name = resultname
        # check if there is a Panel Requested
        if args.Panel:
            panel = cfg.panelsPath + args.Panel + '.xlsx'
            result = panelmdl.panelrun(panel, result)
            if len(result == 0):
                logger.error('No variants present after running selected Panel on selected Data')
                exit(1)
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
        filemgr.df_to_excel(result, outpath)
        logger.info('Trios Analisis Complete')
        return 0


def main():
    cfg.setup_logging()
    Parser()


if __name__ == "__main__":
    main()
