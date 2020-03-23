import json
import logging
import os
import shlex
import shutil
import subprocess

import xmltodict
import yaml

from MODApy import cfg

logger = logging.getLogger(__name__)
logger2 = logging.getLogger('Pipeline Module')
hdlr = logging.FileHandler(cfg.rootDir + '/logs/pipe_run.log')
formatter = logging.Formatter("%(asctime)s %(name)-25s %(levelname)-8s %(message)s")
hdlr.setFormatter(formatter)
logger2.addHandler(hdlr)
logger2.setLevel(logging.DEBUG)


# TODO: restructure pipeline, being less open, so less prone to errors


class PipeStep(object):
    '''
    Class defined for each step of the Pipeline
    '''

    def __init__(self, name, command, subcommand, version, inputfile, outputfile, args):
        self.name = name
        self.command = command
        self.subcommand = subcommand
        self.version = version
        self.inputfile = inputfile
        self.outputfile = outputfile
        self.args = args

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    def stepinfo(self):
        print('Name:', self.name)
        print('Command', self.command + self.version, self.subcommand, self.args)
        print('Input File:', self.inputfile)
        print('Output File:', self.outputfile)


class Pipeline(object):
    '''
    General Class for Pipelines
    '''

    def __init__(self, name, reference, url='', description=''):
        self.name = name
        self.url = url
        self.description = description
        self.reference = reference
        self.required_files = []
        self.steps = []

    def add_steps(self, step):
        self.steps.append(step)

    def add_req_files(self, path):
        self.required_files.append(path)

    @staticmethod
    def _buildpipe(pipedict):
        """
         Private Class Method to build pipeline from loaded json,xml or yaml
         Parameters
         --------------------
         pipefile File loaded by from_json, from_xml or from_yaml
        """
        name = pipedict['INFO']['name']
        url = pipedict['INFO']['url']
        description = pipedict['INFO']['description']
        reference = pipedict['INFO']['reference']

        newpipe = Pipeline(name, reference, url, description)

        for x in pipedict['INFO']['required_files']:
            newpipe.add_req_files(x)

        steps = list(pipedict['STEPS'].values())

        for i in range(0, len(steps)):
            name = steps[i]['name']
            command = steps[i]['command']
            subcommand = steps[i]['subcommand']
            version = steps[i]['version']
            inputfile = steps[i]['input']
            outputfile = steps[i]['output']
            args = steps[i]['args']
            newstep = PipeStep(name, command, subcommand,
                               version, inputfile, outputfile, args)
            newpipe.add_steps(newstep)

        return newpipe

    @classmethod
    def from_json(cls, jsonpath):
        """
        Class Method to read json and build the Pipeline
        Parameters
        ----------
        jsonpath
            Path to the json file used to create the Pipeline.

        Returns Pipeline object
        """
        with open(jsonpath) as f:
            pipedict = json.load(f)

        builtpipe = cls._buildpipe(pipedict)

        return builtpipe

    @classmethod
    def from_yaml(cls, yamlpath):
        """
        Class Method to read yaml and build the Pipeline
        Parameters
        ----------
        jsonpath
            Path to the yaml file used to create the Pipeline.

        Returns Pipeline object
        """
        with open(yamlpath) as f:
            pipedict = yaml.load(f)
        builtpipe = cls._buildpipe(pipedict)

        return builtpipe

    @classmethod
    def from_xml(cls, xmlpath):
        """
        Class Method to parse xml and build the Pipeline
        Parameters
        ----------
        xmlpath
            Path to the xml file used to create the Pipeline.

        Returns Pipeline object
        """
        with open(xmlpath) as f:
            pipedict = xmltodict.parse(f.read())['root']
        builtpipe = cls._buildpipe(pipedict)

        return builtpipe

    def runpipeline(self, fastq1: str, fastq2=None, keeptmp=False, startStep=0, endStep=0):
        """
        Method to run the Pipeline
        Parameters
        ----------
        fastq1
            Path to the first fastq file.
        fastq2
            Path to the second fastq file, in case of paired reads.
        """

        '''
        Method to run selected Pipeline on fastq files
        '''
        logger.info(self.steps)
        logger.info("Nro de Pasos: %s" % str(len(self.steps)))
        logger2.info(self.steps)
        logger2.info("Nro de Pasos: %s" % str(len(self.steps)))
        patientname = fastq1.split('/')[-1].split('.')[0].split('_')[0]
        ref = cfg.referencesPath + self.reference + '/' + self.reference + '.fa'
        pipedir = "".join(x for x in self.name if x.isalnum())
        if cfg.testFlag:
            tmpdir = cfg.testPath + 'Pipelines/' + patientname + '/' + pipedir + '/tmp/'
        else:
            tmpdir = cfg.resultsPath + 'Pipelines/' + patientname + '/' + pipedir + '/tmp/'

        os.makedirs(tmpdir, exist_ok=True)
        samplename = patientname + '_MODApy'
        logger2.info('Running ' + str(self.name) + ' pipeline on patient: ' + str(patientname))
        # bool to check if first step
        first = True
        if endStep == 0:
            endStep = len(self.steps) + 1
        for step in self.steps[startStep:endStep]:
            if first == True:
                logger2.debug('First Step')
                first = False
                if type(step.inputfile) == list:
                    if any(x in y for y in [fastq1, fastq2] for x in ['.fastq', '.fastq.gz', '.fq', '.fq.gz']):
                        if (fastq1 is not None) & (fastq2 is not None):
                            inputfile = fastq1 + ' ' + fastq2
                        else:
                            print(
                                'WARNING: This pipeline was designed for Pair End and you are running it as Single End')
                            inputfile = fastq1
                elif type(step.inputfile) == str:
                    if any(x in y for y in [fastq1] for x in ['.fastq', '.fastq.gz', '.fq', '.fq.gz']):
                        if (fastq1 is not None) & (fastq2 is not None):
                            print(
                                'WARNING: This pipeline was designed for Single End and you are running it as Pair End')
                            inputfile = fastq1 + ' ' + fastq2
                    else:
                        inputfile = fastq1
            # If it's not first step, input depends on output of previous step + patientname
            else:
                if type(step.inputfile) == list:
                    inputfile1 = step.inputfile[0].replace(
                        'patientname', tmpdir + patientname + '/' + patientname)
                    inputfile2 = step.inputfile[1].replace(
                        'patientname', tmpdir + patientname + '/' + patientname)
                    inputfile = inputfile1 + ' ' + inputfile2
                elif type(step.inputfile) == str:
                    inputfile = step.inputfile.replace(
                        'patientname', tmpdir + patientname)
            # replaces patient name in outputfiles
            if type(step.outputfile) == str:
                outputfile = step.outputfile.replace(
                    'patientname', tmpdir + patientname)
            else:
                return 'Error Parsing output file. It should be a string.'

            logger2.info(step.name)
            args = step.args.replace(
                'patientname', tmpdir + patientname).replace('reference', ref).replace('samplename', patientname)
            cmdver = step.version.replace('.', '_')
            javacmds = ['GATK', 'picard', 'SnpSift', 'snpEff']
            if any(javacmd in step.command for javacmd in javacmds):
                cmd = 'java -jar -Xmx12G -Djava.io.tmpdir=%s ' % tmpdir + cfg.binPath + step.command + '/' + step.command + '_' + cmdver \
                      + '.jar ' + step.subcommand
            else:
                cmd = cfg.binPath + step.command + '/' + \
                      step.command + '_' + cmdver + ' ' + step.subcommand
            if 'HaplotypeCaller' in cmd:
                cmdstr = cmd + ' ' + args + ' ' + ' -I ' + inputfile + ' ' + outputfile
            else:
                cmdstr = cmd + ' ' + args + ' ' + ' ' + inputfile + ' ' + outputfile
                print(cmd)
            cmd = shlex.split(cmdstr)

            logging.info('Subprocess: ' + cmdstr)
            logger2.info('Subprocess: ' + cmdstr)
            stdcmds = ['bwa', 'bedtools', 'snpEff', 'SnpSift']
            print(cmd)
            try:
                if any(stdcmd in s for s in cmd for stdcmd in stdcmds):
                    output = cmd[-1]
                    del cmd[-1]
                    with open(output, 'w+') as out:
                        cmdrun = subprocess.Popen(
                            cmd, stderr=subprocess.PIPE, stdout=out, universal_newlines=True)
                else:
                    cmdrun = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE,
                                              universal_newlines=True)
                out, err = cmdrun.communicate()
                if out is not None:
                    logging.debug('Output: ' + out.strip())
                    logger2.debug('Output: ' + out.strip())
                if err is not None:
                    logging.debug('Stderr: ' + err.strip())
                    logger2.debug('Stderr: ' + err.strip())
                if cmdrun.returncode != 0:
                    logging.error(
                        'Subprocess failed with error code: ' + str(cmdrun.returncode))
                    logger2.error(
                        'Subprocess failed with error code: ' + str(cmdrun.returncode))
                    logging.error('Check log for more details')
                    logger2.error('Check log for more details')
                    exit(cmdrun.returncode)

            except (OSError, subprocess.CalledProcessError) as exception:
                logging.debug('Subprocess failed')
                logging.debug('Exception ocurred: ' + str(exception))
                logger2.debug('Subprocess failed')
                logger2.debug('Exception ocurred: ' + str(exception))
                logging.info(
                    'There was an error when running the pipeline. Please check logs for more info')
                logger2.info(
                    'There was an error when running the pipeline. Please check logs for more info')
                exit(1)
            else:
                logging.info('Subprocess finished')
                logger2.info('Subprocess finished')
        if cfg.testFlag:
            if os.path.exists(tmpdir + patientname + "_MODApy.final.vcf"):
                shutil.move(tmpdir + patientname + "_MODApy.final.vcf",
                            cfg.testPath + patientname+'_MODApy/' + patientname + "_MODApy.final.vcf")
            if os.path.exists(tmpdir + patientname + "_realigned_reads_recal.bam"):
                shutil.move(tmpdir + patientname + "_realigned_reads_recal.bam",
                            cfg.testPath + patientname+'_MODApy/' + patientname + "MODApy_realigned_reads_recal.bam")
        else:
            if os.path.exists(tmpdir + patientname + "_MODApy.final.vcf"):
                shutil.move(tmpdir + patientname + "_MODApy.final.vcf",
                            cfg.patientPath + patientname+'_MODApy/' + patientname + "_MODApy.final.vcf")
            if os.path.exists(tmpdir + patientname + "_realigned_reads_recal.bam"):
                shutil.move(tmpdir + patientname + "_realigned_reads_recal.bam",
                            cfg.patientPath + patientname+'_MODApy/' + patientname + "MODApy_realigned_reads_recal.bam")
        if keeptmp is False:
            shutil.rmtree(tmpdir)

    def pipelineinfo(self):
        '''
        Method to print Pipeline Info
        '''
        print('Name:', self.name)
        print('Reference:', self.reference)
        print('URL:', self.url)
        print('Description:', self.description)
        print('Additional Files Required:', self.required_files)
        print('Steps:', self.steps)
