import json
import logging
import os
import shlex
import shutil
import subprocess

from MODApy import cfg


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
    def _buildpipe(pipefile):
        """
         Private Class Method to build pipeline from loaded json,xml or yaml
         Parameters
         --------------------
         pipefile File loaded by from_json, from_xml or from_yaml
        """
        name = pipefile['INFO']['name']
        url = pipefile['INFO']['url']
        description = pipefile['INFO']['description']
        reference = pipefile['INFO']['reference']

        newpipe = Pipeline(name, reference, url, description)

        for x in pipefile['INFO']['required_files']:
            newpipe.add_req_files(x)

        for x in range(1, len(pipefile['STEPS']) + 1):
            name = pipefile['STEPS'][str(x)]['name']
            command = pipefile['STEPS'][str(x)]['command']
            subcommand = pipefile['STEPS'][str(x)]['subcommand']
            version = pipefile['STEPS'][str(x)]['version']
            inputfile = pipefile['STEPS'][str(x)]['input']
            outputfile = pipefile['STEPS'][str(x)]['output']
            args = pipefile['STEPS'][str(x)]['args']
            newstep = PipeStep(name, command, subcommand, version, inputfile, outputfile, args)
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
            pipefile = json.load(f)

        builtpipe = cls._buildpipe(pipefile)

        return builtpipe

    def runpipeline(self, fastq1: str, fastq2=None, keeptmp=False, startStep=0):
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
        patientname = fastq1.split('/')[-1].split('.')[0].split('_')[0]
        ref = cfg.referencesPath + self.reference + '/' + self.reference + '.fa'
        pipedir = "".join(x for x in self.name if x.isalnum())
        tmpdir = cfg.resultsPath + 'Pipelines/' + patientname + '/' + pipedir + '/tmp/'
        os.makedirs(tmpdir, exist_ok=True)
        print('Running', self.name, 'pipeline on patient:', patientname)
        # bool to check if first step
        if startStep == 0:
            first = True
        else:
            first = False
        for step in self.steps[startStep:]:
            # Checks if first step, we should input the exact input as inputfiles
            if first is True:
                first = False
                if type(step.inputfile) == list:
                    if (fastq1 is not None) & (fastq2 is not None):
                        inputfile = fastq1 + ' ' + fastq2
                    else:
                        print('WARNING: This pipeline was designed for Pair End and you are running it as Single End')
                        inputfile = fastq1
                elif type(step.inputfile) == str:
                    if (fastq1 is not None) & (fastq2 is not None):
                        print('WARNING: This pipeline was designed for Single End and you are running it as Pair End')
                        inputfile = fastq1 + ' ' + fastq2
                    else:
                        inputfile = fastq1
                else:
                    return 'Error Parsing input file. It should be a string or list of strings.'
            # If it's not first step, input depends on output of previous step + patientname
            else:
                inputfile = step.inputfile.replace('patientname', tmpdir + patientname)

            # replaces patient name in outputfiles
            if type(step.outputfile) == str:
                outputfile = step.outputfile.replace('patientname', tmpdir + patientname)
            else:
                return 'Error Parsing output file. It should be a string.'

            print(step.name)
            args = step.args.replace('patientname', tmpdir + patientname).replace('reference', ref)
            cmdver = step.version.replace('.', '_')
            javacmds = ['GATK', 'picard', 'SnpSift', 'snpEff']
            if any(javacmd in step.command for javacmd in javacmds):
                cmd = 'java -jar ' + cfg.binPath + step.command + '/' + step.command + '_' + cmdver \
                      + '.jar ' + step.subcommand
            else:
                cmd = cfg.binPath + step.command + '/' + step.command + '_' + cmdver + ' ' + step.subcommand

            cmdstr = cmd + ' ' + args + ' ' + ' ' + inputfile + ' ' + outputfile

            cmd = shlex.split(cmdstr)

            # logging stuff
            formatter2 = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            logging.basicConfig(level=logging.DEBUG, filename='log', format=formatter2)
            console = logging.StreamHandler()
            console.setLevel(logging.INFO)
            formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
            console.setFormatter(formatter)
            logging.getLogger('').addHandler(console)
            # done logging config
            logging.info('Subprocess: ' + cmdstr)
            stdcmds = ['bwa', 'snpEff', 'SnpSift']
            try:
                if any(stdcmd in s for s in cmd for stdcmd in stdcmds):
                    output = cmd[-1]
                    del cmd[-1]
                    with open(output, 'w+') as out:
                        cmdrun = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=out, universal_newlines=True)
                else:
                    cmdrun = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE,
                                              universal_newlines=True)
                out, err = cmdrun.communicate()
                if out is not None:
                    logging.debug('Output: ' + out.strip())
                if err is not None:
                    logging.debug('Stderr: ' + err.strip())
                if cmdrun.returncode != 0:
                    logging.error('Subprocess failed with error code: ' + str(cmdrun.returncode))
                    logging.error('Check log for more details')
                    exit(cmdrun.returncode)

            except (OSError, subprocess.CalledProcessError) as exception:
                logging.debug('Subprocess failed')
                logging.debug('Exception ocurred: ' + str(exception))
                logging.info('There was an error when running the pipeline. Please check logs for more info')
                exit(1)
            else:
                logging.info('Subprocess finished')
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
