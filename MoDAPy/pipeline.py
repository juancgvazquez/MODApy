import json, datetime, shlex, logging
from subprocess import Popen, STDOUT, PIPE, CalledProcessError
from os import path

'''
Class defined for each step of the Pipeline
'''


class PipeStep(object):
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

	def stepInfo(self):
		print('Name:', self.name)
		print('Command', self.command + self.version, self.subcommand, self.args)
		print('Input File:', self.inputfile)
		print('Output File:', self.outputfile)


'''
General Class for Pipelines
'''


class Pipeline(object):
	def __init__(self, name, referencepath, url='', description=''):
		self.name = name
		self.url = url
		self.description = description
		self.referencepath = referencepath
		self.required_files = []
		self.steps = []

	def add_steps(self, step):
		self.steps.append(step)

	def add_req_files(self, path):
		self.required_files.append(path)

	'''
	Class Method to create Pipeline from Json
	'''

	@classmethod
	def from_json(cls, jsonpath):
		with open(jsonpath) as f:
			jsf = json.load(f)

		name = jsf['INFO']['name']
		url = jsf['INFO']['url']
		description = jsf['INFO']['description']
		reference = jsf['INFO']['reference']

		newpipe = Pipeline(name, reference, url, description)

		for x in jsf['INFO']['required_files']:
			newpipe.add_req_files(x)

		for x in jsf['STEPS']:
			name = jsf['STEPS'][x]['name']
			command = jsf['STEPS'][x]['command']
			subcommand = jsf['STEPS'][x]['subcommand']
			version = jsf['STEPS'][x]['version']
			inputfile = jsf['STEPS'][x]['input']
			outputfile = jsf['STEPS'][x]['output']
			args = jsf['STEPS'][x]['args']
			newstep = PipeStep(name, command, subcommand, version, inputfile, outputfile, args)
			newpipe.add_steps(newstep)

		return newpipe

	'''
	Method to run selected Pipeline on fastq files
	'''

	def runPipeline(self, fastq1: str, fastq2=None):
		patientname = fastq1.split('/')[-1].split('.')[0].split('_')[0]
		ref = self.referencepath
		print('Running', self.name, 'pipeline on patient:', patientname)
		# bool to check if first step
		first = True
		for step in self.steps:
			# Checks if first step, we should input the exact input as inputfiles
			if first == True:
				first = False
				if type(step.inputfile) == list:
					if (fastq1 != None) & (fastq2 != None):
						inputfile = fastq1 + fastq2
					else:
						print('WARNING: This pipeline was designed for Pair End and you are running it as Single End')
						inputfile = fastq1
				elif type(step.inputfile) == str:
					if (fastq1 != None) & (fastq2 != None):
						print('WARNING: This pipeline was designed for Single End and you are running it as Pair End')
						inputfile = fastq1 + fastq2
					else:
						inputfile = fastq1
				else:
					return 'Error Parsing input file. It should be a string or list of strings.'
			# If it's not first step, input depends on output of previous step + patientname
			else:
				inputfile = step.inputfile.replace('patientname', patientname)

			# replaces patient name in outputfiles
			if type(step.outputfile) == str:
				outputfile = step.outputfile.replace('patientname', patientname)
			else:
				return 'Error Parsing output file. It should be a string.'

			print(step.name)
			rootdir = path.dirname(path.abspath(__file__))
			args = step.args.replace('patientname', patientname).replace('reference', ref)
			cmdver = step.version.replace('.', '_')
			cmd = rootdir + '/bin/' + step.command + '/' + step.command + '_' + cmdver + ' ' + step.subcommand
			cmdstr = cmd + ' ' + args + ' ' + ' ' + inputfile + ' ' + outputfile
			cmd = shlex.split(cmdstr)
			formatter2 = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
			logging.basicConfig(level=logging.DEBUG, filename='log', format=formatter2)
			console = logging.StreamHandler()
			console.setLevel(logging.INFO)
			formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
			console.setFormatter(formatter)
			logging.getLogger('').addHandler(console)
			logging.info('Subprocess: "' + cmdstr)
			try:
				cmdrun = Popen(cmd, stderr=PIPE, stdout=PIPE, universal_newlines=True)
			except (OSError, CalledProcessError) as exception:
				out, err = cmdrun.communicate()
				if out != '':
					logging.debug('Output: ' + out.strip())
				if err != '':
					logging.debug('Stderr: ' + err.strip())
				logging.debug('Subprocess failed')
				logging.debug('Exception ocurred: ' + str(exception))
				logging.info('There was an error when running the pipeline. Please check logs for more info')
				exit(1)
			else:
				logging.info('Subprocess finished')

	'''
	Method to print Pipeline Info
	'''

	def pipelineInfo(self):
		print('Name:', self.name)
		print('Reference:', self.referencepath)
		print('URL:', self.url)
		print('Description:', self.description)
		print('Additional Files Required:', self.required_files)
		print('Steps:', self.steps)
