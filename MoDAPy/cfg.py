import configparser
import os

cfg = configparser.ConfigParser()
cfgpath = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config.ini')
cfg.read(cfgpath)

patientPath = cfg['PATHS']['PatientPath']
panelsPath = cfg['PATHS']['PanelsPath']
reportsPath = cfg['PATHS']['ReportsPath']
resultsPath = cfg['PATHS']['ResultsPath']

def setcfg(cfgdict:dict):
	#for x in cfgdict:
	#	x = cfg[x]
	pass