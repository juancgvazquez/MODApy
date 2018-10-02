import configparser
import os

cfg = configparser.ConfigParser()
cfgpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.ini0)
cfg.read(cfgpath)

patientPath = cfg['PATHS']['PatientPath']
panelsPath = cfg['PATHS']['PanelsPath']
reportsPath = cfg['PATHS']['ReportsPath']
resultsPath = cfg['PATHS']['ResultsPath']
