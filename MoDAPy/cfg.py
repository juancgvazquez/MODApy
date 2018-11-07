import configparser, os

cfg = configparser.ConfigParser()
cfgPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config.ini')
cfg.read(cfgPath)

rootDir = os.path.dirname(os.path.abspath(__file__))

patientPath = cfg['PATHS']['PatientPath']
panelsPath = cfg['PATHS']['PanelsPath']
reportsPath = cfg['PATHS']['ReportsPath']
resultsPath = cfg['PATHS']['ResultsPath']
pipelinesPath = cfg['PATHS']['PipelinesPath']
referencesPath = cfg['PATHS']['ReferencesPath']

if cfg.has_option('PATHS', 'BinPath'):
	binPath = cfg['PATHS']['BinPath']
else:
	binPath = rootDir + '/bin/'
