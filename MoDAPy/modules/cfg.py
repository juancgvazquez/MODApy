import configparser

cfg = configparser.ConfigParser()
cfg.read('config.ini')

patientPath = cfg['PATHS']['PatientPath']
panelsPath = cfg['PATHS']['PanelsPath']
reportsPath = cfg['PATHS']['ReportsPath']
resultsPath = cfg['PATHS']['ResultsPath']