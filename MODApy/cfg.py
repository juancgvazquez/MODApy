import configparser
import json
import logging
import logging.config
import os

# config parsing from here on, parses paths and things from config.ini
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
testPath = cfg['PATHS']['TestPath']

if cfg.has_option('PATHS', 'BinPath'):
    binPath = cfg['PATHS']['BinPath']
else:
    binPath = rootDir + '/bin/'

if cfg['GENERAL'].getboolean('TestMode'):
    testFlag = True
else:
    testFlag = False


# logging config from here on, parses logcfg.json and all things there
logCfgPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'logcfg.json')
def setup_logging(default_path=logCfgPath, default_level=logging.INFO, env_key='LOG_CFG'):
    """
    Method to configure logging

    """
    path = default_path
    value = os.getenv(env_key, None)
    if value:
        path = value
    if os.path.exists(path):
        with open(path, 'rt') as f:
            config = json.load(f)
        logging.config.dictConfig(config)
    else:
        logging.basicConfig(level=default_level)
