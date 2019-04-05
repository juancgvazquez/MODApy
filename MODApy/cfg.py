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

patientPath = cfg['PATHS']['patientpath']
panelsPath = cfg['PATHS']['panelspath']
reportsPath = cfg['PATHS']['reportspath']
resultsPath = cfg['PATHS']['resultspath']
pipelinesPath = cfg['PATHS']['pipelinespath']
referencesPath = cfg['PATHS']['referencespath']
testPath = cfg['PATHS']['testpath']
variantsDBPath = cfg['PATHS']['dbpath']
tmpDir = cfg['PATHS']['tmppath']

if cfg.has_option('PATHS', 'binpath'):
    binPath = cfg['PATHS']['binpath']
else:
    binPath = rootDir + '/bin/'

if cfg['GENERAL'].getboolean('testmode'):
    testFlag = True
else:
    testFlag = False


def setConfig(section, key, value):
    if section in cfg.sections():
        if key in cfg[section]:
            cfg[section][key] = value
            logging.info('Changed %s to %s' % (key, value))
        else:
            logging.error('Key error: Requested %s key is not present in config file' % key)
    else:
        logging.error('Section error: Section %s not present in config file' % section)
    with open(cfgPath, 'w') as cfgfile:
        cfg.write(cfgfile)


def setup_logging():
    """
    Method to configure logging

    """

    def _touch(path):
        basedir = os.path.dirname(path)
        if not os.path.exists(basedir):
            os.makedirs(basedir)
        with open(path, 'a'):
            os.utime(path, None)

    _touch(rootDir + '/logs/currentrun.log')
    _touch(rootDir + '/logs/info.log')
    _touch(rootDir + '/logs/errors.log')
    if not os.path.exists(rootDir + '/logs/downloads.log'):
        with open(rootDir + '/logs/downloads.log', 'w') as dlog:
            json.dump({}, dlog)
    logCfg = {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "simple": {
                "format": "%(asctime)s %(name)-25s %(levelname)-8s %(message)s",
                "datefmt": "%m-%d %H:%M"
            },
            "console": {
                "format": "%(name)-25s: %(levelname)-8s %(message)s"
            },
            "current_run": {
                "format": "%(message)s"
            }
        },
        "handlers": {
            "console": {
                "class": "logging.StreamHandler",
                "level": "INFO",
                "formatter": "console",
                "stream": "ext://sys.stdout"
            },
            "current_run": {
                "class": "logging.FileHandler",
                "level": "INFO",
                "formatter": "current_run",
                "filename": rootDir + "/logs/currentrun.log",
                "mode": 'w'
            },
            "info_file_handler": {
                "class": "logging.handlers.RotatingFileHandler",
                "level": "INFO",
                "formatter": "simple",
                "filename": rootDir + "/logs/info.log",
                "maxBytes": 10485760,
                "backupCount": 10,
                "encoding": "utf8"
            },
            "error_file_handler": {
                "class": "logging.handlers.RotatingFileHandler",
                "level": "DEBUG",
                "formatter": "simple",
                "filename": rootDir + "/logs/errors.log",
                "maxBytes": 10485760,
                "backupCount": 10,
                "encoding": "utf8"

            }
        },
        "root": {
            "level": "DEBUG",
            "handlers": [
                "console",
                "current_run",
                "info_file_handler",
                "error_file_handler"]
        }
    }

    logging.config.dictConfig(logCfg)
