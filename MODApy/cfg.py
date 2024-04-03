import configparser
import json
import logging
import logging.config
import os

from redis import Redis

from rq import Queue


default_cfg = """
[GENERAL]
testmode = False
cores = 1
processing_mode = local

[PATHS]
patientpath = ./Patients/
mitopatientpath = ./mitocondrial/Patients/
panelspath = ./Panels/
reportspath = ./Reports/
resultspath = ./Results/
pipelinespath = ./Pipelines/
referencespath = ./References/
testpath = ./test/
dbpath = ./VariantsDB/variantsDB.csv
tmppath = ./tmp/
binpath = ./bin/

[OUTPUT]
columnsorder = GENE_NAME, AMINOCHANGE, HGVS.P, HGVS.C, RSID, IMPACT, EFFECT,
               FEATURE_ID, VARDB_FREQ, ALLELE_FREQ, 1000GP3_AF, 1000GP3_AFR_AF,
               1000GP3_AMR_AF, 1000GP3_EAS_AF, 1000GP3_EUR_AF, 1000GP3_SAS_AF,
               ESP6500_MAF_EA, ESP6500_MAF_AA, ESP6500_MAF_ALL, CLINVAR_CLNSIG,
               CLINVAR_CLNDSDB,CLINVAR_CLNDSDBID, CLINVAR_CLNDBN, CLINVAR_CLNREVSTAT,
               CLINVAR_CLNACC, POLYPHEN_PRED,POLYPHEN_SCORE, VENN, ZIGOSITY, CHROM,
               POS, REF, ALT, VARTYPE
"""


# config parsing from here on, parses paths and things from config.ini
cfgPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "config.ini")
cfg = configparser.ConfigParser()
if os.path.exists(cfgPath):
    cfg.read(cfgPath)
else:
    cfg.read_string(default_cfg)
    with open(cfgPath, "w") as cfgfile:
        cfg.write(cfgfile)

rootDir = os.path.dirname(os.path.abspath(__file__))

patientPath = cfg["PATHS"]["patientpath"]
mitopatientPath = cfg["PATHS"]["mitopatientpath"]
panelsPath = cfg["PATHS"]["panelspath"]
reportsPath = cfg["PATHS"]["reportspath"]
resultsPath = cfg["PATHS"]["resultspath"]
pipelinesPath = cfg["PATHS"]["pipelinespath"]
referencesPath = cfg["PATHS"]["referencespath"]
testPath = cfg["PATHS"]["testpath"]
variantsDBPath = cfg["PATHS"]["dbpath"]
tmpDir = cfg["PATHS"]["tmppath"]

if cfg.has_option("PATHS", "binpath"):
    binPath = cfg["PATHS"]["binpath"]
else:
    binPath = rootDir + "/bin/"

if cfg["GENERAL"].getboolean("testmode"):
    testFlag = True
else:
    testFlag = False

# queues config
if cfg.has_option("REDIS", "host"):
    redis_conn = Redis(cfg["REDIS"]["host"])
else:
    redis_conn = Redis()
short_queue = Queue(name="short_queue", default_timeout=-1, connection=redis_conn)
long_queue = Queue(name="long_queue", default_timeout=-1, connection=redis_conn)


def setConfig(section, key, value):
    if section in cfg.sections():
        if key in cfg[section]:
            cfg[section][key] = value
            logging.info("Changed %s to %s" % (key, value))
        else:
            logging.error(
                "Key error: Requested %s key is not present in config file" % key
            )
    else:
        logging.error("Section error: Section %s not present in config file" % section)
    with open(cfgPath, "w") as cfgfile:
        cfg.write(cfgfile)


def setup_logging():
    """
    Method to configure logging

    """

    def _touch(path):
        basedir = os.path.dirname(path)
        if not os.path.exists(basedir):
            os.makedirs(basedir)
        with open(path, "a"):
            os.utime(path, None)

    _touch(rootDir + "/logs/current_run.log")
    _touch(rootDir + "/logs/info.log")
    _touch(rootDir + "/logs/errors.log")
    _touch(rootDir + "/logs/pipe_run.log")
    if not os.path.exists(rootDir + "/logs/downloads.log"):
        with open(rootDir + "/logs/downloads.log", "w") as dlog:
            json.dump({}, dlog)
    logCfg = {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "simple": {
                "format": "%(asctime)s %(name)-25s %(levelname)-8s %(message)s",
                "datefmt": "%m-%d %H:%M",
            },
            "console": {"format": "%(name)-25s: %(levelname)-8s %(message)s"},
            "current_run": {"format": "%(message)s"},
        },
        "handlers": {
            "console": {
                "class": "logging.StreamHandler",
                "level": "INFO",
                "formatter": "console",
                "stream": "ext://sys.stdout",
            },
            "pipe_run": {
                "class": "logging.FileHandler",
                "level": "INFO",
                "formatter": "current_run",
                "filename": rootDir + "/logs/pipe_run.log",
                "mode": "w",
            },
            "current_run": {
                "class": "logging.FileHandler",
                "level": "INFO",
                "formatter": "current_run",
                "filename": rootDir + "/logs/current_run.log",
                "mode": "w",
            },
            "info_file_handler": {
                "class": "logging.handlers.RotatingFileHandler",
                "level": "INFO",
                "formatter": "simple",
                "filename": rootDir + "/logs/info.log",
                "maxBytes": 10485760,
                "backupCount": 10,
                "encoding": "utf8",
            },
            "error_file_handler": {
                "class": "logging.handlers.RotatingFileHandler",
                "level": "DEBUG",
                "formatter": "simple",
                "filename": rootDir + "/logs/errors.log",
                "maxBytes": 10485760,
                "backupCount": 10,
                "encoding": "utf8",
            },
        },
        "root": {
            "level": "DEBUG",
            "handlers": [
                "console",
                "current_run",
                "info_file_handler",
                "error_file_handler",
            ],
        },
    }

    logging.config.dictConfig(logCfg)
