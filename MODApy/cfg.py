import configparser
import json
import logging
import logging.config
import os

from redis import Redis

from rq import Queue


class Config:
    def __init__(self):
        # Default configuration
        self.default_cfg = """
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
        """

        # Config file path
        self.cfgPath = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), "config.ini"
        )

        # config parsing
        self.cfg = configparser.ConfigParser()
        if os.path.exists(self.cfgPath):
            self.cfg.read(self.cfgPath)
        else:
            self.cfg.read_string(self.default_cfg)
            with open(self.cfgPath, "w") as cfgfile:
                self.cfg.write(cfgfile)

        self.rootDir = os.path.dirname(
            (os.path.abspath(__file__))
        )  # Adjusted rootDir calculation

        self.patientPath = self.cfg["PATHS"]["patientpath"]
        self.mitopatientPath = self.cfg["PATHS"]["mitopatientpath"]
        self.panelsPath = self.cfg["PATHS"]["panelspath"]
        self.reportsPath = self.cfg["PATHS"]["reportspath"]
        self.resultsPath = self.cfg["PATHS"]["resultspath"]
        self.pipelinesPath = self.cfg["PATHS"]["pipelinespath"]
        self.referencesPath = self.cfg["PATHS"]["referencespath"]
        self.testPath = self.cfg["PATHS"]["testpath"]
        self.variantsDBPath = self.cfg["PATHS"]["dbpath"]
        self.tmpDir = self.cfg["PATHS"]["tmppath"]

        if self.cfg.has_option("PATHS", "binpath"):
            self.binPath = self.cfg["PATHS"]["binpath"]
        else:
            self.binPath = self.rootDir + "/bin/"

        if self.cfg["GENERAL"].getboolean("testmode"):
            self.testFlag = True
        else:
            self.testFlag = False

        # Redis configuration
        self.redis_conn = self._configure_redis()

        # Queues config
        self.short_queue = Queue(
            name="short_queue", default_timeout=-1, connection=self.redis_conn
        )
        self.long_queue = Queue(
            name="long_queue", default_timeout=-1, connection=self.redis_conn
        )

    def _configure_redis(self):
        if self.cfg.has_option("REDIS", "host"):
            return Redis(self.cfg["REDIS"]["host"])
        else:
            return Redis()

    def setConfig(self, section, key, value):
        if section in self.cfg.sections():
            if key in self.cfg[section]:
                self.cfg[section][key] = value
                logging.info("Changed %s to %s" % (key, value))
            else:
                logging.error(
                    "Key error: Requested %s key is not present in config file" % key
                )
        else:
            logging.error(
                "Section error: Section %s not present in config file" % section
            )
        with open(self.cfgPath, "w") as cfgfile:
            self.cfg.write(cfgfile)

    def setup_logging(self):
        """
        Method to configure logging
        """

        def _touch(path):
            basedir = os.path.dirname(path)
            if not os.path.exists(basedir):
                os.makedirs(basedir)
            with open(path, "a"):
                os.utime(path, None)

        logs_dir = os.path.join(self.rootDir, "logs")
        _touch(os.path.join(logs_dir, "current_run.log"))
        _touch(os.path.join(logs_dir, "info.log"))
        _touch(os.path.join(logs_dir, "errors.log"))
        _touch(os.path.join(logs_dir, "pipe_run.log"))
        if not os.path.exists(os.path.join(logs_dir, "downloads.log")):
            with open(os.path.join(logs_dir, "downloads.log"), "w") as dlog:
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
                    "filename": os.path.join(logs_dir, "pipe_run.log"),
                    "mode": "w",
                },
                "current_run": {
                    "class": "logging.FileHandler",
                    "level": "INFO",
                    "formatter": "current_run",
                    "filename": os.path.join(logs_dir, "current_run.log"),
                    "mode": "w",
                },
                "info_file_handler": {
                    "class": "logging.handlers.RotatingFileHandler",
                    "level": "INFO",
                    "formatter": "simple",
                    "filename": os.path.join(logs_dir, "info.log"),
                    "maxBytes": 10485760,
                    "backupCount": 10,
                    "encoding": "utf8",
                },
                "error_file_handler": {
                    "class": "logging.handlers.RotatingFileHandler",
                    "level": "DEBUG",
                    "formatter": "simple",
                    "filename": os.path.join(logs_dir, "errors.log"),
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


configuration = Config()
