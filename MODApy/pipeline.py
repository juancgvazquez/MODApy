import json
import logging
import os
import shlex
import shutil
import subprocess

from MODApy import configuration, vcfmgr

import xmltodict

import yaml


logger = logging.getLogger(__name__)
logger2 = logging.getLogger("Pipeline Module")
os.makedirs(configuration.rootDir + "/logs", exist_ok=True)
hdlr = logging.FileHandler(configuration.rootDir + "/logs/pipe_run.log")
formatter = logging.Formatter("%(asctime)s %(name)-25s %(levelname)-8s %(message)s")
hdlr.setFormatter(formatter)
logger2.addHandler(hdlr)
logger2.setLevel(logging.DEBUG)


# TODO: restructure pipeline, being less open, so less prone to errors


class PipeStep(object):
    """
    Class defined for each step of the Pipeline
    """

    def __init__(self, name, command, subcommand, version, inputfile, outputfile, args):
        """
        Initializes a PipeStep object.

        Parameters
        ----------
        name : str
            Name of the step.
        command : str
            Command to be executed for the step.
        subcommand : str
            Subcommand to be executed for the step.
        version : str
            Version of the command.
        inputfile : str or list
            Path to the input file(s).
        outputfile : str
            Path to the output file.
        args : str
            Additional arguments for the command.
        """
        self.name = name
        self.command = command
        self.subcommand = subcommand
        self.version = version
        self.inputfile = inputfile
        self.outputfile = outputfile
        self.args = args

    def __str__(self):
        """
        Returns the name of the step as a string.

        Returns
        -------
        str
            Name of the step.
        """
        return self.name

    def __repr__(self):
        """
        Returns the name of the step as a string.

        Returns
        -------
        str
            Name of the step.
        """
        return self.name

    def stepinfo(self):
        """
        Prints information about the step.
        """
        logger.info("Name:", self.name)
        logger.info("Command", self.command + self.version, self.subcommand, self.args)
        logger.info("Input File:", self.inputfile)
        logger.info("Output File:", self.outputfile)


class Pipeline(object):
    """
    General Class for Pipelines
    """

    def __init__(self, name, reference, url="", description=""):
        """
        Initializes a Pipeline object.

        Parameters
        ----------
        name : str
            Name of the pipeline.
        reference : str
            Reference for the pipeline.
        url : str, optional
            URL for the pipeline, by default "".
        description : str, optional
            Description of the pipeline, by default "".
        """
        self.name = name
        self.url = url
        self.description = description
        self.reference = reference
        self.required_files = []
        self.steps = []

    def add_steps(self, step):
        """
        Adds a step to the pipeline.

        Parameters
        ----------
        step : PipeStep
            Step to be added to the pipeline.
        """
        self.steps.append(step)

    def add_req_files(self, path):
        """
        Adds a required file to the pipeline.

        Parameters
        ----------
        path : str
            Path to the required file.
        """
        self.required_files.append(path)

    @staticmethod
    def _buildpipe(pipedict):
        """
        Private Class Method to build pipeline from loaded json,xml or yaml

        Parameters
        ----------
        pipedict : dict
            Dictionary representing the pipeline.
        """
        name = pipedict["INFO"]["name"]
        url = pipedict["INFO"]["url"]
        description = pipedict["INFO"]["description"]
        reference = pipedict["INFO"]["reference"]

        newpipe = Pipeline(name, reference, url, description)

        for x in pipedict["INFO"]["required_files"]:
            newpipe.add_req_files(x)

        steps = list(pipedict["STEPS"].values())

        for i in range(0, len(steps)):
            name = steps[i]["name"]
            command = steps[i]["command"]
            subcommand = steps[i]["subcommand"]
            version = steps[i]["version"]
            inputfile = steps[i]["input"]
            outputfile = steps[i]["output"]
            args = steps[i]["args"]
            newstep = PipeStep(
                name, command, subcommand, version, inputfile, outputfile, args
            )
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
            pipedict = json.load(f)

        builtpipe = cls._buildpipe(pipedict)

        return builtpipe

    @classmethod
    def from_yaml(cls, yamlpath):
        """
        Class Method to read yaml and build the Pipeline

        Parameters
        ----------
        jsonpath
            Path to the yaml file used to create the Pipeline.

        Returns Pipeline object
        """
        with open(yamlpath) as f:
            pipedict = yaml.load(f)
        builtpipe = cls._buildpipe(pipedict)

        return builtpipe

    @classmethod
    def from_xml(cls, xmlpath):
        """
        Class Method to parse xml and build the Pipeline

        Parameters
        ----------
        xmlpath
            Path to the xml file used to create the Pipeline.

        Returns Pipeline object
        """
        with open(xmlpath) as f:
            pipedict = xmltodict.parse(f.read())["root"]
        builtpipe = cls._buildpipe(pipedict)

        return builtpipe

    def runpipeline(
        self,
        fastq1: str,
        fastq2=None,
        keeptmp=False,
        startStep=0,
        endStep=0,
        patientPath=None,
    ):
        """
        Method to run the Pipeline

        Parameters
        ----------
        fastq1
            Path to the first fastq file.
        fastq2
            Path to the second fastq file, in case of paired reads.
        """

        try:
            logger.info(self.steps)
            logger.info("Nro de Pasos: %s" % str(len(self.steps)))
            logger2.info(self.steps)
            logger2.info("Nro de Pasos: %s" % str(len(self.steps)))
            patientname = fastq1.split("/")[-1].split(".")[0].split("_1")[0]
            ref = (
                configuration.referencesPath
                + self.reference
                + "/"
                + self.reference
                + ".fa"
            )
            pipedir = "".join(x for x in self.name if x.isalnum())
            if patientPath is None:
                patientPath = configuration.patientPath
            if configuration.testFlag:
                tmpdir = (
                    configuration.testPath
                    + "Pipelines/"
                    + patientname
                    + "/"
                    + pipedir
                    + "/tmp/"
                )
            else:
                tmpdir = (
                    configuration.resultsPath
                    + "Pipelines/"
                    + patientname
                    + "/"
                    + pipedir
                    + "/tmp/"
                )

            os.makedirs(tmpdir, exist_ok=True)
            samplename = patientname + "_MODApy"
            logger2.info(
                "Running "
                + str(self.name)
                + " pipeline on patient: "
                + str(patientname)
            )
            # bool to check if first step
            first = True
            if endStep == 0:
                endStep = len(self.steps) + 1
            if startStep > 0:
                first = False
            logger2.info(f"STARTSTEP {startStep}")
            for step in self.steps[startStep:endStep]:
                if first is True:
                    logger2.debug("First Step")
                    first = False
                    if isinstance(step.inputfile, list):
                        if any(
                            x in y
                            for y in [fastq1, fastq2]
                            for x in [".fastq", ".fastq.gz", ".fq", ".fq.gz"]
                        ):
                            if (fastq1 is not None) & (fastq2 is not None):
                                inputfile = fastq1 + " " + fastq2
                            else:
                                logger.warn(
                                    "WARNING: This pipeline was designed for Pair End \
                                    and you are running it as Single End"
                                )
                                inputfile = fastq1
                    elif isinstance(step.inputfile, str):
                        if any(
                            x in y
                            for y in [fastq1]
                            for x in [".fastq", ".fastq.gz", ".fq", ".fq.gz"]
                        ):
                            if (fastq1 is not None) & (fastq2 is not None):
                                logger.warn(
                                    "WARNING: This pipeline was designed for Single \
                                        End and you are running it as Pair End"
                                )
                                inputfile = fastq1 + " " + fastq2
                        else:
                            inputfile = fastq1
                # If it's not first step, input depends on output of previous
                # step + patientname
                else:
                    if isinstance(step.inputfile, str):
                        inputfile1 = step.inputfile[0].replace(
                            "patientname",
                            tmpdir + patientname + "/" + patientname,
                        )
                        inputfile2 = step.inputfile[1].replace(
                            "patientname",
                            tmpdir + patientname + "/" + patientname,
                        )
                        inputfile = inputfile1 + " " + inputfile2
                    elif isinstance(step.inputfile, str):
                        inputfile = step.inputfile.replace(
                            "patientname", tmpdir + patientname
                        )
                # replaces patient name in outputfiles
                if isinstance(step.outputfile, str):
                    outputfile = step.outputfile.replace(
                        "patientname", tmpdir + patientname
                    )
                else:
                    return "Error Parsing output file. It should be a string."

                logger2.info(step.name)
                args = (
                    step.args.replace("patientname", tmpdir + patientname)
                    .replace("reference", ref)
                    .replace("samplename", samplename)
                )
                cmdver = step.version.replace(".", "_")
                javacmds = ["GATK", "picard", "SnpSift", "snpEff"]
                if any(javacmd in step.command for javacmd in javacmds):
                    cmd = (
                        "java -jar -Xmx12G -Djava.io.tmpdir=%s " % "~/.tmp"
                        + configuration.binPath
                        + step.command
                        + "/"
                        + step.command
                        + "_"
                        + cmdver
                        + ".jar "
                        + step.subcommand
                    )
                else:
                    cmd = (
                        configuration.binPath
                        + step.command
                        + "/"
                        + step.command
                        + "_"
                        + cmdver
                        + " "
                        + step.subcommand
                    )
                if "HaplotypeCaller" in cmd:
                    cmdstr = cmd + " " + args + " " + inputfile + " " + outputfile
                else:
                    cmdstr = cmd + " " + args + " " + " " + inputfile + " " + outputfile
                    logger.info(f"Command is: {cmd}")
                cmd = shlex.split(cmdstr)

                logging.info("Subprocess: " + cmdstr)
                logger2.info("Subprocess: " + cmdstr)
                stdcmds = ["bwa", "bedtools", "snpEff", "SnpSift"]
                logger.info(f"Command is: {cmd}")
                try:
                    if any(stdcmd in s for s in cmd for stdcmd in stdcmds):
                        output = cmd[-1]
                        del cmd[-1]
                        with open(output, "w+") as out:
                            cmdrun = subprocess.Popen(
                                cmd,
                                stderr=subprocess.PIPE,
                                stdout=out,
                                universal_newlines=True,
                            )
                    else:
                        cmdrun = subprocess.Popen(
                            cmd,
                            stderr=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            universal_newlines=True,
                        )
                    out, err = cmdrun.communicate()
                    if out is not None:
                        logging.debug("Output: " + out.strip())
                        logger2.debug("Output: " + out.strip())
                    if err is not None:
                        logging.debug("Stderr: " + err.strip())
                        logger2.debug("Stderr: " + err.strip())
                    if cmdrun.returncode != 0:
                        logging.error(
                            "Subprocess failed with error code: "
                            + str(cmdrun.returncode)
                        )
                        logger2.error(
                            "Subprocess failed with error code: "
                            + str(cmdrun.returncode)
                        )
                        logging.error("Check log for more details")
                        logger2.error("Check log for more details")
                        exit(cmdrun.returncode)

                except (OSError, subprocess.CalledProcessError) as exception:
                    logging.debug("Subprocess failed")
                    logging.debug("Exception ocurred: " + str(exception))
                    logger2.debug("Subprocess failed")
                    logger2.debug("Exception ocurred: " + str(exception))
                    logging.info(
                        "There was an error when running the pipeline. Please \
                            check logs for more info"
                    )
                    logger2.info(
                        "There was an error when running the pipeline. Please \
                            check logs for more info"
                    )
                    exit(1)
                else:
                    logging.info("Subprocess finished")
                    logger2.info("Subprocess finished")
            if configuration.testFlag:
                if os.path.exists(tmpdir + patientname + "_MODApy.final.vcf"):
                    file = (
                        configuration.testPath
                        + patientname
                        + "_MODApy/"
                        + patientname
                        + "_MODApy.final.vcf"
                    )
                    os.makedirs(
                        configuration.testPath + patientname + "_MODApy", exist_ok=True
                    )
                    shutil.move(tmpdir + patientname + "_MODApy.final.vcf", file)
                    logger2.info("Parsing final VCF file")
                    logging.info("Parsing final VCF file")
                    vcfmgr.ParsedVCF.from_vcf(file).to_csv(
                        file.split(".vcf")[0] + ".csv", index=False
                    )
                if os.path.exists(tmpdir + patientname + "_realigned_reads_recal.bam"):
                    shutil.move(
                        tmpdir + patientname + "_realigned_reads_recal.bam",
                        configuration.testPath
                        + patientname
                        + "_MODApy/"
                        + patientname
                        + "MODApy_realigned_reads_recal.bam",
                    )
                if os.path.exists(tmpdir + patientname + "_realigned_reads_recal.bai"):
                    shutil.move(
                        tmpdir + patientname + "_realigned_reads_recal.bai",
                        configuration.testPath
                        + patientname
                        + "_MODApy/"
                        + patientname
                        + "MODApy_realigned_reads_recal.bai",
                    )
            else:
                if os.path.exists(tmpdir + patientname + "_MODApy.final.vcf"):
                    os.makedirs(patientPath + patientname + "_MODApy", exist_ok=True)
                    file = (
                        patientPath
                        + patientname
                        + "_MODApy/"
                        + patientname
                        + "_MODApy.final.vcf"
                    )
                    shutil.move(tmpdir + patientname + "_MODApy.final.vcf", file)
                    logger2.info("Parsing final VCF file")
                    logging.info("Parsing final VCF file")
                    vcfmgr.ParsedVCF.from_vcf(file).to_csv(
                        file.split(".vcf")[0] + ".csv", index=False
                    )
                if os.path.exists(tmpdir + patientname + "_realigned_reads_recal.bai"):
                    shutil.move(
                        tmpdir + patientname + "_realigned_reads_recal.bai",
                        patientPath
                        + patientname
                        + "_MODApy/"
                        + patientname
                        + "MODApy_realigned_reads_recal.bai",
                    )
                if os.path.exists(tmpdir + patientname + "_realigned_reads_recal.bam"):
                    shutil.move(
                        tmpdir + patientname + "_realigned_reads_recal.bam",
                        patientPath
                        + patientname
                        + "_MODApy/"
                        + patientname
                        + "MODApy_realigned_reads_recal.bam",
                    )
            if keeptmp is False:
                shutil.rmtree(tmpdir)
            logger2.info("Pipeline completed!")
            logging.info("Pipeline completed!")
        except Exception as error:
            logger2.info("Pipeline Failed!")
            logging.info("Pipeline Failed!")
            logger2.debug(f"Pipeline error: {error}", exc_info=True)
            logging.debug(f"Pipeline error: {error}", exc_info=True)

    def pipelineinfo(self):
        """
        Method to print Pipeline Info
        """
        logger.info("Name:", self.name)
        logger.info("Reference:", self.reference)
        logger.info("URL:", self.url)
        logger.info("Description:", self.description)
        logger.info("Additional Files Required:", self.required_files)
        logger.info("Steps:", self.steps)
