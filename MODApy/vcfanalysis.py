import logging

from MODApy import vcfmgr, cfg
from MODApy.utils import checkFile
import os

logger = logging.getLogger()


def single(patient, panel):
    try:
        checkFile(patient, ".vcf")
        checkFile(panel, ".xlsx")
        logger.info("Running %s on patient %s" % (str(panel), str(patient)))
        result = vcfmgr.ParsedVCF.from_vcf(patient).panel(panel)
        outpath = (
                cfg.patientPath
                + result.name
                + "/Panels/"
                + result.name
                + "_"
                + panel
                + ".xlsx"
        )
        os.makedirs(os.path.dirname(outpath), exist_ok=True)
        result.vcf_to_excel(outpath)
        logger.info("Single Analisis Complete")
        logger.info("File available at:%s" % outpath)
        return outpath
    except Exception as err:
        logger.error("Single analysis Failed")
        logger.debug(f"Error was: {err}", exc_info=True)
        raise Exception


def duos(patient1, patient2, VennPlace=None, Panel=None, Filter=None):
    try:
        checkFile(patient1, ".vcf")
        checkFile(patient2, ".vcf")
        logger.info(
            "Running Duos Study on %s and %s" % (str(patient1), str(patient2)))
        pvcfs = vcfmgr.ParsedVCF.mp_parser(patient1, patient2)
        result = pvcfs[0].duos(pvcfs[1], VENNPLACE=VennPlace)
        resultname = result.name
        outpath = (
                cfg.resultsPath
                + "Duos/"
                + result.name.replace(":", "_")
                + "/"
                + result.name.replace(":", "_")
        )
        result.name = resultname
        if VennPlace is not None:
            outpath = outpath + "_Venn" + VennPlace.replace(":", "_")
        if Panel is not None:
            logger.info("Running panel {}".format(Panel))
            panel = Panel
            checkFile(panel, ".xlsx")
            result = result.panel(panel)
            result.name = resultname
            outpath = outpath + "_P" + Panel
        if Filter[0] is not None:
            for x in Filter:
                if (len(x.split())) != 2:
                    logger.error(
                        "--Filter accepts only two arguments. \
                            Usage: --Filter COLUMN_NAME TEXT_TO_FILTER"
                    )
                    exit(1)
                else:
                    x = x.split()
                    if x[1] == "Empty":
                        result = result[result[x[0]] != ""]
                    else:
                        result = result[~result[x[0]].str.contains(x[1])]
                result.name = resultname
                outpath = outpath + "_F" + str(x[0]) + str(x[1])
        outpath = outpath + ".xlsx"
        logger.info("Writing Result File")
        result.vcf_to_excel(outpath)
        logger.info("Duos Analisis Complete")
        logger.info("File available at:%s" % outpath)
    except Exception as err:
        logger.error("Duos Analisis Failed")
        logger.debug(f"Error was: {err}", exc_info=True)
        raise Exception


def Trios(patient1, patient2, patient3, VennPlace=None, Filter=None,
          Panel=None):
    try:
        checkFile(patient1, ".vcf")
        checkFile(patient2, ".vcf")
        checkFile(patient3, ".vcf")
        logger.info(
            "Running Trios Study on %s, %s and %s"
            % (str(patient1), str(patient2), str(patient3))
        )
        pvcfs = vcfmgr.ParsedVCF.mp_parser(patient1, patient2, patient3)
        result = pvcfs[0].duos(pvcfs[1]).duos(pvcfs[2], VENNPLACE=VennPlace)
        resultname = result.name
        outpath = (
                cfg.resultsPath
                + "Trios/"
                + result.name.replace(":", "_")
                + "/"
                + result.name.replace(":", "_")
        )
        result.name = resultname
        if VennPlace is not None:
            outpath = outpath + "_Venn" + VennPlace.replace(":", "_")
        # check if there is a Panel Requested
        if Panel:
            logger.info("Running panel {}".format(Panel))
            panel = Panel
            checkFile(panel, ".xlsx")
            result = result.panel(panel)
            result.name = resultname
            outpath = outpath + "_Panel" + Panel
        # check if there is a Filter Requested
        if Filter[0] is not None:
            for x in Filter:
                if (len(x.split())) != 2:
                    logger.error(
                        "--Filter accepts only two arguments. \
                            Usage: --Filter COLUMN_NAME TEXT_TO_FILTER"
                    )
                    exit(1)
                else:
                    x = x.split()
                    if x[1] == "Empty":
                        result = result[result[x[0]] != ""]
                    else:
                        result = result[~result[x[0]].str.contains(x[1])]
                    result.name = resultname
                    outpath = outpath + "_Filter" + str(x[0]) + str(x[1])
        outpath = outpath + ".xlsx"
        result.vcf_to_excel(outpath)
        logger.info("Trios Analisis Complete")
        logger.info("File available at:%s" % outpath)
    except Exception as err:
        logger.error("Trios Analisis Failed")
        logger.debug(f"Error was: {err}", exc_info=True)
        raise Exception
