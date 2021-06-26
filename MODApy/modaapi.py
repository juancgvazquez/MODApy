from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from MODApy import cfg, pipeline, vcfmgr
from MODApy.utils import checkFile
import logging
import os

logger = logging.getLogger()
app = FastAPI()

origins = [
    "*",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=['*'],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class Single(BaseModel):
    patient: str
    panel: str


@app.post("/modaapi/single")
async def single(data: Single):
    data = data.dict()
    try:
        panel = data['panel']
        patient = data['patient']
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
        return 200
    except Exception as err:
        logger.info("Single analysis Failed")
        logger.debug(f"Error was: {err}")
        raise HTTPException(status_code=404, detail=str(err))


@app.post("/modaapi/duos")
async def duos(
        Patient1: str,
        Patient2: str,
        Panel: str = None,
        VennPlace: str = None,
        Filter: str = None,
):
    try:
        patient1 = Patient1
        patient2 = Patient2
        # Checks file existence and type for patients
        pt1Check = checkFile(patient1, ".vcf")
        pt2Check = checkFile(patient2, ".vcf")
        logger.info(
            "Running Duos Study on %s and %s" % (str(Patient1), str(Patient2)))
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
    except:
        logger.info("Duos Analisis Failed")
    return 0


@app.post("/modaapi/trios")
async def trios(
        Patient1: str,
        Patient2: str,
        Patient3: str,
        Panel: str = None,
        VennPlace: str = None,
        Filter: str = None,
):
    try:
        patient1 = Patient1
        patient2 = Patient2
        patient3 = Patient3
        # Checks file existence and type for patients
        pt1Check = checkFile(patient1, ".vcf")
        pt2Check = checkFile(patient2, ".vcf")
        pt3Check = checkFile(patient3, ".vcf")
        logger.info(
            "Running Trios Study on %s, %s and %s"
            % (str(Patient1), str(Patient2), str(Patient3))
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
    except:
        logger.info("Trios Analisis Failed")
    return 0


@app.post("/modaapi/pipeline")
async def run_pipeline(
        Pipeline: str,
        FQ_1: str,
        FQ_2: str = "",
        startStep: int = 0,
        endStep: int = 0,
        keeptmp: bool = False,
):
    pipe = Pipeline

    checkFile(pipe, Pipeline.split(".")[-1])

    newpipe = pipeline.Pipeline.from_json(pipe)

    if FQ_2 != "":
        fq1 = FQ_1
        fq2 = FQ_2
        checkFile(fq1, "." + fq1.split(".")[-1])
        checkFile(fq2, "." + fq2.split(".")[-1])
        cfg.long_queue.enqueue(
            newpipe.runpipeline,
            args=[fq1, fq2],
            kwargs={
                "keeptmp": keeptmp,
                "startStep": startStep,
                "endStep": endStep,
            },
        )
        return 0
    else:
        fq1 = FQ_1
        checkFile(fq1, "." + fq1.split(".")[-1])
        cfg.long_queue.enqueue(
            newpipe.runpipeline,
            args=[fq1],
            kwargs={
                "keeptmp": keeptmp,
                "startStep": startStep,
                "endStep": endStep,
            },
        )
        return 0
