from fastapi import FastAPI, HTTPException, status
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from typing import Optional
from pydantic import BaseModel
from MODApy import cfg, pipeline, vcfmgr, vcfanalysis
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


class Pipeline(BaseModel):
    Pipeline: str
    FQ_1: str
    FQ_2: Optional[str] = ""
    startStep: Optional[int] = 0
    endStep: Optional[int] = 0
    keeptmp: Optional[bool] = False


class Single(BaseModel):
    patient: str
    panel: str


class Duos(BaseModel):
    patient1: str
    patient2: str
    panel: Optional[str] = None
    vennPlace: Optional[str] = None
    filter: Optional[str] = None


class Trios(BaseModel):
    patient1: str
    patient2: str
    patient3: str
    panel: Optional[str] = None
    vennPlace: Optional[str] = None
    filter: Optional[str] = None


@app.post("/modaapi/single")
async def single(data: Single):
    data = data.dict()
    try:
        panel = data['panel']
        patient = data['patient']
        job_id = cfg.short_queue.enqueue(vcfanalysis.single,
                                         args=[patient, panel])
        return JSONResponse(status_code=status.HTTP_202_ACCEPTED,
                            content=f"Job Queued. Job id is {job_id}")
    except Exception as err:
        logger.error('Api error on Single')
        logger.debug(f"Error was: {err}", exc_info=True)
        raise HTTPException(status_code=404, detail=str(err))


@app.post("/modaapi/duos")
async def duos(data: Duos):
    data = data.dict()
    try:
        patient1 = data['Patient1']
        patient2 = data['Patient2']
        VennPlace = data['vennplace']
        Panel = data['panel']
        Filter = data['filter']
        job_id = cfg.short_queue.enqueue(vcfanalysis.duos,
                                         args=[patient1, patient2],
                                         kwargs={
                                             "VennPlace": VennPlace,
                                             "Panel": Panel,
                                             "Filter": Filter
                                         })
        return JSONResponse(status_code=status.HTTP_202_ACCEPTED,
                            content=f"Job Queued. Job id is {job_id}")
    except Exception as err:
        logger.error('Api error on Duos')
        logger.debug(f"Error was: {err}", exc_info=True)
        raise HTTPException(status_code=404, detail=str(err))


@app.post("/modaapi/trios")
async def trios(data: Trios):
    try:
        patient1 = data['Patient3']
        patient2 = data['Patient3']
        patient3 = data['Patient3']
        VennPlace = data['vennplace']
        Panel = data['panel']
        Filter = data['filter']
        # Checks file existence and type for patients
        job_id = cfg.short_queue.enqueue(vcfanalysis.trios,
                                         args=[patient1, patient2, patient3],
                                         kwargs={
                                             "VennPlace": VennPlace,
                                             "Panel": Panel,
                                             "Filter": Filter
                                         })
        job_id = job_id.job_id
        return JSONResponse(status_code=status.HTTP_202_ACCEPTED,
                            content=f"Job Queued. Job id is {job_id}")
    except Exception as err:
        logger.error('Api error on Trios')
        logger.debug(f"Error was: {err}", exc_info=True)
        raise HTTPException(status_code=404, detail=str(err))


@app.post("/modaapi/pipeline")
async def run_pipeline(data: Pipeline):
    try:
        pipe = data['Pipeline']

        checkFile(pipe, data['Pipeline'].split(".")[-1])

        newpipe = pipeline.Pipeline.from_json(pipe)
        fq1 = data['FQ_1']
        fq2 = data['FQ_2']

        if fq2 != "":

            checkFile(fq1, "." + fq1.split(".")[-1])
            checkFile(fq2, "." + fq2.split(".")[-1])
            job_id = cfg.long_queue.enqueue(
                newpipe.runpipeline,
                args=[fq1, fq2],
                kwargs={
                    "keeptmp": data['keeptmp'],
                    "startStep": data['startStep'],
                    "endStep": data['endStep'],
                },
            )
            job_id = job_id.job_id
            return JSONResponse(status_code=status.HTTP_202_ACCEPTED,
                                content=f"Job Queued. Job id is {job_id}")
        else:
            checkFile(fq1, "." + fq1.split(".")[-1])
            job_id = cfg.long_queue.enqueue(
                newpipe.runpipeline,
                args=[fq1],
                kwargs={
                    "keeptmp": data['keeptmp'],
                    "startStep": data['startStep'],
                    "endStep": data['endStep'],
                },
            )
            job_id = job_id.job_id
            return JSONResponse(status_code=status.HTTP_202_ACCEPTED,
                                content=f"Job Queued. Job id is {job_id}")
    except Exception as err:
        logger.error('Api error on Pipeline')
        logger.debug(f"Error was: {err}", exc_info=True)
        raise HTTPException(status_code=404, detail=str(err))
