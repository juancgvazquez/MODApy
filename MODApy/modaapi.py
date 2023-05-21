import logging
from typing import Optional

from MODApy import cfg, pipeline, vcfanalysis
from MODApy.utils import checkFile

from fastapi import FastAPI, HTTPException, status
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from pydantic import BaseModel

logger = logging.getLogger()
app = FastAPI()

origins = [
    "*",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class Pipeline(BaseModel):
    """
    Represents a pipeline.

    Attributes:
        Pipeline (str): The pipeline name.
        FQ_1 (str): The first FastQ file.
        FQ_2 (str, optional): The second FastQ file. Defaults to "".
        startStep (int, optional): The starting step. Defaults to 0.
        endStep (int, optional): The ending step. Defaults to 0.
        keeptmp (bool, optional): Whether to keep temporary files. Defaults to False.
    """

    Pipeline: str
    FQ_1: str
    FQ_2: Optional[str] = ""
    startStep: Optional[int] = 0
    endStep: Optional[int] = 0
    keeptmp: Optional[bool] = False


class Single(BaseModel):
    """
    Represents single input data.

    Attributes:
        patient (str): The patient identifier.
        panel (str): The panel identifier.
    """

    patient: str
    panel: str


class Duos(BaseModel):
    """
    Represents duos input data.

    Attributes:
        patient1 (str): The first patient identifier.
        patient2 (str): The second patient identifier.
        panel (str, optional): The panel identifier. Defaults to None.
        vennPlace (str, optional): The vennPlace identifier. Defaults to None.
        filter (str, optional): The filter identifier. Defaults to None.
    """

    patient1: str
    patient2: str
    panel: Optional[str] = None
    vennPlace: Optional[str] = None
    filter: Optional[str] = None


class Trios(BaseModel):
    """
    Represents trios input data.

    Attributes:
        patient1 (str): The first patient identifier.
        patient2 (str): The second patient identifier.
        patient3 (str): The third patient identifier.
        panel (str, optional): The panel identifier. Defaults to None.
        vennPlace (str, optional): The vennPlace identifier. Defaults to None.
        filter (str, optional): The filter identifier. Defaults to None.
    """

    patient1: str
    patient2: str
    patient3: str
    panel: Optional[str] = None
    vennPlace: Optional[str] = None
    filter: Optional[str] = None


@app.post("/modaapi/single")
async def single(data: Single):
    """
    Handles single input data.

    Parameters:
        data (Single): The single input data.

    Returns:
        JSONResponse: The response containing the job ID.
    """
    data = data.dict()
    try:
        panel = data["panel"]
        patient = data["patient"]
        job_id = cfg.short_queue.enqueue(vcfanalysis.single, args=[patient, panel])
        job_id = job_id.id
        return JSONResponse(
            status_code=status.HTTP_202_ACCEPTED,
            content=f"Job Queued. Job id is {job_id}",
        )
    except Exception as err:
        logger.error("Api error on Single")
        logger.debug(f"Error was: {err}", exc_info=True)
        raise HTTPException(status_code=404, detail=str(err))


@app.post("/modaapi/duos")
async def duos(data: Duos):
    """
    Handles duos input data.

    Parameters:
        data (Duos): The duos input data.

    Returns:
        JSONResponse: The response containing the job ID.
    """
    data = data.dict()
    try:
        patient1 = data["patient1"]
        patient2 = data["patient2"]
        VennPlace = data["vennPlace"]
        Panel = data["panel"]
        Filter = data["filter"]
        job_id = cfg.short_queue.enqueue(
            vcfanalysis.duos,
            args=[patient1, patient2],
            kwargs={"VennPlace": VennPlace, "Panel": Panel, "Filter": Filter},
        )
        job_id = job_id.id
        return JSONResponse(
            status_code=status.HTTP_202_ACCEPTED,
            content=f"Job Queued. Job id is {job_id}",
        )
    except Exception as err:
        logger.error("Api error on Duos")
        logger.debug(f"Error was: {err}", exc_info=True)
        raise HTTPException(status_code=404, detail=str(err))


@app.post("/modaapi/trios")
async def trios(data: Trios):
    """
    Handles trios input data.

    Parameters:
        data (Trios): The trios input data.

    Returns:
        JSONResponse: The response containing the job ID.
    """
    data = data.dict()
    try:
        patient1 = data["patient1"]
        patient2 = data["patient2"]
        patient3 = data["patient3"]
        VennPlace = data["vennPlace"]
        Panel = data["panel"]
        Filter = data["filter"]
        # Checks file existence and type for patients
        job_id = cfg.short_queue.enqueue(
            vcfanalysis.trios,
            args=[patient1, patient2, patient3],
            kwargs={"VennPlace": VennPlace, "Panel": Panel, "Filter": Filter},
        )
        job_id = job_id.id
        return JSONResponse(
            status_code=status.HTTP_202_ACCEPTED,
            content=f"Job Queued. Job id is {job_id}",
        )
    except Exception as err:
        logger.error("Api error on Trios")
        logger.debug(f"Error was: {err}", exc_info=True)
        raise HTTPException(status_code=404, detail=str(err))


@app.post("/modaapi/pipeline")
async def run_pipeline(data: Pipeline):
    """
    Runs the pipeline.

    Parameters:
        data (Pipeline): The pipeline data.

    Returns:
        JSONResponse: The response containing the job ID.
    """
    try:
        data = data.dict()
        pipe = data["Pipeline"]

        checkFile(pipe, data["Pipeline"].split(".")[-1])

        newpipe = pipeline.Pipeline.from_json(pipe)
        fq1 = data["FQ_1"]
        fq2 = data["FQ_2"]

        if fq2 != "":
            checkFile(fq1, "." + fq1.split(".")[-1])
            checkFile(fq2, "." + fq2.split(".")[-1])
            job_id = cfg.long_queue.enqueue(
                newpipe.runpipeline,
                args=[fq1, fq2],
                kwargs={
                    "keeptmp": data["keeptmp"],
                    "startStep": data["startStep"],
                    "endStep": data["endStep"],
                },
            )
            job_id = job_id.id
            return JSONResponse(
                status_code=status.HTTP_202_ACCEPTED,
                content=f"Job Queued. Job id is {job_id}",
            )
        else:
            checkFile(fq1, "." + fq1.split(".")[-1])
            job_id = cfg.long_queue.enqueue(
                newpipe.runpipeline,
                args=[fq1],
                kwargs={
                    "keeptmp": data["keeptmp"],
                    "startStep": data["startStep"],
                    "endStep": data["endStep"],
                },
            )
            job_id = job_id.id
            return JSONResponse(
                status_code=status.HTTP_202_ACCEPTED,
                content=f"Job Queued. Job id is {job_id}",
            )
    except Exception as err:
        logger.error("Api error on Pipeline")
        logger.debug(f"Error was: {err}", exc_info=True)
        raise HTTPException(status_code=404, detail=str(err))
