import glob
import logging
import os

from MODApy.cfg import cfg, patientPath, variantsDBPath
from MODApy.vcfmgr import ParsedVCF

import cyvcf2

import numpy as np

import pandas as pd

logger = logging.getLogger(__name__)


# TODO: FIND A MORE EFFICIENT WAY TO SUM EMPTY


class ParquetVarDB(pd.DataFrame):
    @property
    def _constructor(self):
        return ParquetVarDB

    @classmethod
    def from_parquetdb(cls, parquetpath):
        if os.path.exists(parquetpath):
            try:
                db = pd.read_parquet(parquetpath)
            except Exception as e:
                logger.error("There was an error parsing Parquet File")
                logger.debug("", exc_info=True)
                logger.debug(str(e))
                raise e
        else:
            logger.error("Path to Parquet file incorrect.")
            raise FileNotFoundError
        db = db.pipe(ParquetVarDB)
        return db

    @classmethod
    def buildDB(
        cls,
        patientPath=patientPath,
        dbpath=variantsDBPath,
        db=None,
        filetype="vcf",
        prioritized=False,
    ):
        def patientLister(db=db, filetype=filetype):
            vcfspath = []
            for dirpath, dirnames, filenames in os.walk(patientPath):
                for filename in [
                    f for f in filenames if f.lower().endswith(f".final.{filetype}")
                ]:
                    vcfspath.append(os.path.join(dirpath, filename))
            final_list = []
            for idx in range(len(vcfspath)):
                splitfn = vcfspath[idx].rsplit("/", maxsplit=1)[-1]
                if "_MODApy" in splitfn:
                    if not (any(splitfn in string for string in final_list)):
                        final_list.append(vcfspath[idx])
                else:
                    if not (
                        any(
                            splitfn.strip(f".final.{filetype}") in string
                            for string in final_list
                        )
                    ):
                        final_list.append(vcfspath[idx])
            for x in final_list:
                fn = x.rsplit("/", maxsplit=1)[-1]
                if any(
                    fn.strip(f".final.{filetype}") + "_MODApy.final.{filetype}"
                    in string
                    for string in final_list
                ):
                    final_list.remove(x)
            vcfspath = final_list

            if db is not None:
                if filetype == "vcf":
                    try:
                        vcfsnames = [cyvcf2.Reader(x).samples[0] for x in vcfspath]
                    except Exception as e:
                        logger.info(
                            "No Sample name in one of the vcfs files. Using File Names Instead"
                        )
                        logger.debug(str(e))
                        vcfsnames = [
                            x.rsplit("/", maxsplit=1)[-1].strip(f".final.{filetype}")
                            for x in vcfspath
                        ]
                elif filetype == "parquet":
                    try:
                        vcfsnames = [
                            pd.read_parquet(x, columns='SAMPLE', nrows=1).columns[
                                'SAMPLE'
                            ][0]
                            for x in vcfspath
                        ]
                    except Exception as e:
                        logger.info(
                            "No Sample name in one of the vcfs files. Using File Names Instead"
                        )
                        logger.debug(str(e))
                        vcfsnames = [
                            x.rsplit("/", maxsplit=1)[-1].strip(f".final.{filetype}")
                            for x in vcfspath
                        ]
                addpatnames = [
                    x
                    for x in vcfsnames
                    if (
                        x not in db.SAMPLE.unique().astype(str).tolist()
                        and x + "_MODApy" not in db.SAMPLE.unique().astype(str).tolist()
                        and x.replace("_MODApy", "")
                        not in db.SAMPLE.unique().astype(str).tolist()
                    )
                ]
                if len(addpatnames) >= 1:
                    logger.info("Adding Patients: {}".format([x for x in addpatnames]))
                else:
                    logger.error("No Patients to Add")
                    exit(1)
                patientslist = [x for x in vcfspath for y in addpatnames if y in x]
            else:
                patientslist = vcfspath

            return patientslist

        def dbbuilder(patientslist, db=db, prioritized=prioritized):
            logger.info("Parsing Patients")
            pvcfs = ParsedVCF.mp_parser(*patientslist, prioritized=prioritized)
            logger.info("Patients Parsed")
            logger.info("Building Database")
            if not os.path.exists(dbpath):
                os.makedirs(dbpath)
            if not os.path.exists(dbpath):
                pvcfs[0].vcf_to_parquet(dbpath, partition_cols=['SAMPLE'], append=False)
                [
                    x.vcf_to_parquet(
                        dbpath,
                        partition_cols=['SAMPLE'],
                        append=True,
                    )
                    for x in pvcfs[1:]
                ]
            else:
                [
                    x.vcf_to_parquet(
                        dbpath,
                        partition_cols=['SAMPLE'],
                        append=True,
                    )
                    for x in pvcfs
                ]
            logger.info("Database Built")

        try:
            logger.info("Checking DB File")
            patientslist = patientLister(db=db)
        except Exception as e:
            logger.debug(str(e))
            exit()
        sublists = [
            patientslist[i : i + int(cfg["GENERAL"]["cores"])]
            for i in range(0, len(patientslist), int(cfg["GENERAL"]["cores"]))
        ]
        for lista in sublists:
            db = dbbuilder(lista, db, prioritized=prioritized)
        return db
