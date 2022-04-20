import os
import logging

logger = logging.getLogger(__name__)


def checkFile(filePath, extension):
    if os.path.isfile(filePath):
        fileName, fileExtension = os.path.splitext(filePath)
        if extension == fileExtension:
            return True
    else:
        error = f"""{filePath} couldn't be found. "
        Please check if file exists and that it's extension is {extension}"""
        logger.error(error)
        raise FileNotFoundError(error)
