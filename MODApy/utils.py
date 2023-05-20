import logging
import os

logger = logging.getLogger(__name__)


class InvalidFileError(Exception):
    pass


def checkFile(filePath, extension):
    if os.path.isfile(filePath):
        fileName, fileExtension = os.path.splitext(filePath)
        if extension == fileExtension:
            return True
    error = f"""{filePath} couldn't be found. "
    Please check if the file exists and that its extension is {extension}"""
    logger.error(error)
    raise InvalidFileError(error)
