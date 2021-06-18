import os
import logging

logger = logging.getLogger(__name__)


def checkFile(filePath, extension):
    if os.path.isfile(filePath):
        fileName, fileExtension = os.path.splitext(filePath)
        if extension == fileExtension:
            return True
    else:
        logger.error(
            "%s couldn't be found. Please check if file exists and that it's extension is %s"
            % (filePath, extension)
        )
        exit(1)
