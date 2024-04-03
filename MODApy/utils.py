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


def aminoChange(value: str):
    """
    Given a string `value`, extract the amino acid change from the
    `HGVS.P` format.

    Parameters
    ----------
        value (str): A string representing the `HGVS.P` format.

    Returns
    -------
        str: A string representing the amino acid change or
        "." if not applicable.
    """
    try:
        value = value.replace("p.", "")
        if value[:3] != value[-3:]:
            return "CHANGE"
        else:
            return "."
    except Exception:
        return "."


def divide(x, y):
    """
    Method to divide x on y, needed for dividing freqs.
    Parameters
    ----------
    x
        The dividend
    y
        The divisor
    Returns result or x.
    """
    try:
        return float(x) / y
    except Exception:
        return x
