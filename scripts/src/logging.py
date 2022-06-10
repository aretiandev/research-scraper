import logging
from datetime import datetime
from pathlib import Path


def configLogger(
    name=None,
    logfile=None,
    level="debug",
    file_level="debug",
    console_level="info",
    fmt="%(asctime)s %(name)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
):
    """Configure Root Logger and return named logger."""

    # Get root logger
    logger = logging.getLogger()
    logger.setLevel(level.upper())

    # Configure console output
    sh = logging.StreamHandler()
    sh_formatter = logging.Formatter("%(message)s")
    sh.setFormatter(sh_formatter)
    sh.setLevel(console_level.upper())
    logger.addHandler(sh)

    # Configure file output
    if logfile is not None:
        fh = logging.FileHandler(filename=logfile, mode="a")
        fh_formatter = logging.Formatter(fmt=fmt, datefmt=datefmt)
        fh.setFormatter(fh_formatter)
        fh.setLevel(file_level.upper())
        logger.addHandler(fh)

    # Return named logger
    if name is not None:
        return logging.getLogger(name)
    else:
        return logger
