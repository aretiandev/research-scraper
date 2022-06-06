import logging
from datetime import datetime
from pathlib import Path


def configLogger(
    name=None,
    filename=None,
    level="debug",
    file_level="debug",
    console_level="info",
    fmt="%(asctime)s %(name)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
):
    """Configure Logger

    If name=None, configures the root logger."""

    # Get root logger
    logger = logging.getLogger(name=name)
    logger.setLevel(level.upper())

    # Configure console output
    sh = logging.StreamHandler()
    sh_formatter = logging.Formatter("%(message)s")
    sh.setFormatter(sh_formatter)
    sh.setLevel(console_level.upper())
    logger.addHandler(sh)

    # Configure file output
    if filename is not None:
        fh = logging.FileHandler(filename=filename, mode="a")
        fh_formatter = logging.Formatter(fmt=fmt, datefmt=datefmt)
        fh.setFormatter(fh_formatter)
        fh.setLevel(file_level.upper())
        logger.addHandler(fh)

    return logger
