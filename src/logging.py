import logging


def create_logger(
        name,
        logfile='log/snakemake.log',
        level=logging.INFO,
        fmt='%(asctime)s %(name)s %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %I:%M:%S%p'):

    logger = logging.getLogger(name)
    logger.setLevel(level)

    fh = logging.FileHandler(logfile, mode="w")
    fh_formatter = logging.Formatter(fmt=fmt, datefmt=datefmt)
    fh.setFormatter(fh_formatter)
    # fh.setLevel(logging.INFO)

    sh = logging.StreamHandler()
    sh_formatter = logging.Formatter('%(message)s')
    sh.setFormatter(sh_formatter)
    # sh.setLevel(logging.WARNING)

    logger.addHandler(fh)
    logger.addHandler(sh)

    return logger
