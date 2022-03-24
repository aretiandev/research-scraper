import logging


def create_logger(
        name,
        logfile='log/snakemake.log',
        level=None,
        file_level=None,
        stream_level=None,
        fmt='%(asctime)s %(name)s %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %I:%M:%S%p'):

    logger = logging.getLogger(name)

    level_dict = {
        'debug': logging.DEBUG,
        'info': logging.INFO,
        'warning': logging.WARNING,
        'error': logging.ERROR,
        'critical': logging.CRITICAL,
    }

    logging_level = logging.INFO
    if level:
        logging_level = level_dict[level]
    logger.setLevel(logging_level)

    fh = logging.FileHandler(logfile, mode="w")
    fh_formatter = logging.Formatter(fmt=fmt, datefmt=datefmt)
    fh.setFormatter(fh_formatter)
    if file_level:
        file_logging_level = level_dict[file_level]
        fh.setLevel(file_logging_level)

    sh = logging.StreamHandler()
    sh_formatter = logging.Formatter('%(message)s')
    sh.setFormatter(sh_formatter)
    if stream_level:
        stream_logging_level = level_dict[stream_level]
        sh.setLevel(stream_logging_level)

    logger.addHandler(fh)
    logger.addHandler(sh)

    return logger
