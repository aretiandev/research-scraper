import sqlite3
import subprocess
import sys

from src.logging import configLogger

log = configLogger(snakemake.rule, snakemake.config["logfile"])  # type: ignore # noqa


def main(db):
    subprocess.run(["rm", "-rf", "output"])
    log.info("Deleted 'output/'.")
    subprocess.run(["rm", "-rf", "log"])
    log.info("Deleted 'log/'.")
    subprocess.run(["rm", "saopaulo.db"])
    log.info("Deleted 'saopaulo.db'.")


if __name__ == "__main__":
    db = "saopaulo.db"
    if len(sys.argv) > 1:
        db = sys.argv[1]
    main(db)
