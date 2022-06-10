#!/user/bin/env python
"""
Output Edges

Output: "output/gephi/edges_{institution}.csv"
"""
import sqlite3
from contextlib import closing

import pandas as pd

from src.logging import configLogger

log = configLogger(snakemake.rule, snakemake.config["logfile"])  # type: ignore # noqa


def main(db="saopaulo.db"):
    out_file = snakemake.output[0]  # type: ignore # noqa

    with closing(sqlite3.connect(db)) as conn:
        with conn:
            try:
                result = conn.execute(
                    """
                    SELECT * FROM edges LIMIT ?;
                    """,
                    (snakemake.config["n_edges"],),  # type: ignore # noqa
                )
            except Exception:
                log.exception("Error.")
                import pdb

                pdb.set_trace()
            rows = result.fetchall()

    df = pd.DataFrame(rows)
    df = df[[1, 2]]
    df.columns = ["Source", "Target"]

    df.to_csv(out_file, index=False)
    log.info(f"Saved {out_file}.")


if __name__ == "__main__":
    main()
