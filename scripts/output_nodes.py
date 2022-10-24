#!/user/bin/env python
"""
Output Nodes

Output: "output/gephi/nodes_{institution}.csv"
"""
import itertools
import sqlite3
import time
from contextlib import closing

import pandas as pd

from src.logging import configLogger

log = configLogger(snakemake.rule, snakemake.config["logfile"])  # type: ignore # noqa


def main():
    institution = snakemake.wildcards.institution  # type: ignore # noqa
    db = f"saopaulo_{institution}.db"
    out_file = snakemake.output[0]  # type: ignore # noqa

    with closing(sqlite3.connect(db)) as conn:
        log.info("Querying database.")
        with conn:
            result = conn.execute(
                """
                SELECT orcids, authors, subjects FROM papers 
                WHERE institution=?
                AND orcids IS NOT NULL
                AND authors IS NOT NULL
                AND subjects IS NOT NULL;
                """,
                (institution,),
            )

            rows = result.fetchall()

    # Create orcids set
    log.info("Creating set of orcids.")
    try:
        orcids_g = (eval(row[0]) for row in rows)
    except Exception:
        log.exception("Error.")
        import pdb

        pdb.set_trace()
    orcids = set(itertools.chain.from_iterable(orcids_g))
    orcids.remove("")

    new_rows = []
    for row in rows:
        lst = [eval(lst) for lst in row]
        new_rows.append(lst)

    t1 = time.perf_counter()
    # Attach subject to each orcid key
    log.info("Getting names and subjects.")
    nodes = {}
    for orcid in orcids:
        for row in new_rows:
            if orcid in row[0]:
                # Get author name
                idx = row[0].index(orcid)
                name_lst = row[1]
                nodes[orcid] = [name_lst[idx]]

                # Get subjects
                nodes[orcid].append(row[2][0])

                break
    t2 = time.perf_counter()
    log.info(f"Done in {t2-t1:.0f} seconds.")

    # Creating dataframe
    df = pd.DataFrame.from_dict(nodes, orient="index").reset_index()
    df.columns = ["id", "label", "subject"]
    df.to_csv(out_file, index=False)
    log.info(f"Saved {out_file}.")


if __name__ == "__main__":
    main()
