import sqlite3
import sys
from contextlib import closing

import yaml

from src.logging import configLogger

log = configLogger("create_tables")  # type: ignore # noqa


def main(db):
    with closing(sqlite3.connect(db)) as conn:
        with conn:
            # conn.execute(
            #     """CREATE TABLE IF NOT EXISTS urls (
            #         id INTEGER PRIMARY KEY,
            #         items TEXT,
            #         url_stem TEXT UNIQUE,
            #         date_created TEXT,
            #         current INTEGER,
            #         url_scraped INTEGER
            #         )"""
            # )

            conn.execute(
                # """
                # CREATE TABLE IF NOT EXISTS papers (
                #     id INTEGER PRIMARY KEY,
                #     authors TEXT,
                #     orcids TEXT,
                #     author_urls TEXT,
                #     subjects TEXT,
                #     institution TEXT,
                #     date_created TEXT
                # )
                # """
                """
                CREATE TABLE IF NOT EXISTS papers(
                    id INTEGER PRIMARY KEY,
                    paper_id INTEGER,
                    institution TEXT,
                    url_stem TEXT UNIQUE,
                    date_created TEXT,
                    current INT,
                    url_scraped INT,
                    authors TEXT,
                    orcids TEXT,
                    author_urls TEXT,
                    subjects TEXT
                    );
                """
            )

            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS edges (
                    id INTEGER PRIMARY KEY,
                    Source TEXT,
                    Target TEXT,
                    institution TEXT,
                    date_created TEXT,
                    UNIQUE(Source,Target)
                )
                """
            )

    log.info(f"Created {db} and tables: 'papers', 'edges'.")


if __name__ == "__main__":
    # Get institution from config.yml
    with open('config.yml') as f:
        config = yaml.safe_load(f)
        try:
            institution = config.get('institutions')[0]
            db = f"saopaulo_{institution}.db"
        except TypeError:
            pass

    # Override if supplied in command line
    if len(sys.argv) > 1:
        db = sys.argv[1]

    # If not set, use saopaulo.db as default
    try: 
        db
    except NameError:
        db = "saopaulo.db"

    main(db)
