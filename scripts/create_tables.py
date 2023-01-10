import sqlite3
import sys
from contextlib import closing

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

    log.info(f"Created '{db}' and tables: 'papers', 'edges'.")


if __name__ == "__main__":
    db = "recerca.db"
    if len(sys.argv) > 1:
        db = sys.argv[1]
    main(db)
