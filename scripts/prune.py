import sqlite3
import subprocess
import sys


def main(db):
    subprocess.run(["rm", "-rf", "output"])
    subprocess.run(["rm", "-rf", "log"])
    subprocess.run(["rm", "saopaulo.db"])

    conn = sqlite3.connect(db)
    c = conn.cursor()

    c.execute(
        """CREATE TABLE IF NOT EXISTS urls (
            id INTEGER PRIMARY KEY,
            items TEXT,
            url_stem TEXT UNIQUE,
            date_created TEXT,
            current INTEGER,
            url_scraped INTEGER
            )"""
    )

    conn.commit()
    conn.close()


if __name__ == "__main__":
    db = "saopaulo.db"
    if len(sys.argv) > 1:
        db = sys.argv[1]
    main(db)
