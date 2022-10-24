#!/user/bin/env python
"""
Fetch Catalog

Fetches all pages from the search results in Repositorio.

Output: "output/{target}/{url_number}.html"
"""
import asyncio
import random
import sqlite3
import time
import traceback
from requests_html import HTMLSession
from contextlib import closing
from pathlib import Path

import aiohttp

from src.logging import configLogger

log = configLogger(snakemake.rule, snakemake.config["logfile"])  # type: ignore # noqa


async def fetch(session, url, output_folder):
    """Fetch HTML from url and write [paper_id].html"""
    async with session.get(url[1]) as response:
        try:
            html_doc = await response.text()
        except Exception:
            log.exception(f"UnicodeDecodeError in fetch(): url: {url[1]}, file: {url[0]}.html.")
            html_doc = traceback.format_exc()
            with open(f"{output_folder}/{url[2]}/errors.txt", "a") as f:
                f.write(f"{url[0]}.html\n")
    with open(f"{output_folder}/{url[2]}/{url[0]}.html", "w") as f:
        log.info(f"Writing {output_folder}/{url[2]}/{url[0]}.html.")
        f.write(html_doc)


async def fetch_with_sem(session, url, output_folder, sem, rate):
    """Fetch HTML from url with Semaphore"""
    async with sem:
        t1 = time.perf_counter()
        sleep_task = asyncio.create_task(asyncio.sleep(rate + random.random()))
        await fetch(session, url, output_folder)
        await sleep_task
        t2 = time.perf_counter()
        log.debug(f"Fetched in {t2-t1:.0f} seconds.")


async def fetch_all(urls, output_folder, limit, rate):
    """ "Fetch list of HTMLs from list of urls"""
    sem = asyncio.Semaphore(limit)
    async with aiohttp.ClientSession() as session:
        tasks = (fetch_with_sem(session, url, output_folder, sem, rate) for url in urls)
        result = await asyncio.gather(*tasks)
        return result


def main():
    rule = snakemake.rule  # type: ignore # noqa
    date = snakemake.config["date"]  # type: ignore # noqa
    limit = snakemake.config["limit"]  # type: ignore # noqa
    rate = snakemake.config["rate"]  # type: ignore # noqa
    # institution = snakemake.wildcards.institution  # type: ignore # noqa
    institution_list = snakemake.config["institutions"]  # type: ignore # noqa
    target = rule.split("fetch_")[1]

    if rule == "fetch_catalog":
        n_catalog_urls = snakemake.config["n_catalog_urls"]  # type: ignore # noqa
        url_template = "https://repositorio.usp.br/result.php?filter%5B0%5D=unidadeUSP%3A%22{}%22&fields%5B0%5D=name&fields%5B1%5D=author.person.name&fields%5B2%5D=authorUSP.name&fields%5B3%5D=about&fields%5B4%5D=description&fields%5B5%5D=unidadeUSP&page={}"

        session = HTMLSession()

        r = session.get(url_template.format(institution, 1))
        print(r.html)
        raise Exception('My exception!')
        

        urls = []
        for institution in institution_list:
            urls_institution = [
                (i, url_template.format(institution, i + 1), institution) for i in range(n_catalog_urls)
            ]
            urls = urls + urls_institution

    elif rule == "fetch_papers":
        n_papers = snakemake.config.get("n_papers")  # type: ignore # noqa
        url_root = "https://repositorio.usp.br/"

        # Fill paper_id columns
        counts = []
        for institution in institution_list:
            with closing(sqlite3.connect(f"saopaulo_{institution}.db")) as conn:
                with conn:
                    count_cursor = conn.execute("SELECT COUNT(*) FROM papers where institution=(?)", (institution,))
                    count = count_cursor.fetchall()[0][0]
            log.info(f"Count of papers in {institution}: {count}")
            counts.append(count)
            log.info(f"Counts: {counts}")
        paper_ids = [i for count in counts for i in range(count)]
        ids = [i + 1 for i in range(sum(counts))]
        paper_ids = list(zip(paper_ids, ids))
        log.info("length")
        log.info(len(paper_ids))
        log.info(len(ids))
        log.info(f"paper_ids: {paper_ids}")

        log.info(f"Before filling paper_ids:")
        with closing(sqlite3.connect(f"saopaulo_{institution}.db")) as conn:
            with conn:
                result = conn.execute("SELECT * FROM papers LIMIT 10;")
                log.info(result.fetchall())
        with closing(sqlite3.connect(f"saopaulo_{institution}.db")) as conn:
            with conn:
                conn.executemany(
                    "UPDATE papers SET paper_id=(?) WHERE id=(?)",
                    paper_ids,
                )
        log.info(f"After filling paper_ids:")
        with closing(sqlite3.connect(f"saopaulo_{institution}.db")) as conn:
            with conn:
                result = conn.execute("SELECT * FROM papers LIMIT 105;")
                log.info(result.fetchall())

        # Get list of urls
        urls = []
        for institution in institution_list:
            with closing(sqlite3.connect(f"saopaulo_{institution}.db")) as conn:
                with conn:
                    # count_cursor = conn.execute("SELECT COUNT(*) FROM papers where institution=(?)", (institution,))
                    # count = count_cursor.fetchall()[0][0]
                    # log.info(f"Count of papers in {institution}: {count}")
                    # paper_id_institution = [(i, institution) for i in range(count)]
                    # conn.executemany(
                    #     "UPDATE papers SET paper_id=(?) WHERE institution=(?) AND paper_id IS NULL",
                    #     paper_id_institution,
                    # )
                    result = conn.execute(
                        "SELECT paper_id, url_stem FROM papers WHERE institution=(?) LIMIT (?)",
                        (institution, n_papers),
                    )
                    rows = result.fetchall()

            urls_institution = [(row[0], url_root + row[1], institution) for row in rows]
            urls = urls + urls_institution
            log.info(f"Urls for institution {institution}: {urls[:3]}")

    output_folder = f"output/{target}"
    Path(output_folder).mkdir(parents=True, exist_ok=True)

    log.info(f"Fetching {len(urls)} htmls.")
    t1 = time.perf_counter()
    asyncio.run(fetch_all(urls, output_folder, limit, rate))
    t2 = time.perf_counter()
    log.info(f"Done in {t2-t1:.0f} seconds.")


if __name__ == "__main__":
    main()
