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
from contextlib import closing
from pathlib import Path

import aiohttp

from src.logging import configLogger

log = configLogger(snakemake.rule, snakemake.config["logfile"])  # type: ignore # noqa


async def fetch(session, url, output_folder, output_file):
    """Fetch HTML from url"""
    async with session.get(url) as response:
        try:
            html_doc = await response.text()
        except Exception as e:
            html_doc = e
            log.exception(f"Exception in fetch(): url: {url}, file: {output_file}.")
            with open("{output_folder}/errors.txt", "a") as f:
                f.write(f"{output_file}.txt")
        with open(f"{output_folder}/{output_file}.html", "w") as f:
            f.write(html_doc)


async def fetch_with_sem(session, url, output_folder, output_file, sem, rate):
    """Fetch HTML from url with Semaphore"""
    async with sem:
        t1 = time.perf_counter()
        sleep_task = asyncio.create_task(asyncio.sleep(rate + random.random()))
        await fetch(session, url, output_folder, output_file)
        await sleep_task
        t2 = time.perf_counter()
        log.debug(f"Fetched in {t2-t1:.0f} seconds.")


async def fetch_all(urls, output_folder, limit, rate):
    """ "Fetch list of HTMLs from list of urls"""
    sem = asyncio.Semaphore(limit)
    async with aiohttp.ClientSession() as session:
        tasks = (fetch_with_sem(session, url, output_folder, i, sem, rate) for i, url in enumerate(urls))
        result = await asyncio.gather(*tasks)
        return result


def main():
    rule = snakemake.rule  # type: ignore # noqa
    limit = snakemake.config["limit"]  # type: ignore # noqa
    rate = snakemake.config["rate"]  # type: ignore # noqa

    if rule == "fetch_catalog":
        n_catalog_urls = snakemake.config["n_catalog_urls"]  # type: ignore # noqa
        url_root = "https://repositorio.usp.br/result.php?filter%5B0%5D=unidadeUSP%3A%22EP%22&fields%5B0%5D=name&fields%5B1%5D=author.person.name&fields%5B2%5D=authorUSP.name&fields%5B3%5D=about&fields%5B4%5D=description&fields%5B5%5D=unidadeUSP&page="
        urls = [url_root + str(i) for i in range(1, n_catalog_urls + 1)]
        output_folder = "output/catalog"
    elif rule == "fetch_papers":
        n_papers = snakemake.config.get("n_papers")  # type: ignore # noqa
        with closing(sqlite3.connect("saopaulo.db")) as conn:
            with conn:
                if n_papers:
                    result = conn.execute("SELECT * FROM urls LIMIT (?)", (n_papers,))
                else:
                    result = conn.execute("SELECT * FROM urls")
                rows = result.fetchall()
        url_root = "https://repositorio.usp.br/"
        urls = [url_root + row[2] for row in rows]
        output_folder = "output/papers"

    Path(output_folder).mkdir(exist_ok=True)

    log.info(f"Fetching {len(urls)} htmls.")
    t1 = time.perf_counter()
    asyncio.run(fetch_all(urls, output_folder, limit, rate))
    t2 = time.perf_counter()
    log.info(f"Done in {t2-t1:.0f} seconds.")


if __name__ == "__main__":
    main()
