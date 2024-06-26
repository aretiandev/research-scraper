import asyncio
import logging
import os
import sqlite3

import pandas as pd
from src.bot import ping_and_wait
from src.logging import configLogger
from src.scrape import WebsiteDownError, get_urls, scrape

log = logging.getLogger(__name__)

url_root = "https://portalrecerca.csuc.cat/"


def build_urls(items, institution):
    """Build URLs for urls and data rules"""

    # Build URLs for url rules
    if items in ["author_urls", "paper_urls", "group_urls", "project_urls"]:
        while True:
            try:
                urls = get_urls(items=items, institution=institution)
                break
            except WebsiteDownError:
                message = (
                    ":warning: Website is down."
                    + "Pinging URL and resuming when back online."
                )
                ping_and_wait(url_root, slack_msg=message, notifications=1)

    # Build URls for data rules
    elif items in ["author_data", "paper_data", "group_data", "project_data"]:
        try:
            item_urls = pd.read_csv(snakemake.input[0])  # type: ignore # noqa
        except pd.errors.EmptyDataError:
            return []
        item_urls = list(set(item_urls["0"]))
        item_urls.sort()
        if items == "paper_data":
            urls = [url_root + url + "?mode=full" for url in item_urls]
        else:
            urls = [url_root + url for url in item_urls]

    return urls


def main():
    # Configure root logger. All logger inherit this configuration.
    configLogger(filename="log/snakemake.log", file_level="debug", console_level="info")

    # Parameters
    batch_size = snakemake.params.get("batch_size")  # type: ignore # noqa
    timeout = snakemake.params.get("timeout")  # type: ignore # noqa
    database = snakemake.params.get("database")  # type: ignore # noqa
    out_file = snakemake.output[0]  # type: ignore # noqa
    items = out_file.split("/")[-1].split(".")[0][9:].rsplit("_", 1)[0]
    # institution_list = snakemake.params.get("institution_list")  # type: ignore # noqa
    institution = snakemake.wildcards.institution  # type: ignore # noqa

    urls = build_urls(items=items, institution=institution)
    if not urls:
        result_df = pd.DataFrame()
        result_df.to_csv(out_file, index=None)
        log.debug(f"No {items} found. Saved empty dataframe: {out_file}")
        return

    # Batch size and start_pos
    if batch_size is None:
        batch_size = 20
        if items == "paper_urls":
            batch_size = 50
    start_pos = 0

    # Start scraper
    while True:
        try:
            log.info("Building list of URLs to scrape.")

            # Create batch URLs
            batch_urls = [
                urls[i : i + batch_size]
                for i in range(start_pos, len(urls), batch_size)
            ]

            asyncio.run(
                scrape(
                    items=items,
                    urls=batch_urls,
                    out_file=out_file,
                    out_sql=False,
                    timeout=timeout,
                    database=database,
                )
            )
            break

        except WebsiteDownError as e:
            last_pos = urls.index(e.url)
            log.info(
                f"Website down during scrape({items}), url #{last_pos}, url: {e.url}"
            )

            ping_and_wait(url_root, notifications=1)

            # Move existing results to backup folder before continuing
            if os.path.isfile(out_file):
                new_folder = "data/archive"
                log.debug(f"Existing file found! Moving to folder: {new_folder}")
                if not os.path.exists(new_folder):
                    log.debug(f"Creating new folder: {new_folder}")
                    os.makedirs(new_folder)
                out_file_name = out_file.split("/")[-1]
                new_out_file = new_folder + "/" + out_file_name
                os.rename(out_file, new_out_file)
                log.info(f"Moved existing results to: {new_out_file}")

            # Continue where we left off
            start_pos = last_pos

            log.info(f"Restarting scrape({items}) from position {start_pos}.")

    # SQL: reset current variable for old data
    # if items in ["author_urls", "paper_urls", "group_urls", "project_urls"]:
    #     if out_file is not None:
    #         date_today = out_file.split("/")[-1][:8]
    #         conn = sqlite3.connect(database)
    #         with conn:
    #             conn.execute(
    #                 "UPDATE urls SET current = 0 WHERE date_created != ?", (date_today,)
    #             )
    #         conn.close()


if __name__ == "__main__":
    main()
