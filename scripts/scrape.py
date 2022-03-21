import asyncio
import pandas as pd
from src.scrape import scrape, WebsiteDownError
from src.logging import create_logger
from src.bot import ping_and_wait

log = create_logger(__name__, f"log/{__name__}.log")


def main():
    url_root = 'https://portalrecerca.csuc.cat'

    items = snakemake.wildcards.item_name + '_links'
    batch_size = snakemake.params.get('batch_size')
    if batch_size is None:
        batch_size = 20
    out_file = snakemake.output[0]

    if items in ['author_data', 'paper_data', 'group_data', 'project_data']:
        item_urls = pd.read_csv(snakemake.input[0])
        item_urls = list(item_urls['0'])
        if items == 'paper_data':
            urls = [url_root + url + '?mode=full' for url in item_urls]
        else:
            urls = [url_root + url for url in item_urls]
    else:
        urls = None

    if items in ['author_links', 'paper_links', 'group_links', 'project_links']:
        if items == 'paper_links':
            batch_size = 50

    while True:
        try:
            asyncio.run(scrape(items=items, urls=urls, batch_size=batch_size, out_file=out_file))
            break
        except WebsiteDownError:
            log.info(f"Website is down during execution of scrape({items}). Pinging URL...")
            ping_and_wait(url_root)


if __name__ == "__main__":
    main()
