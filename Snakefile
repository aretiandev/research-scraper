import pandas as pd
from datetime import date
from src.helpers import scrape
import asyncio

date_today = date.today().strftime("%Y%m%d")

rule author_links:
    output:
        f'data/author_urls_{date_today}.csv'
    run:
        items = 'author_links'
        batch_size = 2
        out_file = f'data/author_urls_{date_today}.csv'
        asyncio.run(scrape(items=items, batch_size=batch_size, out_file=out_file))

rule authors:
    input:
        f'data/author_urls_{date_today}.csv'
    output:
        f'data/nodes_{date_today}_0.csv'
    run:
        author_urls = pd.read_csv(f'data/author_urls_{date_today}.csv')
        author_urls = list(author_urls['0'])
        url_root = 'https://portalrecerca.csuc.cat'
        urls = [url_root + url for url in author_urls]

        items = 'authors'
        batch_size = 100
        out_file = f'data/nodes_{date_today}_0.csv'
        asyncio.run(scrape(items=items, urls=urls, batch_size=batch_size, out_file=out_file))

rule clean_authors:
    input:
        f'data/nodes_{date_today}.csv'
    shell:
        "rm {input}"
