import pandas as pd
from datetime import date
from src.scrape import scrape
import asyncio

date_today = date.today().strftime("%Y%m%d")

### SCRAPE AUTHOR_LINKS

async def test_scrape_author_links():
    author_urls = pd.read_csv(f'data/author_urls_{date_today}.csv')
    author_urls = list(author_urls['0'])
    url_root = 'https://portalrecerca.csuc.cat'
    urls = [url_root + url for url in author_urls]

    items = 'authors'
    batch_size = 100
    out_file = f'data/nodes_{date_today}.csv'
    return await scrape(items=items, urls=urls, batch_size=batch_size, out_file=out_file)


### SCRAPE AUTHORS

async def test_scrape_authors():
    # Build urls
    author_urls = pd.read_csv(f'./data/author_urls_{date_today}.csv')
    author_urls = list(author_urls['0'])
    url_root = 'https://portalrecerca.csuc.cat'
    urls = [url_root + url for url in author_urls]

    # Get author data in batch
    items = 'authors'
    batch_size = 100
    out_file = f'./data/nodes_{date_today}.csv'
    return await scrape(items=items, urls=urls, batch_size=batch_size, out_file=out_file)


if __name__ == "__main__":
    asyncio.run(test_scrape_authors())