import pandas as pd
from datetime import date
from src.helpers import scrape
import asyncio

date_today = date.today().strftime("%Y%m%d")
url_root = 'https://portalrecerca.csuc.cat'

rule author_links:
    output:
        f'data/author_urls_{date_today}.csv'
    run:
        items = 'author_links'
        batch_size = 2
        asyncio.run(scrape(items=items, batch_size=batch_size, out_file={output}))

rule authors:
    input:
        f'data/author_urls_{date_today}.csv'
    output:
        f'data/nodes_{date_today}.csv'
    run:
        author_urls = pd.read_csv({input})
        author_urls = list(author_urls['0'])
        urls = [url_root + url for url in author_urls]
        items = 'authors'
        batch_size = 100
        asyncio.run(scrape(items=items, urls=urls, batch_size=batch_size, out_file={output}))

rule clean_authors:
    input:
        f'data/nodes_{date_today}.csv'
    shell:
        "rm {input}"

rule project_links:
    output:
        f'data/project_links_{date_today}.csv'
    run:
        items = 'project_links'
        batch_size = 10
        asyncio.run(scrape(items=items, batch_size=batch_size, out_file={output}))

rule projects:
    input:
        f'data/project_links_{date_today}.csv'
    output:
        f'data/projects_{date_today}.csv'
    run:
        project_urls = pd.read_csv({input})
        project_urls = list(project_urls['0'])
        urls = [url_root + url for url in project_urls]
        items = 'projects'
        batch_size = 10
        asyncio.run(scrape(items=items, urls=urls, batch_size=batch_size, out_file={output}))

rule group_links:
    output:
        f'data/group_links_{date_today}.csv'
    run:
        items = 'group_links'
        batch_size = 10
        asyncio.run(scrape(items=items, batch_size=batch_size, out_file={output}))

rule groups:
    input:
        f'data/group_links_{date_today}.csv'
    output:
        f'data/groups_{date_today}.csv'
    run:
        group_urls = pd.read_csv({input})
        group_urls = list(group_urls['0'])
        urls = [url_root + url for url in group_urls]
        items = 'groups'
        batch_size = 10
        asyncio.run(scrape(items=items, urls=urls, batch_size=batch_size, out_file={output}))

rule paper_links:
    output:
        f'data/paper_links_{date_today}.csv'
    run:
        items = 'paper_links'
        batch_size = 10
        asyncio.run(scrape(items=items, batch_size=batch_size, out_file={output}))

rule papers:
    input:
        f'data/paper_links_{date_today}.csv'
    output:
        f'data/papers_{date_today}.csv'
    run:
        paper_urls = pd.read_csv({input})
        paper_urls = list(paper_urls['0'])
        urls = [url_root + url + '?mode=full' for url in paper_urls]
        items = 'papers'
        batch_size = 10
        asyncio.run(scrape(items=items, urls=urls, batch_size=batch_size, out_file={output}))
