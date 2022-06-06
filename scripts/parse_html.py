#!/user/bin/env python
"""
Get Paper URLs

Parse HTML from catalog and extract links to papers.

Input: "output/{target}/{url_number}.html"
Output: rows in "saopaulo.db"
"""
from bs4 import BeautifulSoup

from src.logging import configLogger
from src.sqlite import insert_urls

log = configLogger(snakemake.rule, snakemake.config['logfile'])  # type: ignore # noqa


def get_paper_urls(catalog_html):
    """Parse list of HTMLs and extract paper urls"""
    soup = BeautifulSoup(catalog_html, "html5lib")
    articles = soup.find_all("article")
    paper_urls = [article.p.a["href"] for article in articles]
    return set(paper_urls)


def main():
    input_file = snakemake.input[0]  # type: ignore # noqa
    date = snakemake.config["date"]  # type: ignore # noqa

    with open(input_file) as f:
        catalog_html = f.read()

    paper_urls = get_paper_urls(catalog_html)
    insert_urls(urls=paper_urls, items="paper_urls", date=date)


if __name__ == "__main__":
    main()
