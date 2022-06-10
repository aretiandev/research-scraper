#!/user/bin/env python
"""
Parse HTML

Parse HTML from catalog or papers and extract relevant data.
- From catalog: extract paper urls and write to 'urls' table.
- From papers: extract paper information and write to 'papers' table.

Input: "output/{target}/{url_number}.html"
Output: rows in "saopaulo.db"
"""
from pathlib import Path

from bs4 import BeautifulSoup

from src.logging import configLogger
from src.sqlite import insert_edges, insert_paper, insert_urls

log = configLogger(snakemake.rule, snakemake.config["logfile"])  # type: ignore # noqa


def parse_catalog_html(html_doc):
    """Parse catalog HTML and extract paper urls"""
    soup = BeautifulSoup(html_doc, "html5lib")
    articles = soup.find_all("article")
    paper_urls = [article.p.a["href"] for article in articles]
    return set(paper_urls)


def parse_paper_html(html_doc):
    """Parse paper HTML and extract paper data"""
    soup = BeautifulSoup(html_doc, "html5lib")

    list_elems = soup.select('article article ul li:-soup-contains("Authors: ") ul li')

    orcids = []
    authors = []
    author_urls = []

    for list_elem in list_elems:
        a_tags = list_elem.find_all("a")

        author_name = a_tags[0].text
        authors.append(author_name)

        author_url_raw = a_tags[0]["href"][2:]
        author_url = author_url_raw.split("/")[1].replace(" ", "%20").replace('"', "%22")
        author_urls.append(author_url)

        try:
            orcid = a_tags[1]["href"].split("/")[-1]
        except IndexError:
            orcid = ""
        orcids.append(orcid)

    subject_content = soup.select('article article ul li:-soup-contains("Subject")')
    try:
        subjects = [a_tag.text for a_tag in subject_content[0].find_all("a")]
    except IndexError:
        log.exception(f"Error parsing subjects from paper HTML. File: {snakemake.input[0]}")  # type: ignore # noqa
        subjects = []

    paper = {
        "authors": authors,
        "orcids": orcids,
        "author_urls": author_urls,
        "subjects": subjects,
    }

    return paper


def main():
    rule = snakemake.rule  # type: ignore # noqa
    input_file = snakemake.input[0]  # type: ignore # noqa
    date = snakemake.config["date"]  # type: ignore # noqa
    institution = snakemake.wildcards.institution  # type: ignore # noqa

    with open(input_file) as f:
        html_doc = f.read()

    if rule == "parse_catalog":
        paper_urls = parse_catalog_html(html_doc)
        # log.info(f"Parsed catalog HTML: {paper_urls}")
        log.info(f"Parsed catalog HTML from {institution}: {input_file}")
        insert_urls(paper_urls, institution, date)
    elif rule == "parse_papers":
        paper = parse_paper_html(html_doc)
        paper["paper_id"] = Path(input_file).stem
        paper["institution"] = institution
        insert_paper(paper)
        if len(paper["orcids"]):
            insert_edges(paper["orcids"], institution, date)


if __name__ == "__main__":
    main()
