"""Test Scraper

Run tests with:
    python -m unittest

Create pickle data with:
    python -m tests pickle

For coverage report run:
    coverage run -m unittest
    coverage report
"""
import unittest
from unittest.mock import patch
import re
from requests_html import HTML
from scrape import build_urls
from src.scrape import scrape, get_urls, WebsiteDownError, insert_urls
import logging
import sqlite3

logging.getLogger("asyncio").setLevel("ERROR")

id_regex = {}
id_regex["author_urls"] = re.compile(
    r"^orcid/[0-9X]{4}-[0-9X]{4}-[0-9X]{4}-[0-9X]{4}\b"
)
id_regex["paper_urls"] = re.compile(r"^([0-9]{8}|book/isbn/[0-9-].*|article/doi/.*)\b")
id_regex["project_urls"] = re.compile(r"^cris/project/pj[0-9]{7}$")
id_regex["group_urls"] = re.compile(r"^(cris|sgr)/[a-zA-Z0-9]+.*$")
url_regex = re.compile(
    r"^https://portalrecerca.csuc.cat/simple-search\?query=.*&start=[0-9]*"
)


class TestBuildURLs(unittest.TestCase):
    """
    Test author_urls, paper_urls, group_urls and project_urls
    are built properly before scraping
    """

    n_urls = 30
    items_list = ["author_urls", "paper_urls", "group_urls", "project_urls"]

    def setUp(self):
        """Setup mocks for functions that trigger HTTP requests"""
        patcher = patch("src.scrape.get_max_pages")
        self.mock_get_max_pages = patcher.start()
        self.mock_get_max_pages.return_value = self.n_urls
        self.addCleanup(patcher.stop)

    def test_get_urls(self):
        """Check function get_urls() returns URLs with correct URL format"""
        for items in self.items_list:
            with self.subTest(items=items):
                urls = get_urls(items=items)
                urls_match_regex = all([url_regex.match(url) for url in urls])
                self.assertTrue(urls_match_regex)

    def test_get_urls_length(self):
        """Check function get_urls returns a list with correct length"""
        for items in self.items_list:
            with self.subTest(items=items):
                urls = get_urls(items=items)
                self.assertIsInstance(urls, list)
                self.assertEqual(len(urls), self.n_urls)

    def test_get_urls_websitedown(self):
        """Check if website is down, function get_urls() raises WebsiteDownError"""
        self.mock_get_max_pages.side_effect = WebsiteDownError("url")
        for items in self.items_list:
            with self.subTest(items=items):
                with self.assertRaises(WebsiteDownError):
                    get_urls(items=items)

    def test_build_urls(self):
        """Check function build_urls() returns URLs with correct URL format"""
        for items in self.items_list:
            with self.subTest(items=items):
                urls = build_urls(items=items)
                urls_match_regex = all([url_regex.match(url) for url in urls])
                self.assertTrue(urls_match_regex)

    def test_build_urls_length(self):
        """Check function build_urls() returns a list of URLs with correct length"""
        for items in self.items_list:
            with self.subTest(items=items):
                urls = build_urls(items=items)
                self.assertIsInstance(urls, list)
                self.assertEqual(len(urls), self.n_urls)


def get_html_object(items="author_urls"):
    """Helper function to build HTML object from pre-saved html page"""
    with open(f"tests/{items}.html", "rb") as f:
        raw_html = f.read()
    with open(f"tests/{items}.txt", "r") as f:
        url = f.read()
    return HTML(url=url, html=raw_html)


class TestScrape(unittest.IsolatedAsyncioTestCase):
    """Test data is scraped properly"""

    items_list = ["author_urls", "paper_urls", "project_urls", "group_urls"]

    @patch("src.scrape.retry_url")
    async def test_scrape_author_urls(self, mock_retry_url):
        """Check scrape(author_urls) returns a list of URLs for the authors pages"""
        for items in self.items_list:
            with self.subTest(items=items):
                html = get_html_object(items)
                mock_retry_url.return_value.html = html
                result = await scrape(items=items, urls=[[html.url]])

                self.assertEqual(len(result), 300)

                if items == "author_urls":
                    urls_match_regex = all(
                        [id_regex[items].match(url) for url in result]
                    )
                    self.assertTrue(urls_match_regex)


class TestSQL(unittest.IsolatedAsyncioTestCase):
    """Test SQL pipeline: reading and writing to sqlite database"""

    items_list = ["author_urls", "paper_urls", "project_urls", "group_urls"]

    @patch("src.scrape.retry_url")
    async def test_sql_urls(self, mock_retry_url):
        """Test urls are written to the database"""
        items = "author_urls"
        date_today = "20220427"
        db = "../recerca_test.db"

        for items in self.items_list:
            with self.subTest(items=items):

                # Write to SQLite
                html = get_html_object(items)
                mock_retry_url.return_value.html = html
                result = await scrape(items=items, urls=[[html.url]])
                insert_urls(result, items, date_today, db=db)

                # Read from SQLite
                conn = sqlite3.connect(db)
                with conn:
                    for row in conn.execute(
                        "SELECT * FROM urls WHERE items=?;", (items,)
                    ):
                        self.assertIsInstance(row[0], int)
                        self.assertEqual(row[1], items)
                        try:
                            self.assertTrue(id_regex[items].match(row[2]))
                        except AssertionError:
                            print()
                            print(id_regex[items])
                            print(row[2])
                            print()
                            raise
                        self.assertEqual(row[3], date_today)
                        self.assertIsInstance(row[4], int)
                        self.assertIsInstance(row[5], int)

                    cur_count = conn.execute(
                        "SELECT COUNT(*) FROM urls WHERE items=?;", (items,)
                    )
                    self.assertEqual(cur_count.fetchall()[0][0], 300)

                    cur_date = conn.execute(
                        "SELECT DISTINCT date_created FROM urls WHERE items=?;",
                        (items,),
                    )
                    self.assertEqual(cur_date.fetchall()[0][0], date_today)
                conn.close()


if __name__ == "__main__":
    unittest.main()
