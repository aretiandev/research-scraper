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
from src.scrape import scrape, get_urls, WebsiteDownError
import logging

logging.getLogger("asyncio").setLevel("ERROR")


class TestURLs(unittest.TestCase):
    """Test author_urls, paper_urls, group_urls and project_urls"""

    url_regex = re.compile(
        r"^https://portalrecerca.csuc.cat/simple-search\?query=.*&start=[0-9]*"
    )
    n_urls = 30
    items_list = ["author_urls", "paper_urls", "group_urls", "project_urls"]

    def setUp(self):
        """Setup mocks for functions that trigger HTTP requests"""
        patcher = patch("src.scrape.get_max_pages")
        self.mock_get_max_pages = patcher.start()
        self.mock_get_max_pages.return_value = self.n_urls
        self.addCleanup(patcher.stop)

    def test_get_urls(self):
        """Function get_urls() returns URLs with correct URL format"""
        for items in self.items_list:
            with self.subTest(items=items):
                urls = get_urls(items=items)
                urls_match_regex = all([self.url_regex.match(url) for url in urls])
                self.assertTrue(urls_match_regex)

    def test_get_urls_length(self):
        """Function get_urls returns a list with correct length"""
        for items in self.items_list:
            with self.subTest(items=items):
                urls = get_urls(items=items)
                self.assertIsInstance(urls, list)
                self.assertEqual(len(urls), self.n_urls)

    def test_get_urls_websitedown(self):
        """If website is down, function get_urls() raises WebsiteDownError"""
        self.mock_get_max_pages.side_effect = WebsiteDownError("url")
        for items in self.items_list:
            with self.subTest(items=items):
                with self.assertRaises(WebsiteDownError):
                    get_urls(items=items)

    def test_build_urls(self):
        """Function build_urls() returns URLs with correct URL format"""
        for items in self.items_list:
            with self.subTest(items=items):
                urls = build_urls(items=items)
                urls_match_regex = all([self.url_regex.match(url) for url in urls])
                self.assertTrue(urls_match_regex)

    def test_build_urls_length(self):
        """Function build_urls() returns a list of URLsi with correct length"""
        for items in self.items_list:
            with self.subTest(items=items):
                urls = build_urls(items=items)
                self.assertIsInstance(urls, list)
                self.assertEqual(len(urls), self.n_urls)


class TestScrape(unittest.IsolatedAsyncioTestCase):
    n_urls = 1
    items_list = ["author_urls", "paper_urls", "project_urls", "group_urls"]

    def setUp(self):
        patcher = patch("src.scrape.get_max_pages")
        self.addCleanup(patcher.stop)
        self.mock_get_max_pages = patcher.start()
        self.mock_get_max_pages.return_value = self.n_urls

    async def test_scrape_author_urls(self):
        """Returns a list of URLs for the authors pages"""
        for items in self.items_list:
            with self.subTest(items=items):
                with open(f"tests/{items}.html", "rb") as f:
                    raw_html = f.read()
                with open(f"tests/{items}.txt", "r") as f:
                    url = f.read()

                html = HTML(url=url, html=raw_html)
                batch_urls = [[url]]

                with patch("src.scrape.retry_url") as mock_retry_url:
                    mock_retry_url.return_value.html = html
                    result = await scrape(
                        items=items,
                        urls=batch_urls,
                    )

                    self.assertEqual(len(result), 300)

                    if items == "author_urls":
                        url_regex = re.compile(
                            r"^orcid/[0-9X]{4}-[0-9X]{4}-[0-9X]{4}-[0-9X]{4}\b"
                        )
                        urls_match_regex = all([url_regex.match(url) for url in result])
                        self.assertTrue(urls_match_regex)


if __name__ == "__main__":
    unittest.main()
