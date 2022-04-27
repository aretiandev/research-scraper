import unittest
from unittest.mock import patch
import re

from scripts.scrape import build_urls, get_urls, WebsiteDownError


class TestAuthor(unittest.TestCase):
    """Test Author data, objects and methods"""

    url_regex = re.compile(
        r"^https://portalrecerca.csuc.cat/simple-search\?query=.*&start=[0-9]*"
    )
    n_urls = 30

    def setUp(self):
        """Setup mocks for functions that trigger HTTP requests"""
        patcher = patch("scripts.src.scrape.get_max_pages")
        self.addCleanup(patcher.stop)
        self.mock_get_max_pages = patcher.start()
        self.mock_get_max_pages.return_value = self.n_urls

    def test_get_urls(self):
        """URLs have URL format"""
        urls = get_urls(items="author_urls")
        urls_match_regex = all([self.url_regex.match(url) for url in urls])
        self.assertTrue(urls_match_regex)

    def test_get_urls_length(self):
        """URLs is a list with correct length"""
        urls = get_urls(items="author_urls")
        self.assertIsInstance(urls, list)
        self.assertEqual(len(urls), self.n_urls)

    def test_get_urls_websitedown(self):
        """If website is down, raise WebsiteDownError"""
        self.mock_get_max_pages.side_effect = WebsiteDownError("url")
        with self.assertRaises(WebsiteDownError):
            get_urls(items="author_urls")

    def test_build_urls(self):
        """URLs have URL format"""
        urls = build_urls(items="author_urls")
        urls_match_regex = all([self.url_regex.match(url) for url in urls])
        self.assertTrue(urls_match_regex)

    def test_build_urls_length(self):
        """URLs is a list with correct length"""
        urls = build_urls(items="author_urls")
        self.assertIsInstance(urls, list)
        self.assertEqual(len(urls), self.n_urls)


if __name__ == "__main__":
    unittest.main()
