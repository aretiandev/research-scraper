import unittest
import re

# from unittest.mock import patch
from scripts.src.scrape import get_urls


class TestScraper(unittest.TestCase):
    """Test Scraper methods"""

    def setUp(self):
        """Setup URLs"""
        self.urls = get_urls(items="author_urls")
        pass

    def test_author_urls_format(self):
        """Check all author URLs have URL format."""
        regex = re.compile(
            "^https://portalrecerca.csuc.cat/simple-search\?query=.*&start=[0-9]*"
        )
        result = all([regex.match(url) for url in self.urls])
        self.assertTrue(result)

    def test_author_urls_length(self):
        """Check author URLs is a list with more than 100 items"""
        self.assertIsInstance(self.urls, list)
        self.assertGreater(len(self.urls), 30)


if __name__ == "__main__":
    unittest.main()
