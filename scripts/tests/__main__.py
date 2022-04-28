import sys
import unittest
from .create_pickles import main

# No arguments: run tests
if len(sys.argv) == 1:
    from .test_scraper import *  # noqa

    unittest.main()

# pickle + items: pickle objects
elif len(sys.argv) <= 3:
    if sys.argv[1] == "pickle":
        try:
            items = sys.argv[2]
        except IndexError:
            items = None
        main(items)
    else:
        sys.exit("Wrong arguments.")

else:
    sys.exit("Too many arguments.")
