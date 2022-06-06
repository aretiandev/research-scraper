#!/usr/bin/env python
# coding: utf-8
"""
Utils module:
Collection of utils
"""
import logging

log = logging.getLogger(__name__)


def get_date():
    """Get formatted today's date."""
    import os
    import time
    from datetime import date

    os.environ["TZ"] = "America/New_York"
    time.tzset()
    date_today = date.today().strftime("%Y%m%d")
    return date_today
