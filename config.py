import os

from dotenv import load_dotenv

from scripts.src.process import get_date

load_dotenv()


class Config:
    DATE = os.environ.get("DATE") or get_date()
    INSTITUTION_LIST = [
        "IGTP+",
        # "UPC_CIMNE",
        # "UB",
        # "UPF",
        # "UVic-UCC",
        # "UOC",
        # "Agrotecnio",
        # "CRAG",
        # "UdL",
        # "URV",
        # "UdG",
    ]
    THREADS_MAX = 16
    TIMEOUT = 1
    SLACK_BOT_TOKEN = os.environ.get("SLACK_BOT_TOKEN")
    SLACK_MEMBER_ID = os.environ.get("SLACK_MEMBER_ID")
    SIGNING_SECRET = os.environ.get("SIGNING_SECRET")
