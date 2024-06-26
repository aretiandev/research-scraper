import os

from dotenv import load_dotenv

from scripts.src.process import get_date

load_dotenv()


class Config:
    # DATE = os.environ.get("DATE") or get_date()
    DATE = "20240317"
    INSTITUTION_LIST = [
        # "ICFO",
        # "IDIBELL",
        "UPC",
        # "IGTP+",
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
        # "IRSJD",
        # "URL",
        # "UIC",
    ]
    THREADS_MAX = 16
    DATABASE = "espluges.db"
    SLACK_BOT_TOKEN = os.environ.get("SLACK_BOT_TOKEN")
    SLACK_MEMBER_ID = os.environ.get("SLACK_MEMBER_ID")
    SIGNING_SECRET = os.environ.get("SIGNING_SECRET")
    BATCH_SIZE = 50
    TIMEOUT = 1
