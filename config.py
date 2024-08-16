import os

from dotenv import load_dotenv

from scripts.src.process import get_date

load_dotenv()


class Config:
    # DATE = os.environ.get("DATE") or get_date()
    DATE = "20240317"
    INSTITUTION_LIST = [
        #"ICFO",
        #"IDIBELL",
        #"UPC",
        ##"IGTP+",
        #"UPF",
        #"UVic-UCC",
        #"UOC",
        #"Agrotecnio"
        #"CRAG",
        #"UdL",
        #"URV",
        #"UdG",
        #"IRSJD",
        #"URL",
        #"UIC",
        
        "UB",
        #"UPC-CIMNE",
        #VHIR,
        #IDIBAPS,
        #CSIC,
        #RDR,
        #IRSantPau,
        #IDIAPJGol,
        #IMIM,
        #IRBLleida,
        #ICN2,
        #I3PT,
        #IBEC,
        #IRBBarcelona,
        #IEEC,
        #IISPV,
        #IRTA,
        #IDIBGI,
        #CREAF,
        #IREC,
        #ICIQ,
        #IFAE,
        #ISGlobal,
        #ICRA,
        #IPHES,
        #CTFC,
        #CRG,
        #IrsiCaixa,
        #ICP,
        #CRM,
        #ICAC,
        #TecnoCampus,
        #i2CAT,
        #IJC,
        #CED,
        #INEFC,
        #CTTC,
        #UAO,
        #VHIO,
        #IBEI,
        #ICRPC,
        #CREI
    ]
    THREADS_MAX = 16
    DATABASE = "espluges.db"
    SLACK_BOT_TOKEN = os.environ.get("SLACK_BOT_TOKEN")
    SLACK_MEMBER_ID = os.environ.get("SLACK_MEMBER_ID")
    SIGNING_SECRET = os.environ.get("SIGNING_SECRET")
    BATCH_SIZE = 50
    TIMEOUT = 1
