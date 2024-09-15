import os

from dotenv import load_dotenv

from scripts.src.process import get_date

load_dotenv()


class Config:
    # DATE = os.environ.get("DATE") or get_date()
    DATE = "20240905"
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
        #"UB",
        
        #"VHIR",
        #"IRSantPau",
        
        #"IMIM",
        #"IRBLleida",
        #"ICN2",
        #"I3PT",
        #"IDIBAPS",
        #"CSIC",
        #"RDR",
        #"IDIAPJGol",
        #"IBEC",
        #"IRBBarcelona",
        #"IEEC",
        #"IISPV",
        #"IRTA",
        #"IDIBGI",
        
        #"CREAF",
        #"IREC",
        #"ICIQ",
        
         #"IFAE",
         #"ISGlobal",
         #"ICRA",
         #"IPHES",
         #"CTFC",
        
        #"CRG",
        
        #"IrsiCaixa",
        
        #"ICP",
        #"CRM",
        
        #"ICAC",
        #"TecnoCampus",
        #"i2CAT",
        #"IJC",
        #"CED",
        #"INEFC",
        
        #"CTTC",
        
        #"UAO",
        #"VHIO",
        #"IBEI",
        #"ICRPC",
        #"CREI"
        
        "UPC-CIMNE",
    ]
    
    THREADS_MAX = 16
    DATABASE = "espluges.db"
    SLACK_BOT_TOKEN = os.environ.get("SLACK_BOT_TOKEN")
    SLACK_MEMBER_ID = os.environ.get("SLACK_MEMBER_ID")
    SIGNING_SECRET = os.environ.get("SIGNING_SECRET")
    BATCH_SIZE = 50
    TIMEOUT = 1
