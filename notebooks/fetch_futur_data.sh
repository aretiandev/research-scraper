#!/usr/bin/env bash
#===============================================================================
# fetch_data.sh
# Description:
#   Fetches data from researchers from the UPC website (https://futur.upc.edu/)
# Usage:
#   bash fetch_data.sh EETAC   -->  Get full list of authors from EETAC School
#   bash fetch_data.sh 181812  -->  Get information from author with id 181812
#===============================================================================

curl 'https://futur.upc.edu/api/items/search' \
  -H 'Accept: application/json, text/plain, */*' \
  -H 'Accept-Language: en-US,en;q=0.9,es;q=0.8,ko;q=0.7' \
  -H 'Client-ID: FUTUR-20031533201504' \
  -H 'Connection: keep-alive' \
  -H 'Content-Type: application/json;charset=UTF-8' \
  -H 'Cookie: futur.public.language=ca; _ga=GA1.1.1830310843.1704815841; _ga_TE059NZE8L=GS1.2.1704835989.2.0.1704835989.60.0.0; _ga_VVK66TK9MG=GS1.1.1704835984.2.1.1704836007.0.0.0; _ga_P80Z250TEZ=GS1.1.1704835984.2.1.1704836007.37.0.0; _ga_C8L3D4GXB0=GS1.1.1704836047.2.0.1704836047.0.0.0; _ga_EB0WDYRL2H=GS1.1.1704836044.2.0.1704836049.0.0.0; _ga_KD335TP0LJ=GS1.1.1705073741.4.1.1705074580.0.0.0' \
  -H 'Origin: https://futur.upc.edu' \
  -H "Referer: https://futur.upc.edu/'"$1"'" \
  -H 'Sec-Fetch-Dest: empty' \
  -H 'Sec-Fetch-Mode: cors' \
  -H 'Sec-Fetch-Site: same-origin' \
  -H 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36' \
  -H 'sec-ch-ua: "Not_A Brand";v="8", "Chromium";v="120", "Google Chrome";v="120"' \
  -H 'sec-ch-ua-mobile: ?0' \
  -H 'sec-ch-ua-platform: "macOS"' \
  --data-raw '{"query":"accessURL_str:(https\\:\\/\\/futur.upc.edu\\/'"$1"'##* OR https\\:\\/\\/futurdev.upc.edu\\/'"$1"'##*) OR activity_ca_doi_str:\"'"$1"'\" OR activity_ca_isbn_str:\"'"$1"'\" OR id_str:\"'"$1"'\" OR person_ca_orcid_str:\"'"$1"'\" OR person_ca_scopusId_str:\"'"$1"'\" OR journal_ca_issn_str:\"'"$1"'\" OR person_ca_gauss_str:\"'"$1"'\" OR person_ca_isiwok_str:\"'"$1"'\" OR organization_ca_code_str:\"'"$1"'\" OR organization_ca_urlAlias_str:\"'"$1"'\" OR organization_ca_acronym_str:\"'"$1"'\" OR person_ca_fullName_str:\"'"$1"'\"","facets":[],"page":{"pageIdx":0,"pageSize":2},"sort":"subtype desc"}' \
  --compressed