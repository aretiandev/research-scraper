#!/user/bin/env python
"""
Helper functions for Scraping Portal de la Reserca.

Main functions:
    scrape(): Main entrypoint
    scrape_url()
    scrape_author()
    scrape_project()
    scrape_group()
"""

import asyncio
import datetime
import logging
import os
import sys
import time

import pandas as pd
import requests.adapters
from requests_html import AsyncHTMLSession, HTMLSession

from .sqlite import (
    insert_authors,
    insert_groups,
    insert_papers,
    insert_projects,
    insert_urls,
)

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))

log = logging.getLogger(__name__)


class WebsiteDownError(Exception):
    def __init__(self, url):
        self.url = url


def get_max_pages(url):
    """Get max pages from pagination box in footer."""
    log.info("Retrieving number of URLs to scrape:")
    session = HTMLSession()
    r = session.get(url)
    if r.url == "https://portalrecerca.manteniment.csuc.cat":
        raise WebsiteDownError(url)
    pagination_items = r.html.find("div.discovery-result-pagination ul.pagination li")

    if len(pagination_items) == 0:
        return 0
    elif len(pagination_items) == 1:
        return 1

    max_pages_str = (
        pagination_items[-2]
        .text.split("\n")[0]
        .strip()
        .replace(".", "")
        .replace(",", "")
    )
    max_pages = int(max_pages_str)
    log.info(f"{max_pages:,d}.")
    return max_pages


async def retry_url(s, url, attempts=10, selector=None):
    """Retry fetching URL a given number of times."""
    for attempt in range(attempts):
        try:
            log.debug(
                f"""URL: {url}
Attempt: {attempt}"""
            )

            r = await s.get(url)

            if r.url == "https://portalrecerca.manteniment.csuc.cat":
                raise WebsiteDownError(url)

            log.debug(
                f"""URL: {url}
Attempt: {attempt}
Response: {r}"""
            )

            # If selector supplied, break only if it is found
            if selector:
                target = r.html.find(selector, first=True)
                if target:
                    break
            else:
                break
        except WebsiteDownError:
            raise
        except Exception as e:
            print("Exception: ", e)
            time.sleep(1)
            pass
    else:
        log.info(f"No attempts left for url: {url}")
        r = None

    return r


async def scrape_author(s, url, item="profile", attempts=10):
    """Scrape author page.

    Args:
        item: author tab to scrape: [profile, affiliation, project, group]
    Returns:
        result (dict): author information.
    """
    if item == "profile":
        # TODO: automate locale as a url_template to be added to all URLs
        url = url + "?locale=en"
        r = await retry_url(s, url, attempts)

        title = r.html.find("title", first=True).text
        title_error = "Portal de la Recerca de Catalunya: Malformed Request"
        if title_error in title:
            result = {"status_description": "System Error: Malformed Request"}
            return result

        selector = "div#collapseOneresearcherprofile table.table tr"
        rows = r.html.find(selector)
        scrape_dict = {
            "label": "Surnames, Name",
            "id": "ORCID",
            "institution_2": "Universities or Research centres",
        }

        try:
            # TODO: modularize the following code to a separate get_table function
            table = {}
            for row in rows:
                columns = row.find("td")
                columns1_text = columns[1].text.split("\n")
                if len(columns1_text) == 1:
                    table[columns[0].text] = columns1_text[0]
                else:  # if more than one university
                    table[columns[0].text] = columns1_text
            result = {}

            for label, pattern in scrape_dict.items():
                result[label] = table[pattern]
        except KeyError:
            result = {"error": "KeyError"}
        except IndexError:
            result = {"error": "IndexError"}
        except Exception:
            log.exception(f"There was an exception in scrape_author(), url: {url}")
            raise

    elif item == "affiliation":
        url = url + "/researcherdepartaments.html?onlytab=true"
        r = await retry_url(s, url, attempts)
        selector = "table.table tr"
        rows = r.html.find(selector)
        result = {}
        try:
            departments = []
            institutions = []
            for row in rows:
                columns = row.find("td")
                departments.append(columns[0].text)
                institutions.append(columns[1].text)
            result["department"] = departments
            result["institution"] = institutions
        except Exception:  # rows object is empty
            pass

    elif item == "project":
        url = url + "/publicresearcherprojects.html?onlytab=true"
        r = await retry_url(s, url, attempts)
        selector = "table.table tr"
        rows = r.html.find(selector)
        result = []
        for row in rows:
            project = row.find("td a")[0].attrs["href"]
            result.append(project)

    elif item == "group":
        url = url + "/orgs.html?onlytab=true"
        r = await retry_url(s, url, attempts)
        selector = "table.table tr"
        rows = r.html.find(selector)
        result = []
        for row in rows:
            group = row.find("td a")[0].attrs["href"]
            result.append(group)

    return result


async def scrape_project(s, url, tab="information", attempts=10):
    """Scrape project page.

    Args:
        tab: project tab to scrape: [information, researchers]
    Returns:
        project (dict): project information.
    """
    if tab == "information":
        url = url + "?onlytab=true&locale=en"

        r = await retry_url(s, url, attempts)

        tables_lst = pd.read_html(r.text)
        table = tables_lst[0].set_index(0).to_dict()[1]

        attr_keys = {
            "title": "Title",
            "official_code": "Official code",
            "url": "URL",
            "start_date": "Start date",
            "end_date": "End date",
            "institution": "Universities or CERCA centres",
        }

        project = {
            attr_name: table.get(table_key)
            for attr_name, table_key in attr_keys.items()
        }

    elif tab == "researchers":
        url = url + "/researchersprj.html?onlytab=true"
        r = await retry_url(s, url, attempts)

        selector = "table.table"
        tables = r.html.find(selector)

        project = {}

        # Principal researchers
        rows = tables[0].find("tr")
        principal_names = []
        principal_ids = []

        for row in rows[1:]:
            name_object = row.find("td", first=True)
            principal_names.append(name_object.text)
            try:
                orcid = name_object.find("a", first=True).attrs["href"][7:]
                principal_ids.append(orcid)
            except Exception:
                pass

        project["principal_names"] = principal_names
        project["principal_ids"] = principal_ids

        # Researchers
        try:
            rows = tables[1].find("tr")
            res_names = []
            res_ids = []

            for row in rows[1:]:
                name_object = row.find("td", first=True)
                res_names.append(name_object.text)
                try:
                    orcid = name_object.find("a", first=True).attrs["href"][7:]
                    res_ids.append(orcid)
                except Exception:
                    pass

            project["researcher_names"] = res_names
            project["researcher_ids"] = res_ids
        except Exception:
            pass

    return project


# Helper for group page
async def scrape_group(s, url, tab="information", attempts=10):
    """Scrape group page.

    Args:
        tab: group tab to scrape: [information, researchers]
    Returns:
        group (dict): group information
    """
    group = {}
    not_found_msg = "Nothing found"

    if tab == "information":
        url = url + "?onlytab=true&locale=en"
        r = await retry_url(s, url, attempts)

        try:
            title = r.html.find("title", first=True)
            tables_lst = pd.read_html(r.text)
        except ValueError:
            group["status_code"] = r.status_code
            try:
                if not_found_msg in title.text:
                    group["status_description"] = "Title: not found"
            except AttributeError:
                pass

            return group

        table = tables_lst[0].set_index(0).to_dict()[1]

        attr_keys = {
            "name": "Name",
            "acronym": "Acronym",
            "institution": "Universities or CERCA centres",
            "group_url": "URL",
            "sgr_code": "SGR code",
        }

        group = {
            attr_name: table.get(table_key)
            for attr_name, table_key in attr_keys.items()
        }

    elif tab == "researchers":
        r = await retry_url(s, url, attempts)

        try:
            title = r.html.find("title", first=True)
            next_tab_href = r.html.find("#bar-tab-1012900 a", first=True).attrs["href"]
        except AttributeError:
            group["status_code"] = r.status_code
            try:
                if not_found_msg in title.text:
                    group["status_description"] = "Title: not found"
            except Exception:
                pass

            return group

        url_root = "https://portalrecerca.csuc.cat"
        url = url_root + next_tab_href

        r = await retry_url(s, url, attempts)

        selector = "table.table"
        tables = r.html.find(selector)

        # Principal researchers
        rows = tables[0].find("tr")

        principal_names = []
        principal_ids = []

        for row in rows:
            name_object = row.find("td", first=True)
            principal_names.append(name_object.text)
            try:
                orcid = name_object.find("a", first=True).attrs["href"][7:]
            except KeyError:
                orcid = ""
            principal_ids.append(orcid)

        group["principal_names"] = principal_names
        group["principal_ids"] = principal_ids

        # Researchers
        try:
            rows = tables[1].find("tr")
            res_names = []
            res_ids = []

            for row in rows:
                name_object = row.find("td", first=True)
                res_names.append(name_object.text)
                try:
                    orcid = name_object.find("a", first=True).attrs["href"][7:]
                except KeyError:
                    orcid = ""
                res_ids.append(orcid)

            group["researcher_names"] = res_names
            group["researcher_ids"] = res_ids
        except IndexError:
            pass

    return group


# Scrape single URL
async def scrape_url(s, url, items="author_data", attempts=10):
    """Scrape URL.

    Args:
        s: instance of requests_html.AsyncHTMLSession().
        url: URL to scrape.
        items: paper_urls, paper_data, author_urls, author_data,
               project_urls, project_data, group_urls, group_data.
    Returns:
        (dict): item information.
    """

    if items == "author_data":
        result = await asyncio.gather(
            scrape_author(s, url, "profile"),
            scrape_author(s, url, "affiliation"),
            scrape_author(s, url, "project"),
            scrape_author(s, url, "group"),
        )

        author = result[0]  # name and id
        try:
            author["department"] = result[1]["department"]
        except KeyError:
            pass
        try:
            author["institution"] = result[1]["institution"]
        except KeyError:
            pass
        author["projects"] = result[2]
        author["groups"] = result[3]
        author["url"] = url

        return author

    elif items == "project_data":
        result = await asyncio.gather(
            scrape_project(s, url, tab="information", attempts=attempts),
            scrape_project(s, url, tab="researchers", attempts=attempts),
        )

        project = {}
        for d in result:
            project.update(d)

        project["url"] = url
        project["url_stem"] = url[32:].split("?")[0]

        return project

    elif items == "group_data":
        try:
            result = await asyncio.gather(
                scrape_group(s, url, tab="information", attempts=attempts),
                scrape_group(s, url, tab="researchers", attempts=attempts),
            )

            group = {}
            for d in result:
                group.update(d)

            group["url"] = url
            group["url_stem"] = url[32:].split("?")[0]
        except Exception:
            log.exception(f"Exception in scrape_url({items}), url: {url}")
            group = []
            raise

        return group

    elif items == "paper_data":
        paper = {}
        paper["url"] = url
        paper["url_stem"] = url[32:].split("?")[0]

        not_found_msg = "No ha estat possible trobar el que esteu buscant"

        r = await retry_url(s, url, attempts)

        try:
            title = r.html.find("title", first=True)
            table = r.html.find("table.table", first=True)
            rows = table.find("tr")
        except Exception:
            paper["status_code"] = r.status_code
            try:
                if not_found_msg in title.text:
                    paper["status_description"] = "Title: not found"
            except Exception:
                pass

            return paper

        for row in rows:
            # Get columns in row
            columns = row.find("td")

            # Skip if empty row
            if len(columns) == 0:
                continue

            attributes = {
                "dc.contributor.authors": "author",
                "dc.date.issued": "date",
                "dc.publisher": "publisher",
                "dc.identifier.citation": "citation",
                "dc.identifier.issn": "issn",
                "dc.identifier.uri": "uri",
                "dc.identifier.isbn": "isbn",
                "dc.relation.ispartof": "published_in",
                "dc.title": "title",
                "dc.type": "type",
                "dc.identifier.doi": "doi",
                "dc.identifier.sourceid": "sourceid",
                "dc.identifier.sourceref": "sourceref",
                "Appears in Collections:": "appears_in_collections",
            }

            for row_index in attributes.keys():
                if columns[0].text == row_index:
                    paper[attributes[row_index]] = columns[1].text

        # Get authors ORCid'ss
        simple_url = url.split("?")[0] + "?mode=simple"

        r = await retry_url(s, simple_url, attempts)

        try:
            title = r.html.find("title", first=True)
            table = r.html.find("table.table", first=True)
            rows = table.find("tr")
            authors = r.html.find("table.table a.authority.author")
            author_hrefs = []
            for author in authors:
                href = author.attrs["href"]
                if href[:7] == "/orcid/":
                    href = href[7:]
                author_hrefs.append(href)
        except Exception:
            paper["status_code"] = r.status_code
            try:
                if not_found_msg in title.text:
                    paper["status_description"] = "Title: not found"
            except Exception:
                pass

            return paper

        paper["orcids"] = author_hrefs

        return paper

    else:  # paper_urls, author_urls, project_urls or group_urls
        selector = "div.panel.panel-info table.table"
        r = await retry_url(
            s, url, attempts=20, selector=selector
        )  # 20 attempts avoid errors

        try:
            table = r.html.find("div.panel.panel-info table.table", first=True)
            rows = table.find("tr")

            result = []

            for i, row in enumerate(rows):
                if i >= 300:  # Hardcoded: ignore results after row 300
                    continue
                # Get columns in row
                columns = row.find("td")

                # Skip if empty row
                if len(columns) == 0:
                    continue

                if items in ["paper_urls", "project_urls"]:
                    scrape_item = columns[1].find("a")[0].attrs["href"][1:]
                elif items in ["author_urls", "group_urls"]:
                    scrape_item = columns[0].find("a")[0].attrs["href"][1:]

                # Append to results_list
                result.append(scrape_item)

        except Exception:
            log.exception(f"Exception in scrape_url({items}), url: {url}")
            result = []

        return result


def get_urls(items, institution, n_pages=None):
    """Get list of URLs to pass to scrape()

    Args:
        items: author_urls, paper_urls, project_urls, group_urls.
        n_pages: max pages to scrape.
    Returns:
        urls (list): list of URLs.
    """


    # testing
    
    url_template = (
        "https://portalrecerca.csuc.cat/simple-search?"
        + "query="
        + "&location={location}"
        + "&filter_field_1={filter_field_1}"
        + "&filter_type_1={filter_type_1}"
        + "&filter_value_1={filter_value_1}"
        + "&filter_field_2={filter_field_2}"
        + "&filter_type_2={filter_type_2}"
        + "&filter_value_2={filter_value_2}"
        + "&sort_by={sort_by}"
        + "&order={order}"
        + "&rpp={rpp}"
        + "&etal=0"
        + "&filtername=location.coll"
        + "&filtertype=equals"
        + "&filterquery={filterquery}" 
        + "&start="
    )
    
    institution_search_fields = {
        "ICFO": {"filterquery": "359"},
        "IDIBELL": {"filterquery": "320"},
        "UPC": {"filterquery": "305"},
        "IGTP+": {"filterquery": "349"},
        "UPC-CIMNE": {"filterquery": "337"},
        "UB": {"filterquery": "299"},
        "UPF": {"filterquery": "306"},
        "UVic-UCC": {"filterquery": "302"},
        "UOC": {"filterquery": "304"},
        "Agrotecnio": {"filterquery": "321"},
        "CRAG": {"filterquery": "342"},
        "UdL": {"filterquery": "301"},
        "URV": {"filterquery": "308"},
        "UdG": {"filterquery": "300"},
        "IRSJD": {"filterquery": "362"},
        "URL": {"filterquery": "307"},
        "UIC": {"filterquery": "303"},
    }

    base_search_fields = {
        "author_urls": {
            "location": "crisrp",
            "filter_field_1": "resourcetype",
            "filter_type_1": "equals",
            "filter_value_1": "Researchers",
            "filter_field_2": "",
            "filter_type_2": "",
            "filter_value_2": "",
            "sort_by": "crisrp.fullName_sort",
            "order": "asc",
            "rpp": "300",
        },
        "paper_urls": {
            "location": "publications",
            "filter_field_1": "resourcetype",
            "filter_type_1": "equals",
            "filter_value_1": "Items",
            "filter_field_2": "itemtype",
            "filter_type_2": "notequals",
            "filter_value_2": "Phd+Thesis",
            "sort_by": "dc.contributor.authors_sort",
            "order": "asc",
            "rpp": "300",
        },
        "project_urls": {
            "location": "crisproject",
            "filter_field_1": "resourcetype",
            "filter_type_1": "equals",
            "filter_value_1": "Projects",
            "filter_field_2": "",
            "filter_type_2": "",
            "filter_value_2": "",
            "sort_by": "crisrp.title_sort",
            "order": "asc",
            "rpp": "300",
        },
        "group_urls": {
            "location": "crisou",
            "filter_field_1": "resourcetype",
            "filter_type_1": "equals",
            "filter_value_1": "OrgUnits",
            "filter_field_2": "",
            "filter_type_2": "",
            "filter_value_2": "",
            "sort_by": "crisou.name_sort",
            "order": "asc",
            "rpp": "300",
        },
    }

    # Get the base search fields for the items
    search_fields = base_search_fields.get(items, {}).copy()
    
    # Update the search fields with the institution-specific filterquery
    if institution in institution_search_fields:
        search_fields.update(institution_search_fields[institution])

    url_root = url_template.format(**search_fields)

    if not n_pages:
        n_pages = get_max_pages(url_root + "0")

    urls = [url_root + str(page * 300) for page in range(n_pages)]

    return urls



# Main Entrypoint
async def scrape(
    items,
    urls,
    batch_start=0,
    out_file=None,
    out_sql=False,
    timeout=None,
    database="recerca.db",
):
    """Main entry function to scrape Portal de la Reserca.

    Args:
        items:
            author_data, paper_data, project_data, group_data,
            author_urls, paper_urls, project_urls, group_urls.
        urls: list of urls.
        start_pos: URL starting position.
        batch_start_pos: Batch starting position.
        n_pages: max pages to scrape.
        batch_size: batch size.
        out_file: output file.
    Returns:
        result (list of dicts): scraping results.
        out_file (csv): writes DataFrame to csv.
    """

    urls_flat = [url for batch in urls for url in batch]
    batch_size = len(urls[0])

    log.info(
        f"""Scraping items: {items}
    URL count: {len(urls_flat):,d}
    Batch count: {len(urls):,d}
    Batch length: {batch_size:,d}
    Starting batch: {batch_start}
    Output: {out_file}"""
    )

    if out_file:
        result_df = pd.DataFrame()

    for i, batch in enumerate(urls[batch_start:]):
        s = AsyncHTMLSession()
        # https://stackoverflow.com/a/18845952/10688326
        adapter = requests.adapters.HTTPAdapter(pool_connections=100, pool_maxsize=100)
        s.mount("https://", adapter)

        t1 = time.perf_counter()

        tasks = (scrape_url(s, url, items=items) for url in batch)
        try:
            batch_result = await asyncio.gather(*tasks)
        except Exception:
            log.exception(f"Exception in scrape({items}) for batch number {i}:")
            raise

        t2 = time.perf_counter()

        if items in ["author_urls", "paper_urls", "project_urls", "group_urls"]:
            batch_result = [
                i for sublist in batch_result for i in sublist
            ]  # Flatten result

        if out_file:
            result_df = pd.concat(
                [result_df, pd.DataFrame(batch_result)], ignore_index=True
            )
            result_df.to_csv(out_file, index=None)
            log.debug(f"Saved all data up to batch {i}: {out_file}")

        if out_sql:
            log.debug(f"Writing to database. Batch: {i}. # URLS: {len(batch_result)}.")
            date_today = out_file.split("/")[-1][:8] if out_file else None
            if items in ["paper_urls", "author_urls", "project_urls", "group_urls"]:
                insert_urls(batch_result, items, date_today, db=database)
            elif items == "paper_data":
                insert_papers(batch_result, date_today)
            elif items == "author_data":
                insert_authors(batch_result, date_today)
            elif items == "group_data":
                insert_groups(batch_result, date_today)
            elif items == "project_data":
                insert_projects(batch_result, date_today)

        # Log estimated time left
        seconds_left = (len(urls) - i) * (t2 - t1)
        batch_time = str(datetime.timedelta(seconds=round(t2 - t1)))
        time_left = str(datetime.timedelta(seconds=round(seconds_left)))

        print(
            f"Progress: {(i+1)/len(urls)*100:.0f}% ({i+1}/{len(urls):,d}). "
            + f"URLs: {i*batch_size}-{(i+1)*batch_size-1}. "
            + f"Batch time: {batch_time}. Time left: {time_left}.  ",
            end="\r",
        )

        if timeout:
            time.sleep(timeout)

    log.info("\nDone.")
    return batch_result
