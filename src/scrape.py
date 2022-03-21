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

import pandas as pd
from requests_html import HTMLSession, AsyncHTMLSession
import datetime
import time
import asyncio
import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))

from src.logging import create_logger

log = create_logger(__name__, f"log/{__name__}.log")


class WebsiteDownException(Exception):
    pass


def get_max_pages(url):
    """Get max pages from pagination box in footer."""
    log.info("Retrieving number of URLs to scrape:")
    session = HTMLSession()
    r = session.get(url)
    if r.url == 'https://portalrecerca.manteniment.csuc.cat':
        raise WebsiteDownException("Website is down.")
    pagination_items = r.html.find('div.discovery-result-pagination ul.pagination li')
    max_pages_str = pagination_items[-2].text.split("\n")[0].strip().replace('.', '').replace(',', '')
    max_pages = int(max_pages_str)
    log.info(f"{max_pages:,d}.")
    return max_pages


async def retry_url(s, url, attempts=10, selector=None, debug=False):
    """Retry fetching URL a given number of times."""
    for _ in range(attempts):
        try:
            r = await s.get(url)

            if r.url == 'https://portalrecerca.manteniment.csuc.cat':
                raise WebsiteDownException("Website is down.")

            log.debug(f"{r} Attempt: {_}. URL: {url}")

            # If selector supplied, break only if it is found
            if selector:
                target = r.html.find(selector, first=True)
                if target:
                    break
            else:
                break
        except WebsiteDownException:
            raise
        except Exception:
            time.sleep(1)
            pass
    else:
        log.info(f"No attempts left for url: {url}")
        r = None

    return r


async def scrape_author(s, url, item='profile', attempts=10):
    """Scrape author page.

    Args:
        item: author tab to scrape: [profile, affiliation, project, group]
    Returns:
        result (dict): author information.
    """
    if item == 'profile':
        # TODO: automate locale as a url_template to be added to all URLs
        url = url + "?locale=en"
        r = await retry_url(s, url, attempts)

        title = r.html.find('title', first=True).text
        title_error = 'Portal de la Recerca de Catalunya: Malformed Request'
        if title_error in title:
            result = {'status description': 'System Error: Malformed Request'}
            return result

        selector = 'div#collapseOneresearcherprofile table.table tr'
        rows = r.html.find(selector)
        scrape_dict = {
            'label': 'Surnames, Name',
            'id': 'ORCID',
            'institution_2': 'Universities or CERCA centres'
        }

        try:
            # TODO: modularize the following code to a separate get_table function
            table = {}
            for row in rows:
                columns = row.find('td')
                columns1_text = columns[1].text.split("\n")
                if len(columns1_text) == 1:
                    table[columns[0].text] = columns1_text[0]
                else:  # if more than one university
                    table[columns[0].text] = columns1_text
            result = {}

            for label, pattern in scrape_dict.items():
                result[label] = table[pattern]
        except IndexError:
            result = None

    elif item == 'affiliation':
        url = url + '/researcherdepartaments.html?onlytab=true'
        r = await retry_url(s, url, attempts)
        selector = 'table.table tr'
        rows = r.html.find(selector)
        result = {}
        try:
            departments = []
            institutions = []
            for row in rows:
                columns = row.find('td')
                departments.append(columns[0].text)
                institutions.append(columns[1].text)
            result['department'] = departments
            result['institution'] = institutions
        except Exception:  # rows object is empty
            pass

    elif item == 'project':
        url = url + '/publicresearcherprojects.html?onlytab=true'
        r = await retry_url(s, url, attempts)
        selector = 'table.table tr'
        rows = r.html.find(selector)
        result = []
        for row in rows:
            project = row.find('td a')[0].attrs['href']
            result.append(project)

    elif item == 'group':
        url = url + '/orgs.html?onlytab=true'
        r = await retry_url(s, url, attempts)
        selector = 'table.table tr'
        rows = r.html.find(selector)
        result = []
        for row in rows:
            group = row.find('td a')[0].attrs['href']
            result.append(group)

    return result


async def scrape_project(s, url, tab='information', attempts=10):
    """Scrape project page.

    Args:
        tab: project tab to scrape: [information, researchers]
    Returns:
        project (dict): project information.
    """
    if tab == 'information':
        url = url + '?onlytab=true&locale=en'

        r = await retry_url(s, url, attempts)

        tables_lst = pd.read_html(r.text)
        table = tables_lst[0].set_index(0).to_dict()[1]

        attr_keys = {
            'title': 'Title',
            'official code': 'Official code',
            'url': 'URL',
            'start date': 'Start date',
            'end date': 'End date',
            'institution': 'Universities or CERCA centres'}

        project = {attr_name:table.get(table_key) for attr_name,table_key in attr_keys.items()}

    elif tab == 'researchers':
        url = url + '/researchersprj.html?onlytab=true'
        r = await retry_url(s, url, attempts)

        selector = 'table.table'
        tables = r.html.find(selector)

        project = {}

        # Principal researchers
        rows = tables[0].find('tr')
        principal_names = []
        principal_ids = []

        for row in rows[1:]:
            name_object = row.find('td', first=True)
            principal_names.append(name_object.text)
            try:
                orcid = name_object.find('a', first=True).attrs['href'][7:]
                principal_ids.append(orcid)
            except Exception:
                pass

        project['principal names'] = principal_names
        project['principal ids'] = principal_ids

        # Researchers
        try:
            rows = tables[1].find('tr')
            res_names = []
            res_ids = []

            for row in rows[1:]:
                name_object = row.find('td', first=True)
                res_names.append(name_object.text)
                try:
                    orcid = name_object.find('a', first=True).attrs['href'][7:]
                    res_ids.append(orcid)
                except Exception:
                    pass

            project['researcher names'] = res_names
            project['researcher ids'] = res_ids
        except Exception:
            pass

    return project


# Helper for group page
async def scrape_group(s, url, tab='information', attempts=10):
    """Scrape group page.

    Args:
        tab: group tab to scrape: [information, researchers]
    Returns:
        group (dict): group information
    """
    if tab == 'information':
        url = url + '?onlytab=true&locale=en'
        r = await retry_url(s, url, attempts)

        tables_lst = pd.read_html(r.text)
        table = tables_lst[0].set_index(0).to_dict()[1]

        attr_keys = {
            'name': 'Name',
            'acronym': 'Acronym',
            'institution': 'Universities or CERCA centres',
            'group_url': 'URL',
            'sgr_code': 'SGR code'}

        group = {attr_name: table.get(table_key)
                 for attr_name, table_key in attr_keys.items()}

    elif tab == 'researchers':
        r = await retry_url(s, url, attempts)
        next_tab_href = r.html.find('#bar-tab-1012900 a', first=True).attrs['href']
        url_root = 'https://portalrecerca.csuc.cat' 
        url = url_root + next_tab_href

        r = await retry_url(s, url, attempts)

        selector = 'table.table'
        tables = r.html.find(selector)

        group = {}

        # Principal researchers
        rows = tables[0].find('tr')

        principal_names = []
        principal_ids = []

        for row in rows:
            name_object = row.find('td', first=True)
            principal_names.append(name_object.text)
            try:
                orcid = name_object.find('a', first=True).attrs['href'][7:]
            except KeyError:
                orcid = ""
            principal_ids.append(orcid)

        group['principal names'] = principal_names
        group['principal ids'] = principal_ids

        # Researchers
        try:
            rows = tables[1].find('tr')
            res_names = []
            res_ids = []

            for row in rows:
                name_object = row.find('td', first=True)
                res_names.append(name_object.text)
                try:
                    orcid = name_object.find('a', first=True).attrs['href'][7:]
                except KeyError:
                    orcid = ""
                res_ids.append(orcid)

            group['researcher names'] = res_names
            group['researcher ids'] = res_ids
        except IndexError:
            pass

    return group


# Scrape single URL
async def scrape_url(s, url, items='author', attempts=10):
    """Scrape URL.

    Args:
        s: instance of requests_html.AsyncHTMLSession().
        url: URL to scrape.
        items: paper_links, paper, author_links, author,
               project_links, project, group_links, group.
    Returns:
        (dict): item information.
    """

    if items == 'author':
        result = await asyncio.gather(
                scrape_author(s, url, 'profile'),
                scrape_author(s, url, 'affiliation'),
                scrape_author(s, url, 'project'),
                scrape_author(s, url, 'group')
            )

        author = result[0]  # name and id
        try:
            author['department'] = result[1]['department']
        except KeyError:
            pass
        try:
            author['institution'] = result[1]['institution']
        except KeyError:
            pass
        author['projects'] = result[2]
        author['groups'] = result[3]
        author['url'] = url

        return author

    elif items == 'project':
        result = await asyncio.gather(
                scrape_project(s, url, tab='information', attempts=attempts),
                scrape_project(s, url, tab='researchers', attempts=attempts),
            )

        project = {}
        for d in result:
            project.update(d)

        project['url']    = url
        project['url_id'] = url[31:].split('?')[0]

        return project

    elif items == 'group':
        result = await asyncio.gather(
                scrape_group(s, url, tab='information', attempts=attempts),
                scrape_group(s, url, tab='researchers', attempts=attempts),
            )

        group = {}
        for d in result:
            group.update(d)

        group['url']    = url
        group['url_id'] = url[31:].split('?')[0]

        return group

    elif items == 'paper':
        paper = {}
        paper['url']         = url
        paper['url_id']      = url[31:].split('?')[0]

        not_found_msg = 'No ha estat possible trobar el que esteu buscant'

        r = await retry_url(s, url, attempts)

        try:
            title = r.html.find('title', first=True)
            table = r.html.find('table.table', first=True)
            rows = table.find('tr')
        except Exception:
            paper['status code'] = r.status_code
            try:
                if not_found_msg in title.text:
                    paper['status description'] = 'Title: not found'
            except Exception:
                pass

            return paper

        for row in rows:
            # Get columns in row
            columns = row.find('td')

            # Skip if empty row
            if len(columns) == 0:
                continue

            attributes = {
                'dc.contributor.authors' : 'author',
                'dc.date.issued'         : 'date',
                'dc.publisher'           : 'publisher',
                'dc.identifier.citation' : 'citation',
                'dc.identifier.issn'     : 'issn',
                'dc.identifier.uri'      : 'uri',
                'dc.identifier.isbn'     : 'isbn',
                'dc.relation.ispartof'   : 'published_in',
                'dc.title'               : 'title',
                'dc.type'                : 'type',
                'dc.identifier.doi'      : 'doi',
                'dc.identifier.sourceid' : 'sourceid',
                'dc.identifier.sourceref': 'sourceref',
                'Appears in Collections:': 'appears_in_collections'}

            for row_index in attributes.keys():
                if columns[0].text == row_index:
                    paper[attributes[row_index]] = columns[1].text

        # Get authors ORCid'ss
        simple_url = url.split('?')[0] + '?mode=simple'

        r = await retry_url(s, simple_url, attempts)

        try:
            title = r.html.find('title', first=True)
            table = r.html.find('table.table', first=True)
            rows = table.find('tr')
            authors = r.html.find('table.table a.authority.author')
            author_hrefs = []
            for author in authors:
                href = author.attrs['href']
                if href[:7] == '/orcid/':
                    href = href[7:]
                author_hrefs.append(href)
        except Exception:
            paper['status code'] = r.status_code
            try:
                if not_found_msg in title.text:
                    paper['status description'] = 'Title: not found'
            except Exception:
                pass

            return paper

        paper['orcids'] = author_hrefs

        return paper

    else:  # paper_links, author_links, project_links or group_links
        selector = 'div.panel.panel-info table.table'
        r = await retry_url(s, url, attempts=20, selector=selector)  # 20 attempts avoid errors

        try:
            table = r.html.find('div.panel.panel-info table.table', first=True)
            rows = table.find('tr')

            result = []

            for i, row in enumerate(rows):
                if i >= 300:  # Hardcoded: ignore results after row 300
                    continue
                # Get columns in row
                columns = row.find('td')

                # Skip if empty row
                if len(columns) == 0:
                    continue

                if items in ['paper_links', 'project_links']:
                    scrape_item = columns[1].find('a')[0].attrs['href']
                elif items in ['author_links', 'group_links']:
                    scrape_item = columns[0].find('a')[0].attrs['href']

                # Append to results_list
                result.append(scrape_item)

        except Exception:
            log.exception(f"Exception in scrape_url({items}), url: {url}")
            result = []

        return result


# Main Entrypoint
async def scrape(
        items='author',
        urls=None,
        start_pos=0,
        batch_start_pos=0,
        n_pages=None,
        batch_size=None,
        out_file=None):
    """Main entry function to scrape Portal de la Reserca.

    Args:
        items:
            author, paper, project, group,
            author_links, paper_links, project_links, group_links.
        urls: list of urls. Not needed if items in [author_links, paper_links].
        start_pos: URL starting position.
        batch_start_pos: Batch starting position.
        n_pages: max pages to scrape.
        batch_size: batch size.
        out_file: output file.
    Returns:
        result (list of dicts): scraping results.
        out_file (csv): writes DataFrame to csv.
    """

    # Get list of URLs to scrape hyperlinks
    if items in ['author_links', 'paper_links', 'project_links', 'group_links']:

        if not urls:
            url_template =                                        \
                'https://portalrecerca.csuc.cat/simple-search?' + \
                'query='                                        + \
                '&location={location}'                          + \
                '&filter_field_1={filter_field_1}'              + \
                    '&filter_type_1={filter_type_1}'            + \
                    '&filter_value_1={filter_value_1}'          + \
                '&filter_field_2={filter_field_2}'              + \
                    '&filter_type_2={filter_type_2}'            + \
                    '&filter_value_2={filter_value_2}'          + \
                '&sort_by={sort_by}'                            + \
                    '&order={order}'                            + \
                '&rpp={rpp}'                                    + \
                '&etal=0'                                       + \
                '&start='

            if items == 'author_links':
                search_fields = {
                    'location'      : 'crisrp',
                    'filter_field_1': 'resourcetype',
                    'filter_type_1' : 'equals',
                    'filter_value_1': 'Researchers',
                    'filter_field_2': '',
                    'filter_type_2' : '',
                    'filter_value_2': '',
                    'sort_by'       : 'crisrp.fullName_sort',
                    'order'         : 'asc',
                    'rpp'           : '300'
                }

            elif items == 'paper_links':
                search_fields = {
                    'location'      : 'publications',
                    'filter_field_1': 'resourcetype',
                    'filter_type_1' : 'equals',
                    'filter_value_1': 'Items',
                    'filter_field_2': 'itemtype',
                    'filter_type_2' : 'notequals',
                    'filter_value_2': 'Phd+Thesis',
                    'sort_by'       : 'dc.contributor.authors_sort',
                    'order'         : 'asc',
                    'rpp'           : '300'
                }

            elif items == 'project_links':
                search_fields = {
                    'location'      : 'crisproject',
                    'filter_field_1': 'resourcetype',
                    'filter_type_1' : 'equals',
                    'filter_value_1': 'Projects',
                    'filter_field_2': '',
                    'filter_type_2' : '',
                    'filter_value_2': '',
                    'sort_by'       : 'crisrp.title_sort',
                    'order'         : 'asc',
                    'rpp'           : '300'
                }

            elif items == 'group_links':
                search_fields = {
                    'location'      : 'crisou',
                    'filter_field_1': 'resourcetype',
                    'filter_type_1' : 'equals',
                    'filter_value_1': 'OrgUnits',
                    'filter_field_2': '',
                    'filter_type_2' : '',
                    'filter_value_2': '',
                    'sort_by'       : 'crisou.name_sort',
                    'order'         : 'asc',
                    'rpp'           : '300'
                }

            url_root = url_template.format(**search_fields)

            if not n_pages:
                n_pages = get_max_pages(url_root + '0')

            urls = [url_root + str(page*300) for page in range(n_pages)]

    if not batch_size:
        batch_size = len(urls)

    batch_urls = [urls[i:i+batch_size] for i in range(start_pos, len(urls), batch_size)]

    log.info(f"""Scraping items: {items}
    URL count: {len(urls)-start_pos:,d}
    Batch count: {len(batch_urls):,d}
    Batch length: {batch_size:,d}
    Starting URL: {start_pos:,d}
    Starting batch: {batch_start_pos}
    Output: {out_file}""")

    if out_file:
        result_df = pd.DataFrame()

    result = []
    for i, batch in enumerate(batch_urls[batch_start_pos:]):

        s = AsyncHTMLSession()

        t1 = time.perf_counter()
        tasks = (scrape_url(s, url, items=items) for url in batch)
        try:
            batch_result = await asyncio.gather(*tasks)
        except Exception:
            log.exception(f"Exception in scrape({items}) for batch number {i}:")
            raise
        t2 = time.perf_counter()

        # Flatten result
        if items in ['author_links', 'paper_links', 'project_links', 'group_links']:
            batch_result = [i for sublist in batch_result for i in sublist]

        result.extend(batch_result)

        if out_file:
            result_df = result_df.append(batch_result, ignore_index=True)
            result_df.to_csv(out_file, index=None)

        # Log estimated time left
        seconds_left = (len(batch_urls)-i)*(t2-t1)
        batch_time = str(datetime.timedelta(seconds=round(t2-t1)))
        time_left = str(datetime.timedelta(seconds=round(seconds_left)))

        print(
            f"Progress: {(i+1)/len(batch_urls)*100:.0f}% ({i+1}/{len(batch_urls):,d}). " +
            f"URLs: {i*batch_size}-{(i+1)*batch_size-1}. " +
            f"Batch time: {batch_time}. Time left: {time_left}.  ", end="\r")

    log.info("\nDone.")

    return result


def main():
    # parser = argparse.ArgumentParser(description='Scrapes Portal del la Reserca')
    # parser.add_argument('items', type=str, help='Items to be scraped: author_links, paper_links, group_links, project_links, author, paper, group, project')
    # parser.add_argument('-b', '--batch-size', default=20, type=int, help='Number of URLs to be scraped in each batch')
    # parser.add_argument('-o', '--out', type=str, help='Output file')
    # args = parser.parse_args()

    # if args.out is None:
    #     date_today = get_date()
    #     args.out = f'{date_today}_{args.items}.csv'

    # asyncio.run(scrape(items=args.items, batch_size=args.batch_size, out_file=args.out))

    items = snakemake.rule
    batch_size = snakemake.params.batch_size
    out_file = snakemake.output[0]

    print(items)
    print(batch_size)
    print(out_file)

    asyncio.run(scrape(items=items, batch_size=batch_size, out_file=out_file))


if __name__ == "__main__":
    main()

