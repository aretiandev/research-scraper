#!/usr/bin/env python
# coding: utf-8

# # Scraping Portal de la Reserca
# 
# This notebook asynchrnously scrapes information in Portal de la Reserca. It can download the following items:
# 
# - Links to author portals
# - Author information from the author portals.

# # Import modules

# In[ ]:


import pandas as pd
from requests_html import HTMLSession, AsyncHTMLSession
import time
import asyncio
import ast


# # Helper functions

# ## General helper functions

# In[ ]:


def get_max_pages(url):
    """
    Get max pages from pagination box in footer.
    """
    print("Retrieving number of URLs to scrape:", end=" ", flush=True)
    session = HTMLSession()
    r = session.get(url)
    pagination_items = r.html.find('div.discovery-result-pagination ul.pagination li')
    max_pages_str = pagination_items[-2].text.split("\n")[0].strip().replace('.','').replace(',','')
    max_pages = int(max_pages_str)
    print(f"{max_pages:,d}.")
    return max_pages


async def retry_url(s, url, attempts=10):
    """Retry fetching URL a given number of times."""
    for _ in range(attempts):
        try:
            r = await s.get(url)
            break
        except:
            time.sleep(1)
            pass
    else:
        r = {}
    return r


# ## Helper for author pages

# In[ ]:


async def scrape_author(s, url, item='name', attempts=10):
    if item == 'name':
        r = await retry_url(s, url, attempts)
        selector = 'div#fullNameDiv span'
        scraped_objects = r.html.find(selector)
        try:
            result = scraped_objects[0].text
        except IndexError:
            result = None
    
    elif item == 'id':
        r = await retry_url(s, url, attempts)
        selector = 'div#orcidDiv a span'
        scraped_objects = r.html.find(selector)
        try:
            result = scraped_objects[0].text
        except IndexError:
            result = None
    
    elif item == 'affiliation':
        url = url + '/researcherdepartaments.html?onlytab=true'
        r = await retry_url(s, url, attempts)
        selector = 'table.table tr td'
        scraped_objects = r.html.find(selector)
        result = {}
        try:
            result['department'] = scraped_objects[0].text
            result['institution'] = scraped_objects[1].text
        except:
            pass
    
    elif item == 'projects':
        url = url + '/publicresearcherprojects.html?onlytab=true'
        r = await retry_url(s, url, attempts)
        selector = 'table.table tr'
        scraped_objects = r.html.find(selector)
        project_list = []
        for i in range(1,len(scraped_objects)):
            project = scraped_objects[i].find('td a')[0].attrs['href']
            project_list.append(project)
        result = project_list
    
    elif item == 'groups':
        url = url + '/orgs.html?onlytab=true'
        r = await retry_url(s, url, attempts)
        selector = 'table.table tr'
        scraped_objects = r.html.find(selector)
        group_list = []
        for i in range(1,len(scraped_objects)):
            group = scraped_objects[i].find('td a')[0].attrs['href']
            group_list.append(group)
        result = group_list
    
    return result


# ## Helper for project and group pages

# In[ ]:


async def scrape_project(s, url, tab='information', attempts=10):
    if tab == 'information':
        r = await retry_url(s, url, attempts)
        
        selector = 'div#collapseOneprimarydata'
        table = r.html.find(selector, first=True)
        
        project = {}
        
        attributes = {
            '#titleDiv' : 'title', 
            '#oficialcodeDiv': 'official code',
            '#programDiv': 'program',
            '#startdateDiv': 'start date',
            '#expdateDiv': 'end date',
            '#universityDiv': 'institution',
        }
        
        for selector in attributes.keys():
            try:
                attribute_value = table.find(selector, first=True).text
            except AttributeError:
                continue
            
            attribute       = attributes[selector]
            project[attribute] = attribute_value
            
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
            except:
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
                except:
                    pass
                
            project['researcher names'] = res_names
            project['researcher ids'] = res_ids
        except:
            pass
        
                
    return project


# In[ ]:


async def scrape_group(s, url, tab='information', attempts=10):
    if tab == 'information':
        r = await retry_url(s, url, attempts)
        
        selector = 'div#collapseOneorgcard'
        table = r.html.find(selector, first=True)
        
        group = {}
        
        attributes = {
            '#nameDiv' : 'name', 
            '#acronymDiv': 'acronym',
            '#sgrDiv': 'sgr',
            '#urlDiv': 'url',
            '#universityDiv': 'institution',
            # '#startdateDiv': 'start date',
            # '#expdateDiv': 'end date',
        }
        
        for selector in attributes.keys():
            try:
                attribute_value = table.find(selector, first=True).text
            except AttributeError:
                continue
            
            attribute       = attributes[selector]
            group[attribute] = attribute_value
            
            
    elif tab == 'researchers':
        r = await retry_url(s, url, attempts)
        
        # Get next tab url
        selector = 'div#tabs ul li'
        next_tab_url = r.html.find(selector)[1].find('a', first=True).attrs['href']
        url_root = 'https://portalrecerca.csuc.cat'
        next_tab_url = url_root + next_tab_url
                                                                
        r = await retry_url(s, next_tab_url, attempts)
        
        selector = 'table.table'
        tables = r.html.find(selector)
        
        group = {}
        
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
            except:
                pass
            
        group['principal names'] = principal_names
        group['principal ids'] = principal_ids
        
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
                except:
                    pass
                
            group['researcher names'] = res_names
            group['researcher ids'] = res_ids
        except:
            pass
        
                
    return group


# ## Scrape single URL

# In[ ]:


async def scrape_url(s, url, items='authors', attempts=10):
    """
    Async scrape URL. 
    Items = ['authors', 'papers', 'paper_links', 'author_links']
    """
    
    if items == 'authors':
        result = await asyncio.gather(
                scrape_author(s, url, 'name'),
                scrape_author(s, url, 'id'),
                scrape_author(s, url, 'affiliation'),
                scrape_author(s, url, 'projects'),
                scrape_author(s, url, 'groups')
            )
        
        author = {}
        author['name'] = result[0]
        author['id'] = result[1]
        try:
            author['department'] = result[2]['department']
        except KeyError:
            pass
        try:
            author['institution'] = result[2]['institution']
        except KeyError:
            pass
        author['projects'] = result[3]
        author['groups'] = result[4]

        return author
        
        
    elif items == 'projects':    
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
    
    
    elif items == 'groups':    
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
    
        
    elif items == 'papers':
        paper = {}
        paper['url']         = url
        paper['url_id']      = url[31:].split('?')[0]
        
        not_found_msg = 'No ha estat possible trobar el que esteu buscant'
        
        r = await retry_url(s, url, attempts)
        
        try:
            title = r.html.find('title', first=True)
            table = r.html.find('table.table', first=True)
            rows = table.find('tr')
        except:
            paper['status code'] = r.status_code
            try:
                if not_found_msg in title.text:
                    paper['status description'] = 'Title: not found'
            except:
                pass
            
            return paper
            
        for row in rows:
            # Get columns in row
            columns = row.find('td')
            
            # Skip if empty row
            if len(columns) == 0:
                continue
            
            attributes = {
                'dc.contributor.authors' : 'authors', 
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
        except:
            paper['status code'] = r.status_code
            try:
                if not_found_msg in title.text:
                    paper['status description'] = 'Title: not found'
            except:
                pass
            
            return paper
            
        paper['orcids'] = author_hrefs
        
        return paper
            
        
    else: # paper_links, author_links, project_links or group_links
        r = await retry_url(s, url, attempts)
        
        try:
            table = r.html.find('div.panel.panel-info table.table', first=True)
            rows = table.find('tr')
            
            result = []

            for row in rows:
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
            
        except:
            result = []
            
        return result


# ## Main Entrypoint

# In[ ]:


async def scrape(
    items='authors', 
    urls=None, 
    start_pos=0, 
    n_pages=None, 
    batch_size=None, 
    out_file=None):
    """
    Main entry function to scrape Portal de la Reserca.
    Options:
        urls: list of urls. Not needed if items in ['author_links', 'paper_links'].
        items: [authors, papers, author_links, paper_links]
        n_pages: max pages to scrape.
        start_pos: starting position.
        batch_size: batch size.
        out_file: output file .
    """
    
    print(f"Scraping {items} from Portal de la Reserca.")
    if out_file:
        print(f"Saving results to {out_file}.")
    
    
    # Get list of URLs to scrape hyperlinks
    if items in ['author_links', 'paper_links', 'project_links', 'group_links']:
        
        if not urls:
            url_template =                                                        'https://portalrecerca.csuc.cat/simple-search?' +                 'query='                                        +                 '&location={location}'                          +                 '&filter_field_1={filter_field_1}'              +                     '&filter_type_1={filter_type_1}'            +                     '&filter_value_1={filter_value_1}'          +                 '&filter_field_2={filter_field_2}'              +                     '&filter_type_2={filter_type_2}'            +                     '&filter_value_2={filter_value_2}'          +                 '&sort_by={sort_by}'                            +                     '&order={order}'                            +                 '&rpp={rpp}'                                    +                 '&etal=0'                                       +                 '&start='
            
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
    
    print(f"Scraping {items} from {len(urls)-start_pos:,d} URLs in {len(batch_urls):,d} batches of {batch_size:,d}, starting at {start_pos:,d}.")

    if out_file:
        result_df = pd.DataFrame()

    result = []
    for i, batch in enumerate(batch_urls):

        s = AsyncHTMLSession()

        t1 = time.perf_counter()
        tasks = (scrape_url(s, url, items=items) for url in batch)
        batch_result = await asyncio.gather(*tasks)
        t2 = time.perf_counter()

        # Flatten result
        if items in ['author_links', 'paper_links', 'project_links', 'group_links']:
            batch_result = [i for sublist in batch_result for i in sublist]

        result.extend(batch_result)

        if out_file:
            result_df = result_df.append(batch_result, ignore_index=True)
            result_df.to_csv(out_file, index=None)

        # Print estimated time left
        seconds_left = (len(batch_urls)-i)*(t2-t1)
        m, s = divmod(seconds_left, 60)
        h, m = divmod(m, 60)
        
        print(
            f"Progress: {(i+1)/len(batch_urls)*100:.0f}% ({i+1}/{len(batch_urls):,d}). URLs: {i*batch_size}-{(i+1)*batch_size-1}. " +
            f"Batch time: {t2-t1:.2f}s. Time left: {h:.0f}h{m:.0f}m{s:.0f}s.", end="\r")
        
    print("\nDone.")
        
    return result


# # Run Scraper

# ## Nodes: Scrape authors

# ### Get links to author pages

# In[ ]:


items = 'author_links'
batch_size = 10
out_file = './data/author_urls_test.csv'
# author_urls = await scrape(items=items, out_file=out_file)
author_urls = await scrape(items=items, batch_size=batch_size, out_file=out_file)


# ### Scrape author pages

# In[ ]:


# Build urls
author_urls = pd.read_csv('./data/author_urls.csv')
author_urls = list(author_urls['author_urls'])
url_root = 'https://portalrecerca.csuc.cat'
urls = [url_root + url for url in author_urls]

# Get author data in batch
items = 'authors'
batch_size = 10
out_file = './data/nodes_test.csv'

author_data = await scrape(items=items, urls=urls, batch_size=batch_size, out_file=out_file)


# ## Edgelist: scrape publications

# ### Get links to papers

# In[ ]:


items = 'paper_links'
batch_size = 10
out_file = './data/paper_links_test.csv'

paper_urls = await scrape(items=items, batch_size=batch_size, out_file=out_file)


# ### Get coauthors in paper pages

# In[ ]:


# Build urls
paper_urls = pd.read_csv('./data/paper_links.csv')
paper_urls = list(paper_urls['0'])
url_root = 'https://portalrecerca.csuc.cat'
urls = [url_root + url + '?mode=full' for url in paper_urls]

# Run in batch
items = 'papers'
batch_size = 10
start_pos = 56040 # old starting position with problems
# start_pos = 0
out_file = './data/papers_test.csv'

papers = await scrape(items=items, urls=urls, start_pos=start_pos, batch_size=batch_size, out_file=out_file)


# ## Nodes: Scrape projects

# ### Get links to projects from Portal de la Reserca

# In[ ]:


items = 'project_links'
batch_size = 10
out_file = './data/project_links_test.csv'
project_urls = await scrape(items=items, batch_size=batch_size, out_file=out_file)


# ### Alternative: get links to projects from nodes

# In[ ]:


# nodes_df = pd.read_csv("./data/nodes.csv")
# projects = set()
# for projects_string in nodes_df['projects']:
#     projects_list = ast.literal_eval(projects_string)
#     projects.update(projects_list)


# ### Build URLS

# In[ ]:


project_urls = pd.read_csv('./data/project_links.csv')
project_urls = list(project_urls['0'])
url_root = 'https://portalrecerca.csuc.cat'
urls = [url_root + url for url in project_urls]


# ### Scrape projects

# In[ ]:


items = 'projects'
batch_size = 10
out_file = './data/projects_test.csv'
projects = await scrape(items=items, urls=urls, batch_size=batch_size, out_file=out_file)


# ## Node: Scrape groups

# ### Get links to groups from Portal de la Reserca

# In[ ]:


items = 'group_links'
batch_size = 10
out_file = './data/group_links_test.csv'
group_urls = await scrape(items=items, batch_size=batch_size, out_file=out_file)


# ### Alternative: get links to groups from nodes

# In[ ]:


# nodes_df = pd.read_csv("./data/nodes.csv")
# groups = set()
# for groups_string in nodes_df['groups']:
#     groups_list = ast.literal_eval(groups_string)
#     groups.update(groups_list)


# ### Build URLS

# In[ ]:


group_urls = pd.read_csv('./data/group_links.csv')
group_urls = list(group_urls['0'])
url_root = 'https://portalrecerca.csuc.cat'
urls = [url_root + url for url in group_urls]


# ###  Scrape groups

# In[ ]:


items = 'groups'
batch_size = 10
out_file = './data/groups_test.csv'
groups = await scrape(items=items, urls=urls, batch_size=batch_size, out_file=out_file)

