#!/usr/bin/env python
# coding: utf-8

# # Scraping Portal de la Reserca
# 
# This notebook asynchrnously scrapes information in Portal de la Reserca. It can download the following items:
# 
# - Links to author portals
# - Author information from the author portals.

# # Import modules

# In[15]:


import pandas as pd
from requests_html import HTMLSession, AsyncHTMLSession
import time
import asyncio


# # Helper functions

# ## General helper functions

# In[16]:


def get_max_pages(url):
    """
    Get max pages from pagination box in footer.
    """
    session = HTMLSession()
    r = session.get(url)
    pagination_box = r.html.find('ul.pagination.pull-right')
    pagination_items = pagination_box[0].find('li')
    max_pages = pagination_items[-2].text.replace('.','').replace(',','')
    return int(max_pages)


# ## Helper functions for individual and groups of pages

# In[29]:


async def scrape_url(s, url, items='author_links'):
    """
    Async scrape URL. 
    Items = ['paper_links', 'author_links', 'papers', 'authors']
    """
    # print(f"Scraping url: {url}")
    
    if items == 'authors':
        result = await asyncio.gather(
                scrape_author_page(s, url, 'name'),
                scrape_author_page(s, url, 'id'),
                scrape_author_page(s, url, 'institution'),
                scrape_author_page(s, url, 'projects'),
                scrape_author_page(s, url, 'groups')
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
        
    elif items == 'papers':
        while True:
            try:
                r = await s.get(url)
                table = r.html.find('table.table', first=True)
                rows = table.find('tr')
            except AttributeError:
                with open('errors.txt', 'a') as f:
                    f.write(f"Could not find row in url: {url}")
                time.sleep(1)
                continue
            break

        paper = {}
        paper['url']    = url
        paper['url_id'] = url[31:].split('?')[0]
        
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
                'dc.identifier.isbn'    : 'isbn',
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
        while True:
            try:
                r = await s.get(simple_url)
                authors = r.html.find('table.table a.authority.author')
                author_hrefs = []
                for author in authors:
                    href = author.attrs['href']
                    if href[:7] == '/orcid/':
                        href = href[7:]
                    author_hrefs.append(href)
            except AttributeError:
                with open('errors.txt', 'a') as f:
                    f.write(f"Could not find row in url: {url}")
                time.sleep(1)
                continue
            break
            
        paper['orcids'] = author_hrefs
        
        return paper
            
        
    else: # paper_links or author_links
        while True:
            try:
                r = await s.get(url)
                table = r.html.find('div.panel.panel-info table.table', first=True)
                rows = table.find('tr')
            except AttributeError:
                with open('errors.txt', 'a') as f:
                    f.write(f"Could not find row in url: {url}")
                time.sleep(1)
                continue
            break

        result = []
            
        for row in rows:
            # Get columns in row
            columns = row.find('td')
            
            # Skip if empty row
            if len(columns) == 0:
                continue

            if items == 'paper_links':
                # Get paper link
                paper_link = columns[1].find('a')[0].attrs['href']
                scrape_item = paper_link

            if items == 'author_links':
                # Get paper link
                author_link = columns[0].find('a')[0].attrs['href']
                scrape_item = author_link
                
            # Append to results_list
            result.append(scrape_item)
            
        return result


async def scrape_urls(s, urls, items='author_links'):
    """Wrapper to scrape a list of urls."""
    tasks = (scrape_url(s, url, items) for url in urls)
    return await asyncio.gather(*tasks)


# ## Helper functions for author pages

# In[18]:


async def scrape_author_tab(s, url, selector):
    while True:
        try:
            r = await s.get(url)
        except (SSLError, MaxRetryError):
            print(f"Failed request to url: {url}.")
            time.sleep(1)
            continue
        break
    result = r.html.find(selector)
    return result
    
    
async def scrape_author_page(s, url, item='name'):
    if item == 'name':
        selector = 'div#fullNameDiv span'
        result = await scrape_author_tab(s, url, selector)
        try:
            result = result[0].text
        except IndexError:
            result = None
    
    elif item == 'id':
        selector = 'div#orcidDiv a span'
        result = await scrape_author_tab(s, url, selector)
        try:
            result = result[0].text
        except IndexError:
            result = None
    
    elif item == 'institution':
        url_dep = url + '/researcherdepartaments.html?onlytab=true'
        selector = 'table.table tr td'
        institution = await scrape_author_tab(s, url_dep, selector)
        result = {}
        try:
            result['department'] = institution[0].text
            result['institution'] = institution[1].text
        except:
            pass
    
    elif item == 'projects':
        url_proj = url + '/publicresearcherprojects.html?onlytab=true'
        selector = 'table.table tr'
        projects = await scrape_author_tab(s, url_proj, selector)
        project_list = []
        for i in range(1,len(projects)):
            project = projects[i].find('td a')[0].attrs['href']
            project_list.append(project)
        result = project_list
    
    elif item == 'groups':
        url_group = url + '/orgs.html?onlytab=true'
        selector = 'table.table tr'
        groups = await scrape_author_tab(s, url_group, selector)
        group_list = []
        for i in range(1,len(groups)):
            group = groups[i].find('td a')[0].attrs['href']
            group_list.append(group)
        result = group_list
    
    return result


# ## Main Entrypoint

# In[19]:


async def scrape(urls=None, items='authors', n_pages=None, start_pos=0, batch_size=None, out_file=None):
    """
    Main entry function to scrape Portal de la Reserca.
    Options:
        urls: list of urls. If items='paper_links' this is not needed.
        items: [authors, papers, author_links, paper_links]
        start_pos: starting position
        batch_size: batch size
        out_file: output file 
    """
    
    if items == 'author_links':
        url_root =                                                        'https://portalrecerca.csuc.cat/simple-search?' +             'query='                                        +             '&location=crisrp'                              +             '&filter_field_1=resourcetype'                  +                 '&filter_type_1=equals'                     +                 '&filter_value_1=Researchers'               +             '&sort_by=crisrp.fullName_sort'                 +                 '&order=asc'                                +             '&rpp=300'                                      +             '&etal=0'                                       +             '&start='
        
        if not n_pages:
            print("Calculating number of pages to scrape.")
            max_pages = get_max_pages(url_root + '0')
            n_pages = max_pages
            
            
    elif items == 'paper_links':
        url_root =                                                        'https://portalrecerca.csuc.cat/simple-search?' +             'query='                                        +             '&location=publications'                        +             '&filter_field_1=resourcetype'                  +                 '&filter_type_1=equals'                     +                 '&filter_value_1=Items'                     +             '&filter_field_2=itemtype'                      +                 '&filter_type_2=notequals'                  +                 '&filter_value_2=Phd+Thesis'                +             '&sort_by=dc.contributor.authors_sort'          +                 '&order=asc'                                +             '&rpp=300'                                      +             '&etal=0'                                       +             '&start='

        if not n_pages:
            print("Calculating number of pages to scrape.")
            max_pages = get_max_pages(url_root + '0')
            n_pages = max_pages
        
            if not urls:
                urls = [url_root + str(page*300) for page in range(n_pages)]
        
    if not batch_size:
        print(f"Scraping {items} from Portal de la Reserca.")

        urls = [url_root + str(page*300) for page in range(n_pages)]

        s = AsyncHTMLSession()

        print(f"Scraping {len(urls)} URLs...")
        t1 = time.perf_counter()
        result = await scrape_urls(s, urls, items=items)
        t2 = time.perf_counter()

        # Gather all results into single list
        result = [href for sublist in result for href in sublist]

        print(f"Scraped {len(result)} items in {t2-t1:.2f} seconds.")
        
        if out_file:
            result_df = pd.DataFrame(result)
            result_df.to_csv(out_file, index=None)
            print(f"Saved results to {out_file}.")

   
    if batch_size: 

        if not urls:
            raise TypeError("Must provide list of urls or set items='paper_links'")

        batch_urls = [urls[i:i+batch_size] for i in range(0, len(urls), batch_size)]

        print(f"Scraping {len(urls)-start_pos} {items} in {len(batch_urls)} batches of {batch_size}.")
        if out_file:
            print(f"Saving results to {out_file}.")
            if items == 'authors':
                result_df = pd.DataFrame(columns=['name', 'id', 'department', 'institution', 'projects', 'groups'])
            elif items == 'paper_links':
                result_df = pd.DataFrame()
            elif items == 'papers':
                result_df = pd.DataFrame()

        result = []
        for i, batch in enumerate(batch_urls):
            print(f"Scraping batch: {i+1}/{len(batch_urls)}. {items}: {i*batch_size}-{(i+1)*batch_size-1}.", end="\r")
            s = AsyncHTMLSession()

            t1 = time.perf_counter()
            batch_result = await scrape_urls(s, batch, items=items)
            t2 = time.perf_counter()

            if items == 'paper_links':
                # Flatten result
                batch_result = [i for sublist in batch_result for i in sublist]

            # Print estimated time left
            seconds_left = (len(batch_urls)-i)*(t2-t1)
            m, s = divmod(seconds_left, 60)
            h, m = divmod(m, 60)

            print(f"Last batch: {t2-t1:.2f} seconds. Estimated time left: {h:.0f}h{m:.0f}m{s:.0f}s.", end=" ")

            result.extend(batch_result)

            if out_file:
                result_df = result_df.append(batch_result, ignore_index=True)
                result_df.to_csv(out_file, index=None)

        print("\nDone.")
        
    return result


# # Run Scraper

# ## Nodelist: Scrape authors

# ### Get links to author pages

# In[7]:


# Get links to author pages (takes 1m30s)
items = 'author_links'
out_file = './data/author_urls_test.csv'
author_urls = await scrape(items=items, out_file=out_file)


# ### Scrape author pages

# In[8]:


# Build urls
author_urls = pd.read_csv('./data/author_urls.csv')
author_urls = list(author_urls['author_urls'])
url_root = 'https://portalrecerca.csuc.cat'
urls = [url_root + url for url in author_urls]

# Get author data in batch
items = 'authors'
batch_size = 20
out_file = './data/nodelist_test.csv'

author_data = await scrape(urls=urls, items=items, batch_size=batch_size, out_file=out_file)


# ## Edgelist: scrape publications

# ### Get links to papers

# In[9]:


# Get links to papers in batch
items = 'paper_links'
batch_size = 20
out_file = './data/paper_links_test.csv'

paper_links = await scrape(items=items, batch_size=batch_size, out_file=out_file)


# ### Get coauthors in paper pages

# In[30]:


# Build urls
paper_urls = pd.read_csv('./data/paper_links.csv')
paper_urls = list(paper_urls['0'])
url_root = 'https://portalrecerca.csuc.cat'
urls = [url_root + url + '?mode=full' for url in paper_urls]

# Run in batch
items = 'papers'
batch_size = 20
out_file = './data/coauthors.csv'

papers = await scrape(urls=urls, items=items, batch_size=batch_size, out_file=out_file)

