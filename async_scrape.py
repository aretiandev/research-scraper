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


# # Helper functions

# ## Helper Functions for links

# In[ ]:


async def scrape_url(s, url, items='author_links'):
    """
    Async scrape URL. 
    Items = ['paper_links', 'author_links', 'papers', 'authors']
    """
    # print(f"Scraping url: {url}")
    r = await s.get(url)
    table = r.html.find('div.panel.panel-info table.table', first=True)
    rows = table.find('tr')
    
    result_list = []
    
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
            
        elif items == 'author_links':
            # Get paper link
            author_link = columns[0].find('a')[0].attrs['href']
            scrape_item = author_link
            
        elif items == 'papers':
            # Get paper data
            paper_date    = columns[0].text
            paper_title   = columns[1].text
            paper_authors = columns[2].text
            paper_type    = columns[3].text

            paper = {}
            paper['date'] = paper_date
            paper['title'] = paper_title
            paper['type'] = paper_type

            paper_authors_list= paper_authors.split(';')
            for i, author in enumerate(paper_authors_list):
                paper[f"author_{i}"] = author

            scrape_item = paper

        elif items == 'authors':
            # Get author data
            author_name = columns[0].text
            author_last = author_name.split(',')[0]

            try:
                author_first = author_name.split(',')[1]
            except IndexError:  # If there is no comma in the name
                author_first = ''

            author_inst = columns[1].text

            author_dict = {
                'Last Name': author_last,
                'First Name': author_first,
                'Institution': author_inst,
                }

            scrape_item = author_dict
            
        # Append to paper links list
        result_list.append(scrape_item)
        
    return result_list


async def main(urls, items='author_links'):
    """
    Async main loop. Run withn await main(urls) in Jupyter Notebook.
    """
    s = AsyncHTMLSession()
    tasks = (scrape_url(s, url, items) for url in urls)
    return await asyncio.gather(*tasks)

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

async def scrape(items='author_links', n_pages=None):
    """
    Scrape Portal de la Reserca.
    Options:
        items   = ['paper_links', 'author_links', 'papers', 'authors']
        n_pages = Number of pages to scrape.
    """
    print(f"Scraping {items} from Portal de la Reserca.")
    
    # Url root for authors
    url_root =         'https://portalrecerca.csuc.cat/simple-search?' +         'query='                                        +         '&location=crisrp'                              +         '&filter_field_1=resourcetype'                  +             '&filter_type_1=equals'                     +             '&filter_value_1=Researchers'               +         '&sort_by=crisrp.fullName_sort'                 +             '&order=asc'                                +         '&rpp=300'                                      +         '&etal=0'                                       +         '&start='

    if not n_pages:
        print("Calculating number of pages to scrape.")
        max_pages = get_max_pages(url_root + '0')
        n_pages = max_pages

    urls = [url_root + str(page*300) for page in range(n_pages)]

    items_to_scrape = 'author_links'

    print(f"Scraping {len(urls)} URLs...")
    t1 = time.perf_counter()
    result = await main(urls, items=items_to_scrape)
    t2 = time.perf_counter()
    
    # Gather all results into single list
    full_list = [href for sublist in result for href in sublist]
    
    print(f"Scraped {len(full_list)} items in {t2-t1:.2f} seconds.")
    
    return full_list
   


# ## Helper functions for Author pages

# In[ ]:


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


async def scrape_author(s, url):
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
        
    
async def scrape_authors(urls):
    s = AsyncHTMLSession()
    tasks = (scrape_author(s, url) for url in urls)
    
    return await asyncio.gather(*tasks)


# In[ ]:


async def scrape_authors_batch(urls, start_pos=0, batch_size=100, out_file=None):
    """Scrape author pages in batches."""
    
    batch_urls = [urls[i:i+batch_size] for i in range(0, len(urls), batch_size)]
    
    print(f"Scraping {len(urls)-start_pos} author pages in {len(batch_urls)} batches of {batch_size}.")
    if out_file:
        print(f"Saving results to {out_file}.")
    
    if out_file:
        result_df = pd.DataFrame(columns=['name', 'id', 'department', 'institution', 'projects', 'groups'])
        
    result = []
    for i, batch in enumerate(batch_urls):
        print(f"Scraping batch: {i+1}/{len(batch_urls)}. Authors: {i*batch_size}-{(i+1)*batch_size-1}.", end="\r")
        
        t1 = time.perf_counter()
        author_result = await scrape_authors(batch)
        t2 = time.perf_counter()
        
        # Print estimated time left
        seconds_left = (len(batch_urls)-i)*(t2-t1)
        m, s = divmod(seconds_left, 60)
        h, m = divmod(m, 60)
        
        print(f"Last batch: {t2-t1:.2f} seconds. Estimated time left: {h:.0f}h{m:.0f}m{s:.0f}s.", end=" ")
        
        result.extend(author_result)
        
        if out_file:
            result_df = result_df.append(author_result, ignore_index=True)
            result_df.to_csv(out_file, index=None)

    print("\nDone.")
    return result


# # Nodelist: Scrape authors

# ## Get links to author pages

# In[ ]:


# Get links to author pages (takes 1m30s)
# author_urls = asyncio.run(scrape('author_links'))

## Save author URLs

# author_urls_df = pd.DataFrame(author_urls, columns=['author_urls'])
# author_urls_df.to_csv('./data/author_urls.csv', index=False)

## Read author URLs

author_urls = pd.read_csv('./data/author_urls.csv')
author_urls = list(author_urls['author_urls'])


# ## Scrape author pages

# In[ ]:


# Build urls
url_root = 'https://portalrecerca.csuc.cat'
urls = [url_root + url for url in author_urls]

# Run in batch
batch_size=100
out_file = './data/nodelist_batch.csv'

# author_data = await scrape_authors_batch(urls, start_pos=0, batch_size=batch_size, out_file=out_file)
author_data = asyncio.run(scrape_authors_batch(urls, start_pos=0, batch_size=batch_size, out_file=out_file))
