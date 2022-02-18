from requests_html import HTMLSession
import pandas as pd

# Start requests-HTML session
session = HTMLSession()

# url = 'https://portalrecerca.csuc.cat/simple-search?filtername=resourcetype&filterquery=Researchers&filtertype=equals&sort_by=crisrp.fullName_sort&order=ASC&location=crisrp'
# url_root = 'https://portalrecerca.csuc.cat/simple-search?query=&location=crisrp&filter_field_1=resourcetype&filter_type_1=equals&filter_value_1=Researchers&sort_by=crisrp.fullName_sort&order=asc&rpp=100&etal=0&start='

# url_root = 'https://portalrecerca.csuc.cat/simple-search?query=&location=publications&filter_field_1=resourcetype&filter_type_1=equals&filter_value_1=Items&filter_field_2=itemtype&filter_type_2=notequals&filter_value_2=Phd+Thesis&sort_by=dc.contributor.authors_sort&order=asc&rpp=300&etal=0&start='


# Download IGTP papers
url_root = 'https://portalrecerca.csuc.cat/simple-search?query=&location=publications&filter_field_1=resourcetype&filter_type_1=equals&filter_value_1=Items&filter_field_2=itemtype&filter_type_2=notequals&filter_value_2=Phd+Thesis&filter_field_3=location.coll&filter_type_3=equals&filter_value_3=349&sort_by=dc.contributor.authors_sort&order=asc&rpp=300&etal=0&start='

# Begin loop
page_number = 0  # Start loop here
# max_pages = 2
# max_pages = 100
# max_pages = 2240
max_pages = 19
papers_list = []

while page_number < max_pages:
    print(f"Progress: {round(page_number/(max_pages)*100)}%. Scraping page: {page_number}/{max_pages}.", end="\r")
    print('')

    # Create URL
    url = url_root + str(300*page_number)
    # Update page counter
    page_number += 1
    # Get page
    r = session.get(url)
    # Get table
    table = r.html.find('div.panel.panel-info table.table', first=True)
    # Read rows from table
    rows = table.find('tr')

    for row in rows:
        columns = row.find('td')
        # Skip if empty row
        if len(columns) == 0:
            continue
        # Get author data
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

        # try:
            # author_first = author_name.split(',')[1]
        # except Exception:  # If there is no comma in the name
            # author_first = ''
        # author_inst = columns[1].text

        # author_dict = {
        #     'Last Name': author_last,
        #     'First Name': author_first,
        #     'Institution': author_inst,
        #     }

        # Append to author list
        papers_list.append(paper)

papers_df = pd.DataFrame.from_records(papers_list)

data_folder =  '../data/network of talent/research networks/Scraping/'
papers_df.to_csv(data_folder + 'raw/papers_IGTP.csv')
print('')
print('Saved output to file papers_IGTP.csv')
