from requests_html import HTMLSession
import pandas as pd

# Start requests-HTML session
session = HTMLSession()

# url = 'https://portalrecerca.csuc.cat/simple-search?filtername=resourcetype&filterquery=Researchers&filtertype=equals&sort_by=crisrp.fullName_sort&order=ASC&location=crisrp'
url_root = 'https://portalrecerca.csuc.cat/simple-search?query=&location=crisrp&filter_field_1=resourcetype&filter_type_1=equals&filter_value_1=Researchers&sort_by=crisrp.fullName_sort&order=asc&rpp=100&etal=0&start='

# Begin loop
page_number = 0 # Start loop here
max_pages = 190
authors_list = []

while page_number <= max_pages:
    print(f"Progress: {round(page_number/(max_pages)*100)}%. Scraping page: {page_number}/{max_pages}.", end="\r")
    print('')

    # Create URL
    url = url_root + str(100*page_number)
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
        author_name = columns[0].text
        author_last = author_name.split(',')[0]

        try:
            author_first = author_name.split(',')[1]
        except Exception:  # If there is no comma in the name
            author_first = ''

        author_inst = columns[1].text

        author_dict = {
            'Last Name': author_last,
            'First Name': author_first,
            'Institution': author_inst,
            }

        # Append to author list
        authors_list.append(author_dict)

authors_df = pd.DataFrame.from_records(authors_list)

authors_df.to_csv('authors.csv')
print('')
print('Saved output to file authors.csv')
