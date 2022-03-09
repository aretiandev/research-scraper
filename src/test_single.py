"""
Test example to show different results between JupyterLab and Command Line.
"""

from requests_html import HTMLSession

def test_scrape(limit=5):

    url1 = 'https://portalrecerca.csuc.cat/simple-search?query=&location=crisrp&filter_field_1=resourcetype&filter_type_1=equals&filter_value_1=Researchers&filter_field_2=&filter_type_2=&filter_value_2=&sort_by=crisrp.fullName_sort&order=asc&rpp=300&etal=0&start=5400'
    url2 = 'https://portalrecerca.csuc.cat/simple-search?query=&location=crisrp&filter_field_1=resourcetype&filter_type_1=equals&filter_value_1=Researchers&filter_field_2=&filter_type_2=&filter_value_2=&sort_by=crisrp.fullName_sort&order=asc&rpp=300&etal=0&start=5700'

    batch = [url1, url2]

    ## TEST

    s = HTMLSession()
    r = s.get(url1)

    # Render Javascript
    # r.html.render(wait=20)

    print("# Request #")
    print(f"Headers: {r.request.headers}.")
    print(f"HTML length: {len(r.text):,.0f}.")
    print("")

    print("# Parsing #")
    table = r.html.find('div.panel.panel-info table.table', first=True)
    # print(f"Table length: {len(table.text):,.0f}.")

    rows = table.find('tr')
    print(f"Number of rows: {len(rows):,.0f}.")

    result = []

    for i, row in enumerate(rows):
        # Difference #1 between notebook and script results
        if i < 200:
            continue

        print(f"Row length: {len(row.text):,.0f}.")

        if limit and i > limit:
            break

        columns = row.find('td')

        # Difference #2 between notebook and script results
        test_item = columns[0].find('a', first=True)
        print(f"Name cell html length: {len(test_item.html):,.0f}.")
        print(f"Name cell text length: {len(test_item.text):,.0f}.")

        try:
            scrape_item = test_item.attrs['href']
            print(f"ORCID obtained from cell: {scrape_item}.")
        except Exception as e:
            print(f"Found error in row: {i}")
            print(f"Row content: {row.text[:100]}")
            print(f"Exception: {e}.")
            break

        result.append(scrape_item)

    print(result)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        limit = sys.argv[1]
        if limit == 'all':
            limit = None
    else:
        limit = 10 
    print(f"Limit: {limit}.")
    test_scrape(limit=limit)
