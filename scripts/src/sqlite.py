#!/user/bin/env python
"""
SQLite Helper functions

Main functions:
    insert_urls()
    insert_authors()
    insert_papers()
    insert_projects()
    insert_groups()
"""

import os
import sys
import sqlite3
import logging

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))

log = logging.getLogger(__name__)


def insert_urls(urls, items, date):
    rows = []
    for url in urls:
        row = {}
        row['items'] = items
        row['url_stem'] = url
        row['current'] = 1
        row['url_scraped'] = 0
        row['date_created'] = date
        rows.append(row)

    conn = sqlite3.connect('recerca.db')
    with conn:
        # If new url: Insert.
        # If url exists (collides with unique url_stem constraint)
        #   update date_created and set current=1
        # There is no need to check for data updates since url_stem is the data itself.
        conn.executemany("""
            INSERT INTO urls
                ( items, url_stem, date_created, current, url_scraped)
            VALUES
                (:items,:url_stem,:date_created,:current,:url_scraped)
            ON CONFLICT DO UPDATE SET current=1,date_created=:date_created
            """, rows)
    conn.close()


def insert_papers(papers, date):
    for paper in papers:
        if_empty_then_None = [
            'url', 'url_stem', 'date', 'publisher', 'title', 'type', 'author',
            'sourceid', 'sourceref', 'orcids', 'citation', 'issn', 'published_in',
            'doi', 'isbn', 'uri', 'status_code', 'status_description']
        for key in if_empty_then_None:
            paper[key] = paper.get(key)
        list_to_string = ['orcids']
        for key in list_to_string:
            paper[key] = str(paper.get(key, '')) or None
        paper['current'] = 1
        paper['date_created'] = date

    conn = sqlite3.connect('recerca.db')
    try:
        with conn:
            # If url exists set current=0. If data has changed, current=0 will remain.
            conn.executemany("""UPDATE papers
                            SET current = 0
                            WHERE url_stem = :url_stem""", papers)

            # If new url: Insert. If url exists but data has changed (collides with unique constraint)
            #   update date_created and set current=1
            conn.executemany("""INSERT INTO papers (
                                url, url_stem, date, publisher, title, type, author,
                                sourceid, sourceref, orcids, citation, issn, published_in,
                                doi, isbn, uri, status_code, status_description,
                                date_created, current )
                            VALUES (
                                :url,:url_stem,:date,:publisher,:title,:type,:author,
                                :sourceid,:sourceref,:orcids,:citation,:issn,:published_in,
                                :doi,:isbn,:uri,:status_code,:status_description,
                                :date_created,:current )
                            ON CONFLICT DO UPDATE SET current=1,date_created=:date_created
                            """, papers)
    except Exception as e:
        print(e)
        import pdb; pdb.set_trace()
    conn.close()


def insert_authors(authors, date):
    for author in authors:
        if_empty_then_None = ['id', 'url', 'label', 'status_description', 'institution_group']
        for key in if_empty_then_None:
            author[key] = author.get(key)
        list_to_string = ['department', 'institution', 'institution_2', 'projects', 'groups']
        for key in list_to_string:
            author[key] = str(author.get(key, '')) or None
        author['current'] = 1
        author['date_created'] = date

    conn = sqlite3.connect('recerca.db')
    with conn:
        # If url exists set current=0. If data has changed, current=0 will remain.
        conn.executemany("""UPDATE authors
                        SET current = 0
                        WHERE url = :url""", authors)
        # If new url: Insert. If url exists but data has changed (collides with unique constraint) 
        #   update date_created and set current=1
        conn.executemany("""INSERT INTO authors (
                                id, url, label, department, institution, institution_2, projects,
                                groups, status_description, institution_group,
                                date_created, current )
                            VALUES (
                                :id,:url,:label,:department,:institution,:institution_2,:projects,
                                :groups,:status_description,:institution_group,
                                :date_created,:current )
                        ON CONFLICT DO UPDATE SET current=1,date_created=:date_created
                        """, authors)
    conn.close()


def insert_groups(groups, date):
    for group in groups:
        group['name'] = group.get('name')
        group['acronym'] = group.get('acronym')
        group['institution'] = group.get('institution')
        group['group_url'] = group.get('group_url')
        group['sgr_code'] = group.get('sgr_code')
        group['url'] = group.get('url')
        group['url_stem'] = group.get('url_stem')
        group['date_created'] = date
        group['current'] = 1
        list_to_string = ['principal_names', 'principal_ids', 'researcher_names', 'researcher_ids']
        for var in list_to_string:
            group[var] = str(group.get(var, '')) or None

    conn = sqlite3.connect('recerca.db')
    with conn:
        # If url exists set current=0. If data has changed, current=0 will remain.
        conn.executemany("""UPDATE groups
                        SET current = 0
                        WHERE url = :url""", groups)

        # If new url: Insert. If url exists but data has changed (collides with unique constraint) 
        #   update date_created and set current=1
        conn.executemany("""INSERT INTO groups (
                                name, acronym, institution, group_url, sgr_code,
                                principal_names, principal_ids, researcher_names, researcher_ids,
                                url, url_stem, date_created, current )
                            VALUES (
                                :name,:acronym,:institution,:group_url,:sgr_code,
                                :principal_names,:principal_ids,:researcher_names,:researcher_ids,
                                :url,:url_stem,:date_created,:current )
                        ON CONFLICT DO UPDATE SET current=1,date_created=:date_created
                        """, groups)

    conn.close()


def insert_projects(projects, date):
    for project in projects:
        project['current'] = 1
        project['date_created'] = date
        fill_keys = ['principal_names', 'principal_ids', 'researcher_names', 'researcher_ids']
        for keys in fill_keys:
            project[keys]= str(project.get(keys, '')) or None

    conn = sqlite3.connect('recerca.db')
    with conn:
        # If url exists set current=0. If data has changed, current=0 will remain.
        conn.executemany("""UPDATE projects
                        SET current = 0
                        WHERE url_stem = :url_stem""", projects)

        # If new url: Insert. If url exists but data has changed (collides with unique constraint)
        #   update date_created and set current=1
        conn.executemany("""INSERT INTO projects (
                                title, official_code, start_date, end_date, institution,
                                principal_names, principal_ids, researcher_names, researcher_ids,
                                url, url_stem, date_created, current )
                            VALUES (
                                :title,:official_code,:start_date,:end_date,:institution,
                                :principal_names,:principal_ids,:researcher_names,:researcher_ids,
                                :url,:url_stem,:date_created,:current )
                        ON CONFLICT DO UPDATE SET current=1,date_created=:date_created
                        """, projects)

    conn.close()
