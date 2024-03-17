#!/usr/bin/env python
# coding: utf-8
"""
Process module:
Collection of helpers to filter data by institution and add variables to nodelist.
"""

import logging
import sqlite3
from ast import literal_eval

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


def insert_nodes(nodes, date, database="recerca.db"):
    for node in nodes:
        node["current"] = 1
        node["date_created"] = date

    conn = sqlite3.connect(f"../{database}")
    with conn:
        # Set current = 0 for existing records
        conn.executemany(
            """UPDATE nodes
                        SET current = 0
                        WHERE url = :url""",
            nodes,
        )

        conn.executemany(
            """INSERT INTO nodes VALUES (
                            :id,:label,:url,:department,:institution,:institution_2,
                            :projects,:groups,:status_description,:institution_group,:institution_table,
                            :n_affiliations,:single_affiliation,:n_projects,:n_groups,
                            :date_created,:current)
                        ON CONFLICT DO UPDATE SET current=1
                        """,
            nodes,
        )

    conn.close()


def convert_to_list(x):
    """Convert string column to list"""
    try:
        result = literal_eval(x)
    except ValueError:
        result = np.nan
    return result


def get_date():
    """Get formatted today's date."""
    import os
    import time
    from datetime import date

    os.environ["TZ"] = "America/New_York"
    time.tzset()
    date_today = date.today().strftime("%Y%m%d")
    return date_today


def assign_to_group(row, institution_list, label):
    """Adds supplied label to list of institutions."""
    if bool(set(row["institution"]) & set(institution_list)):
        inst_groups = row["institution"].copy()
        inst_groups.append(label)
        return inst_groups
    else:
        return row["institution_group"]


def add_institution(row, institution):
    """Extends list of institutions based on institution_2 column"""
    pass


def add_acronym(row):
    """Helper function for completing acronyms in institution names."""
    acronym_dict = {
        "IJC": "Institut de Recerca contra la Leucèmia Josep Carreras",
        "IrsiCaixa": "IrsiCaixa AIDS Research Institute",
        "CRAG": "Centre de Recerca en Agrigenòmica",
    }
    for acronym, name in acronym_dict.items():
        if row["institution_2"] == name:
            row["institution"].append(acronym)
    return row["institution"]


def clean_authors(input, output):
    """Clean author data

    Args:
        input (csv):  data/{date}_author_data.csv
        output (csv): data/{date}_author_clean.csv
    """
    # Load authors
    authors_df = pd.read_csv(input)
    authors_df["institution"] = authors_df["institution"].apply(
        lambda x: convert_to_list(x)
    )
    authors_df["projects"] = authors_df["projects"].apply(lambda x: convert_to_list(x))
    authors_df["groups"] = authors_df["groups"].apply(lambda x: convert_to_list(x))

    # Create groups of institutions

    # Drop duplicate institutions
    authors_df["institution"] = authors_df["institution"].apply(lambda x: list(set(x)))

    # Add missing acronyms to institutions
    authors_df["institution"] = authors_df.apply(add_acronym, axis=1)

    # Create groups of institutions
    authors_df["institution_group"] = authors_df["institution"]

    # Assign to UPC + CIMNE
    institution_list = ["UPC", "CIMNE"]
    label = "UPC_CIMNE"
    authors_df["institution_group"] = authors_df.apply(
        lambda x: assign_to_group(x, institution_list, label), axis=1
    )

    # Assign to IGTP
    institution_list = ["IGTP", "IJC", "IrsiCaixa"]
    label = "IGTP+"
    authors_df["institution_group"] = authors_df.apply(
        lambda x: assign_to_group(x, institution_list, label), axis=1
    )

    # Drop duplicate institution groups
    authors_df["institution_group"] = authors_df["institution_group"].apply(
        lambda x: list(set(x))
    )

    # TODO
    # Add label to IRSJD researchers
    # Pseudo code: if institution_2 == "Sant Joan de Deu", add IRSJD label to institution_group
    selected_institution = "Institut de Recerca Sant Joan de Déu"
    institution_label = "IRSJD"

    authors_df.loc[
        authors_df["institution_2"] == selected_institution, "institution_group"
    ].apply(lambda x: x.append(institution_label))

    # new_institution = authors_df.loc[
    #     authors_df["institution_2"] == selected_institution, "institution"
    # ].apply(lambda x: x.append(institution_label))
    # authors_df.loc[
    #     authors_df["institution_2"] == selected_institution, "institution"
    # ] = new_institution

    # print("len(new_institution")
    # print(len(new_institution))
    # print(new_institution.head())

    # Debug
    # print(authors_df["institution"].head())
    # print(authors_df["institution"].head()[0])

    # def identify_institution(x):
    #     if x == []:
    #         return 0
    #     try:
    #         if "IRSJD" in x:
    #             return 1
    #     except Exception:
    #         print(x)
    #     return 0

    # test_length = authors_df["institution_group"].apply(identify_institution).sum()
    # # 20230129 2pm: Problem description: some of the rows are empty lists therefore raises an exception that Nonetype is not iterable. solution: write a function that skips the empty lists
    # print(f"authors in IRSJD: {test_length}")
    # if test_length == 0:
    #     raise Exception("No authors in IRSJD")

    # Save
    authors_df.to_csv(output, index=None)
    log.info(f"Saved '{output}'.")


def clean_papers(input, output):
    """Clean scraped papers data.

    Args:
        input (csv):  data/{date}_paper_data.csv
        output (csv): data/{date}_paper_clean.csv
    """
    # Load papers
    papers_df = pd.read_csv(input)
    papers_df["orcids"] = papers_df["orcids"].apply(lambda x: convert_to_list(x))

    # Replace NaNs to empty lists to avoid loops from breaking
    mask = papers_df["orcids"].isna()
    papers_df.loc[mask, "orcids"] = pd.Series([[] for _ in range(len(mask))])

    # Save
    papers_df.to_csv(output, index=None)
    log.info(f"Saved '{output}'.")


def clean(items, input, output):
    if items == "author":
        clean_authors(input, output)
    if items == "paper":
        clean_papers(input, output)


def filter_authors(input, output, institution, out_sql=False, database="recerca.db"):
    """
    Filter authors by institution.

    Args:
        input (csv):  data/{date}/{date}_author_clean.csv
        output (csv): data/{date}_nodes_{institution}.csv
        institution (str): name of institution to filter authors.
    """
    log.info(f"Processing institution group: {institution}.")

    # Load authors
    authors_df = pd.read_csv(input)

    # Extract authors from institution
    mask = authors_df["institution_group"].apply(lambda x: institution in x)
    authors_inst_df = authors_df.copy()[mask]

    # Calculate number of affiliations
    authors_inst_df["n_affiliations"] = authors_inst_df.copy()["institution"].apply(len)

    # Calculate if single or multiple affilations
    authors_inst_df["single_affiliation"] = "Multiple affiliations"
    mask = authors_inst_df["n_affiliations"] == 1
    authors_inst_df.loc[mask, "single_affiliation"] = authors_inst_df.loc[
        mask, "institution"
    ].apply(lambda x: x[0])

    # Add projects and groups
    log.info("Adding projects and groups.")
    authors_inst_df["n_projects"] = authors_inst_df["projects"].apply(len)
    authors_inst_df["n_groups"] = authors_inst_df["groups"].apply(len)

    # Debug
    if len(authors_inst_df) == 0:
        print("ERROR: Empty dataframe")
        print(f"Institution: {institution}")
        print(authors_inst_df)
        print(authors_df["institution_group"].apply(lambda x: institution in x).sum())
        print(authors_df["institution_group"].head())
        print(input)
        raise Exception("Empty DataFrame.")

    # Save
    authors_inst_df.to_csv(output, index=None)
    log.info(f"Saved '{output}'.")

    # Save to SQLite
    if out_sql:
        insert_nodes(authors_inst_df.to_dict("records"), database=database)


def add_nodes_stats(
    input_authors,
    input_papers,
    input_groups,
    output,
    institution,
    out_sql=False,
    database="recerca.db",
):
    """
    Filter authors by institution.

    Args:
        input (csv):  data/{date}/{date}_author_clean.csv
                      data/{date}/{date}_papers_{inst}.csv
                      data/{date}/{date}_group_data.csv
        output (csv): data/{date}_nodes_{institution}.csv
        institution (str): name of institution to filter authors.
    """
    log.info(f"Processing institution group: {institution}.")

    # Load authors
    authors_inst_df = pd.read_csv(input_authors, converters={"groups": eval})

    # Add publication stats
    authors_inst_df = authors_inst_df.set_index("id")  # This mutes pandas warning

    # Load papers
    papers_inst_df = pd.read_csv(input_papers, converters={"orcids": eval})
    # Get papers column
    papers = papers_inst_df[["orcids", "type"]].copy()
    papers = papers.reset_index(drop=True)

    authors_inst_df["n_publications"] = 0
    authors_inst_df["n_articles"] = 0
    authors_inst_df["n_chapters"] = 0
    authors_inst_df["n_books"] = 0
    authors_inst_df["n_other"] = 0

    paper_types = {
        "Journal Article": "n_articles",
        "Chapter in Book": "n_chapters",
        "Book": "n_books",
    }

    def add_publication_stats(paper):
        for orcid in paper["orcids"]:
            # Add publication
            try:
                authors_inst_df.loc[orcid, "n_publications"] += 1
            except KeyError:  # author is not in the institution
                continue

            # Assign type
            assigned_type = False
            for paper_type, column in paper_types.items():
                if paper["type"] == paper_type:
                    authors_inst_df.loc[orcid, column] += 1
                    assigned_type = True
                    break
            if not assigned_type:
                authors_inst_df.loc[orcid, "n_other"] += 1

    log.info(
        "Adding publications number and types to nodelist. This might take a while..."
    )
    papers.apply(lambda x: add_publication_stats(x), axis=1)

    # Add groups names
    df_groups = pd.read_csv(input_groups)
    df_groups["url_id"] = df_groups["url"].str[30:]

    def get_group_names(lst, df_groups):
        groups_names = []
        for group_id in lst:
            group_name = df_groups.loc[df_groups["url_id"] == group_id, "name"].values[
                0
            ]
            groups_names.append(group_name)
        return groups_names

    # Merge
    authors_inst_df["groups_names"] = authors_inst_df["groups"].apply(
        lambda x: get_group_names(x, df_groups)
    )

    # Retrieve index
    authors_inst_df = authors_inst_df.reset_index("id")

    # Save
    log.info(f"Saved: {output}")
    authors_inst_df.to_csv(output, index=None)

    # Save to SQLite
    if out_sql:
        insert_nodes(authors_inst_df.to_dict("records"), database=database)


def filter_papers(
    input_authors, input_papers, output_papers, output_papers_2plus, institution
):
    """
    Filter authors by institution.

    Args:
        input_authors (csv): data/{date}_nodes_{institution}.csv
        input_papers (csv):  data/{date}_paper_clean.csv
        output_papers (csv):       data/{date}_papers_{institution}.csv
        output_papers_2plus (csv): data/{date}_papers_{institution}_2plus.csv
        institution (str): name of institution to filter papers.
    """
    # Load authors
    log.info(f"Loading '{input_authors}'.")
    authors_inst_df = pd.read_csv(input_authors)
    # Debug
    print(f"Getting authors from {institution}")
    print(authors_inst_df.head())
    print("From original csv:")
    print(input_authors)
    # Load papers
    log.info(f"Loading '{input_papers}'.")
    papers_df = pd.read_csv(input_papers, converters={"orcids": eval})

    # Extract papers with authors from institution
    log.info(f"Extracting papers of researchers from {institution}.")
    authors_inst = authors_inst_df["id"].unique()  # Get list of authors
    mask_1plus = papers_df["orcids"].apply(lambda x: bool(set(x) & set(authors_inst)))
    log.info("Extracting papers with more than 2 authors.")
    mask_2plus = papers_df.loc[mask_1plus, "orcids"].apply(
        lambda x: len(set(x) & set(authors_inst)) > 1
    )
    papers_inst_1plus_df = papers_df[mask_1plus]
    papers_inst_2plus_df = papers_df[mask_1plus][mask_2plus]

    # Save
    papers_inst_1plus_df.to_csv(output_papers, index=None)
    papers_inst_2plus_df.to_csv(output_papers_2plus, index=None)
    log.info(f"Saved '{output_papers}'.")
    log.info(f"Saved '{output_papers_2plus}'.")


def main():
    pass


if __name__ == "__main__":
    main()
