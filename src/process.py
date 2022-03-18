#!/usr/bin/env python
# coding: utf-8
"""
Process module: collection of helpers to filter data by institution and add varibables to nodelist.
"""

import numpy as np
import pandas as pd
from ast import literal_eval


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
    os.environ['TZ'] = 'America/New_York'
    time.tzset()
    date_today = date.today().strftime("%Y%m%d")
    return date_today


def assign_to_group(row, institution_list, label):
    """Helper function for institution group assignment."""
    if bool(set(row['institution']) & set(institution_list)):
        inst_groups = row['institution'].copy()
        inst_groups.append(label)
        return inst_groups
    else:
        return row['institution_group']


def add_acronym(row):
    """Helper function for completing acronyms in institution names."""
    acronym_dict = {
        'IJC': 'Institut de Recerca contra la Leucèmia Josep Carreras',
        'IrsiCaixa': 'IrsiCaixa AIDS Research Institute' }
    for acronym, name in acronym_dict.items():
        if row['institution_2'] == name:
            row['institution'].append(acronym)
    return row['institution']


def clean_papers(date_today=None):
    """Clean scraped papers data."""

    if not date_today:
        date_today = get_date()

    # Load papers
    papers_df = pd.read_csv(f'./data/{date_today}_paper_data.csv')
    papers_df['orcids'] = papers_df['orcids'].apply(lambda x: convert_to_list(x))

    # Replace NaNs to empty lists to avoid loops from breaking
    mask = papers_df['orcids'].isna()
    papers_df.loc[mask, 'orcids'] = pd.Series([[] for _ in range(len(mask))])

    # Save
    out_file = f"data/{date_today}_paper_clean.csv"
    papers_df.to_csv(out_file, index=None)
    print(f"Saved '{out_file}'.")


def clean_authors(date_today=None):

    if not date_today:
        date_today = get_date()

    # Load authors
    authors_df = pd.read_csv(f'./data/{date_today}_author_data.csv')
    authors_df['institution'] = authors_df['institution'].apply(lambda x: convert_to_list(x))
    authors_df['projects'] = authors_df['projects'].apply(lambda x: convert_to_list(x))
    authors_df['groups'] = authors_df['groups'].apply(lambda x: convert_to_list(x))

    ## Create groups of institutions

    # Drop duplicate institutions
    authors_df['institution'] = authors_df['institution'].apply(lambda x: list(set(x)))

    # Add missing acronyms to institutions
    authors_df['institution'] = authors_df.apply(add_acronym, axis=1)

    # Create groups of institutions
    authors_df['institution_group'] = authors_df['institution']

    # Assign to UPC + CIMNE
    institution_list = ['UPC', 'CIMNE']
    label = 'UPC_CIMNE'
    authors_df['institution_group'] = authors_df.apply(
                            lambda x: assign_to_group(x, institution_list, label), axis=1)

    # Assign to IGTP
    institution_list = ['IGTP', 'IJC', 'IrsiCaixa']
    label = 'IGTP+'
    authors_df['institution_group'] = authors_df.apply(
                            lambda x: assign_to_group(x, institution_list, label), axis=1)

    # Drop duplicate institution groups
    authors_df['institution_group'] = authors_df['institution_group'].apply(
                            lambda x: list(set(x)))
    
    # Save
    out_file = f"data/{date_today}_author_clean.csv"
    authors_df.to_csv(out_file, index=None)
    print(f"Saved '{out_file}'.")


def filter_authors(institution, date_today=None):
    print(f"Processing institution group: {institution}.")

    if not date_today:
        date_today = get_date()
    
    # Load authors
    authors_df = pd.read_csv(f"data/{date_today}_author_clean.csv")

    # Extract authors from institution
    mask = authors_df['institution_group'].apply(lambda x: institution in x)
    authors_inst_df = authors_df.copy()[mask]
    
    # Calculate number of affiliations
    authors_inst_df['n_affiliations'] = authors_inst_df.copy()['institution'].apply(len)

    # Calculate if single or multiple affilations
    authors_inst_df['single_affiliation'] = 'Multiple affiliations'
    mask = authors_inst_df['n_affiliations'] == 1
    authors_inst_df.loc[mask, 'single_affiliation'] = authors_inst_df.loc[mask, 'institution'].apply(lambda x: x[0])
    
    # Add projects and groups
    print("Adding projects and groups.")
    authors_inst_df['n_projects'] = authors_inst_df['projects'].apply(len)
    authors_inst_df['n_groups']   = authors_inst_df['groups'].apply(len)

    # Save
    out_file = f'./data/{date_today}_nodes_{institution}.csv'
    authors_inst_df.to_csv(out_file, index=None)
    print(f"Saved '{out_file}'.")


def filter_papers(institution, date_today=None):
    if not date_today:
        date_today = get_date()

    # Load authors
    authors_inst_df = pd.read_csv(f'./data/{date_today}_nodes_{institution}.csv')
    # Load papers
    papers_df = pd.read_csv(f"data/{date_today}_paper_clean.csv", converters={'orcids':eval})

    # Extract papers with authors from institution
    print(f"Extracting papers of researchers from {institution}.")
    authors_inst = authors_inst_df['id'].unique() # Get list of authors
    mask_1plus = papers_df['orcids'].apply(lambda x: bool(set(x) & set(authors_inst)))
    # Debug:
#     mask_1plus_backup = mask_1plus.copy()
#     mask_1plus = mask_1plus_backup.copy()
    print("Extracting papers with more than 2 authors.")
    mask_2plus = papers_df.loc[mask_1plus, 'orcids'].apply(lambda x: len(set(x) & set(authors_inst)) > 1)
    papers_inst_1plus_df = papers_df[mask_1plus]
    papers_inst_2plus_df = papers_df[mask_1plus][mask_2plus]
    
    # Save
    out_file_1 = f'./data/{date_today}_papers_{institution}.csv'
    out_file_2 = f'./data/{date_today}_papers_{institution}_2plus.csv'
    papers_inst_1plus_df.to_csv(out_file_1, index=None)
    papers_inst_2plus_df.to_csv(out_file_2, index=None)
    print(f"Saved '{out_file_1}'.")
    print(f"Saved '{out_file_2}'.")
    print("")
