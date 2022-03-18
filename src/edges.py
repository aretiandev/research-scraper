#!/usr/bin/env python
# coding: utf-8
"""
Edges module: collection of helper functions to create edgelists.
"""

import pandas as pd
from itertools import combinations


def get_date():
    import os
    import time
    from datetime import date
    os.environ['TZ'] = 'America/New_York'
    time.tzset()
    date_today = date.today().strftime("%Y%m%d")
    return date_today
    

def create_edgelist(institution, date_today=None, save=True):
    """Create edgelist."""
    print(f"Institution: {institution}.")

    if not date_today:
        date_today = get_date()

    # Load papers with coauthors list
    print(f"{institution} - Loading papers.")
    papers_df = pd.read_csv(f'./data/{date_today}_papers_{institution}_2plus.csv', converters = {'orcids': eval})
    papers = papers_df['orcids'].copy()

    # Get unique list of authors from papers
    authors_papers = list(set(papers.sum()))
    authors_papers.sort()

    # Get list of authors from institution
    print(f"{institution} - Loading nodes.")
    authors_inst_df = pd.read_csv(f'./data/{date_today}_nodes_{institution}.csv')
    authors_inst = authors_inst_df['id']
    authors_inst = authors_inst.unique()
    authors_inst.sort()

    # Combine both
    authors_index = list(set(authors_papers) & set(authors_inst))
    authors_index.sort()
    
    # Create df to store collaborations
    print(f"{institution} - Calculting combinations of authors.")
    author_combinations = combinations(authors_index,2)
    collabs_df = pd.DataFrame(list(author_combinations), columns=['Source', 'Target'])
    collabs_df['Weight'] = 0
    collabs_df = collabs_df.set_index(['Source', 'Target'])

    # Calculate collaborations
    print(f"{institution} - Main loop: counting collaborations.")
    for i, paper in enumerate(papers):
        print(f"{institution} - Progress: {i/len(papers)*100:.0f}%. ({i:,.0f}/{len(papers):,.0f}).", end="\r")

        # Store collaboration
        paper = list(set(paper))
        paper.sort()
        author_pairs = combinations(paper, 2)
        for pair in author_pairs:
            try:
                collabs_df.loc[pair] += 1
            except:
                pass
    collabs_df = collabs_df.reset_index()
    
    if save:
        outfile = f'./data/{date_today}_edges_{institution}.csv'
        collabs_df.to_csv(outfile, index=None)
        print(f"{institution} - Done. Saved '{outfile}'.")


if __name__ == "__main__":
    institution_list = ['IGTP+', 'UPC_CIMNE', 'UB', 'UPF', 'UVic-UCC', 'UOC']
    for institution in institution_list:
        create_edgelist(institution)
