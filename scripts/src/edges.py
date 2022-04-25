#!/usr/bin/env python
# coding: utf-8
"""
Edges module: collection of helper functions to create edgelists.
"""

import pandas as pd
from itertools import combinations
import logging

log = logging.getLogger(__name__)


def get_date():
    import os
    import time
    from datetime import date
    os.environ['TZ'] = 'America/New_York'
    time.tzset()
    date_today = date.today().strftime("%Y%m%d")
    return date_today


def create_edges(input_authors, input_papers, output, institution):
    """
    Create edgelist.

    Args:
        input_authors (csv): data/{date}_nodes_{institution}.csv
        input_papers (csv):  data/{date}_papers_{institution}_2plus.csv
        output (csv):        data/{date}_edges_{institution}.csv
        institution (str):   name of institution to filter papers.
    """
    """Create edgelist."""
    log.info(f"Institution: {institution}.")
    # Load papers with coauthors list
    log.info(f"{institution} - Loading papers.")
    papers_df = pd.read_csv(input_papers, converters={'orcids': eval})
    papers = papers_df['orcids'].copy()

    # Get unique list of authors from papers
    authors_papers = list(set(papers.sum()))
    authors_papers.sort()

    # Get list of authors from institution
    log.info(f"{institution} - Loading nodes.")
    authors_inst_df = pd.read_csv(input_authors)
    authors_inst = authors_inst_df['id']
    authors_inst = authors_inst.unique()
    authors_inst = authors_inst[pd.notnull(authors_inst)]
    authors_inst.sort()

    # Combine both
    authors_index = list(set(authors_papers) & set(authors_inst))
    authors_index.sort()

    # Create df to store collaborations
    log.info(f"{institution} - Calculting combinations of authors.")
    author_combinations = combinations(authors_index, 2)
    collabs_df = pd.DataFrame(list(author_combinations), columns=['Source', 'Target'])
    collabs_df['Weight'] = 0
    collabs_df = collabs_df.set_index(['Source', 'Target'])

    # Calculate collaborations
    log.info(f"{institution} - Main loop: counting collaborations.")
    for i, paper in enumerate(papers):
        print(f"{institution} - Progress: {i/len(papers)*100:.0f}%. ({i:,.0f}/{len(papers):,.0f}).", end="\r")

        # Store collaboration
        paper = list(set(paper))
        paper.sort()
        author_pairs = combinations(paper, 2)
        for pair in author_pairs:
            try:
                collabs_df.loc[pair] += 1
            except Exception:
                pass
    collabs_df = collabs_df.reset_index()

    collabs_df = collabs_df.loc[collabs_df['Weight'] > 0]  # Drop zero weights

    collabs_df.to_csv(output, index=None)
    log.info(f"{institution} - Done. Saved '{output}'.")
