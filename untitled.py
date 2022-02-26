#!/usr/bin/env python

import numpy as np
import pandas as pd
from ast import literal_eval
from itertools import combinations

# Helper functions
def convert_to_list(x):
    """Convert string column to list"""
    try:
        result = literal_eval(x)
    except ValueError:
        result = np.nan
    return result

def belongs_to_list(row, df):
    try:
        result = bool(set(row['orcids']) & set(df))
    except TypeError:
        result = False
    return result

# Main function
def calculate_collaborations(institution, papers_df, res_df, save=False, threshold=None):

    # Filter researchers
    res_df_inst = res_df.loc[res_df['institution'] == institution]
    res_inst = res_df_inst['id'].unique()

    # Convert strings to list of coauthors
    papers_df['orcids'] = papers_df['orcids'].apply(lambda x: convert_to_list(x))

    # Identify authors in institution
    institution_df = res_inst
    mask = papers_df.apply(lambda x: belongs_to_list(x, institution_df), axis=1)

    selected_df = papers_df[mask]
    
    # Test with n papers for debugging (max 2800)
    if threshold:
        selected_df = selected_df[:threshold]

    # Get papers column
    papers = selected_df['orcids'].copy()
    papers = papers.reset_index(drop=True)

    # Get unique list of authors from papers
    authors_index = list(set(papers.sum()))

    authors_index.sort()

    # Create boolean matrix with papers
    paper_bool_df = pd.DataFrame(columns=authors_index, index=range(len(papers)))

    for i, paper in enumerate(papers):
        paper_bool_df.loc[i,:] = 0
        for orcid in paper:
            paper_bool_df.loc[i,orcid] = 1

    # Create papers matrix in numpy
    papers_mat = paper_bool_df.to_numpy()

    # Build collaboration vector to store results
    n_authors = len(authors_index)
    collabs_length = int(n_authors*(n_authors+1)/2 - n_authors) 
    collabs = np.zeros(shape=(collabs_length))

    # Store copy of papers_mat for iterative updating
    papers_mat_i = papers_mat

    # Initialize writing position
    start_pos = 0

    for i in range(0, n_authors-1): #last author loop is unnecessary
        print(f"Progress: {i/(n_authors-2)*100:.0f}%.", end="\r")

        # Initialize matrix
        C = np.identity(n_authors-i)
        C = C[:,1:]
        C[0] = 1

        # Main inner product
        result = np.dot(papers_mat_i, C)

        # Calculate number of collaborations
        result = result - 1
        result = result.clip(0)
        collabs_author = result.sum(axis=0)

        # Store in collabs vector
        end_pos = start_pos + n_authors - i - 1

        collabs[start_pos:end_pos] = collabs_author

        # Update start_pos for writing next loop
        start_pos = end_pos 

        # Remove first author from papers_mat for next loop
        papers_mat_i = papers_mat_i[:,1:]

    author_combinations = combinations(authors_index,2)
    collabs_df = pd.DataFrame(list(author_combinations), columns=['source', 'target'])
    collabs_df['value']=collabs

    if save:
        collabs_df.to_csv(f'./data/edgelist_{institution}.csv')

# Execution
# Load edgelist
papers_df = pd.read_csv('./data/papers.csv')
res_df = pd.read_csv('./data/nodelist.csv')

# Select institution
institution_list = ['IGTP']

for institution in institution_list:
    print(f"Institution: {institution}")
    calculate_collaborations(institution, papers_df, res_df, save=True, threshold=100)