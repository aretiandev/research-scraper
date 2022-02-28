#!/usr/bin/env python
# coding: utf-8

# # Create edgetlist
# 
# This notebook creates the edgelist based on scraped data from Portal de la Reserca.

# # Import modules

# In[ ]:


import numpy as np
import pandas as pd
from itertools import combinations

# # Select institution

# In[ ]:


# Select institution
institution_list = ['IGTP', 'UPC', 'UB', 'UPF', 'UVic-UCC', 'UOC']

for institution in institution_list:
    print("")
    print(f"Institution: {institution}.")

    # institution = 'IGTP'
    # institution = 'UPC'
    # institution = 'UB'
    # institution = 'UPF'
    # institution = 'UVic-UCC'
    # institution = 'UOC'

    # # Create DataFrame to store collaborations

    # ## Get authors from papers

    # In[ ]:


    # Load papers with coauthors list
    print(f"{institution} - Loading papers.")
    papers = pd.read_csv(f'./data/papers_{institution}.csv', converters = {'orcids': eval}, usecols=['orcids'])
    papers = papers['orcids']

    # Get unique list of authors from papers
    authors_papers = list(set(papers.sum()))
    authors_papers.sort()

    # ## Get authors from institution

    # In[ ]:


    # Get list of authors from institution
    print(f"{institution} - Loading nodes.")
    authors_inst = pd.read_csv(f'./data/nodes_{institution}.csv', usecols=['id'])
    authors_inst = authors_inst['id']
    authors_inst = authors_inst.unique()
    authors_inst.sort()

    # ## Combine authors

    # In[ ]:


    # Combine both
    authors_index = list(set(authors_papers) & set(authors_inst))
    authors_index.sort()

    # ## Create df to store collaborations

    # In[ ]:


    print(f"{institution} - Calculting combinations of authors.")
    author_combinations = combinations(authors_index,2)
    collabs_df = pd.DataFrame(list(author_combinations), columns=['Source', 'Target'])
    collabs_df['Weight'] = 0
    collabs_df = collabs_df.set_index(['Source', 'Target'])

    # # Calculate collaborations

    # ## Main loop: add collaborations to df

    # In[ ]:


    print(f"{institution} - Main loop: counting collaborations.")
    for i, paper in enumerate(papers):
        print(f"{institution} - Progress: {i/len(papers)*100:.0f}%. ({i:,.0f}/{len(papers):,.0f}).", end="\r")
        paper = list(set(paper))
        paper.sort()
        author_pairs = combinations(paper, 2)
        for pair in author_pairs:
            try:
                collabs_df.loc[pair] += 1
            except:
                pass

    collabs_df = collabs_df.reset_index()

    # ## Save

    # In[ ]:


    collabs_df.to_csv(f'./data/edgelist_{institution}.csv', index=None)
    print(f"{institution} - Done. Saved './data/edgelist_{institution}.csv'.")

# # Check results

# In[ ]:


# Check for pairs with joint publications
# pos_pairs = collabs_df.loc[collabs_df['value'] > 0]

# # Count number of errors
# errors = 0
# count = 0
# # for idx, row in collabs_df.iterrows():
# for idx, row in pos_pairs.iterrows():
#     print(f"{count}/{len(pos_pairs)}", end="\r")
#     sum=0
#     for paper in papers:
#         if row['Source'] in paper and row['Target'] in paper:
#             sum += 1
                
#     result = row['Weight'] == sum
#     if not result:
#         errors += 1
        
#     count += 1
        
# print(f"Done. Found {errors} errors.")
