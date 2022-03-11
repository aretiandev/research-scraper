#!/usr/bin/env python
# coding: utf-8

# # Create edges
# 
# This notebook processes the publications data and achieves the following:
# 
# - Creates the edges between authors

# # Import modules
import pandas as pd
from itertools import combinations


# # Select institution
def create_edgelist(institution, save=True):
    print(f"Institution: {institution}.")
    ## Get authors from papers

    # Set date for file versions
    date_today = '20220309'

    # Load papers with coauthors list
    print(f"{institution} - Loading papers.")
    papers_df = pd.read_csv(f'./data/papers_{institution}_2plus_{date_today}.csv', converters = {'orcids': eval})
    papers = papers_df['orcids'].copy()

    # Get unique list of authors from papers
    authors_papers = list(set(papers.sum()))
    authors_papers.sort()

    ## Get authors from institution

    # Get list of authors from institution
    print(f"{institution} - Loading nodes.")
    authors_inst_df = pd.read_csv(f'./data/nodes_{institution}_{date_today}.csv')

    authors_inst = authors_inst_df['id']
    authors_inst = authors_inst.unique()
    authors_inst.sort()

    ## Combine authors

    # Combine both
    authors_index = list(set(authors_papers) & set(authors_inst))
    authors_index.sort()
    
    ## Create df to store collaborations

    print(f"{institution} - Calculting combinations of authors.")
    author_combinations = combinations(authors_index,2)
    collabs_df = pd.DataFrame(list(author_combinations), columns=['Source', 'Target'])
    collabs_df['Weight'] = 0
    collabs_df = collabs_df.set_index(['Source', 'Target'])

    # Calculate collaborations

    ## Main loop: add collaborations to df

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
        ## Save
        outfile = f'./data/edges_{institution}_{date_today}.csv'
        collabs_df.to_csv(outfile, index=None)
        print(f"{institution} - Done. Saved '{outfile}'.")
        
    return collabs_df


# # MAIN LOOP: Institution
if __name__ == "__main__":
    # Select institution
    institution_list = ['IGTP+', 'UPC_CIMNE', 'UB', 'UPF', 'UVic-UCC', 'UOC']
    for institution in institution_list:
        collabs_df = create_edgelist(institution)
