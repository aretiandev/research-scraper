#!/usr/bin/env python
# coding: utf-8
"""
Groups module: collection of helper functions to create nodes and edges for groups.
"""

import pandas as pd
import logging
import sys

log = logging.getLogger(__name__)


def create_group_networks(
    input_nodes, input_groups, input_edges, output_group_nodes, output_group_edges
):

    # Get Nodes
    log.info(f"Loading: {input_nodes}")
    author_df = pd.read_csv(input_nodes, converters={"groups": eval})
    
    def safe_eval(x):
        if isinstance(x, str):  # Only evaluate if x is a string
            return eval(x)
        return []
    
    author_df = author_df.loc[:, ~author_df.columns.str.contains('^Unnamed')]
    author_df = author_df.loc[:, ~author_df.columns.str.contains('Weight')]
    author_df = author_df.loc[:, ~author_df.columns.str.contains('Target')]
    
    #calculate researcher node degree
    edges_df = pd.read_csv(input_edges)
    degree_df = edges_df.groupby('Source').sum()
    author_df = author_df.merge(degree_df, how='left', left_on='id', right_on='Source')
    author_df = author_df.drop(columns=['Target'])
    #author_df['Weight'] = author_df['Weight'].fillna(0)

    author_df.to_csv(input_nodes, index=False)
    # Filter nodelist for researchers with nonempty research groups
    mask = author_df["groups"].apply(len) > 0
    author_gp_df = author_df[mask]

    # Clean data
    author_gp_lst = list(author_gp_df["id"].unique())

    author_gp_df = author_gp_df.copy()
    
    full_author_gp_df = author_gp_df.explode('groups').reset_index(drop=True)

    full_author_gp_df["url_id"] = full_author_gp_df["groups"].str.slice(1)

    # Get group names
    
    try:
        group_df = pd.read_csv(input_groups)
        group_df["url_id"] = group_df["url"].str[31:]
    except (pd.errors.EmptyDataError, FileNotFoundError) as e:
        log.warning(f"No data found in groups file {input_groups}. Skipping group processing.")
        nodes_df = pd.DataFrame()
        edges_gp_df = pd.DataFrame()
        nodes_df.to_csv(output_group_nodes, index=None)
        edges_gp_df.to_csv(output_group_edges, index=None)
        return

 
    if 'researcher_names' in group_df.columns:
        group_df = group_df[["name", "principal_names", "researcher_names", "url_id"]]    
        group_df['principal_names'] = group_df['principal_names'].apply(safe_eval)
        group_df['researcher_names'] = group_df['researcher_names'].apply(safe_eval)

        # Combine the lists from 'principal_names' and 'researcher_names' into 'total_names'
        group_df['total_names'] = group_df.apply(lambda row: row['principal_names'] + row['researcher_names'], axis=1)

    else:
        group_df = group_df[["name", "principal_names", "url_id"]]    
        group_df = group_df.rename(columns={'principal_names':'total_names'})
        group_df['total_names'] = group_df['total_names'].apply(safe_eval)
    
    group_df['n_researchers'] = group_df['total_names'].apply(len)

    full_author_gp_df = full_author_gp_df.merge(group_df, how="left", on="url_id")
    full_author_gp_df = full_author_gp_df[['id', 'url_id', 'name', 'institution', 'institution_2','institution_group', 'department', 'total_names', 'n_researchers', 'projects', 'n_projects',
            'n_publications', 'n_articles', 'n_chapters', 'n_books', 'n_other', 'Weight']]
    full_author_gp_df['projects'] = full_author_gp_df['projects'].apply(safe_eval)
 
    columns_to_sum = ['projects', 'n_projects', 'n_publications', 'n_articles', 'n_chapters', 'n_books', 'n_other']

    # Create the initial aggregation dictionary
    agg_dict = {col: 'first' for col in full_author_gp_df.columns if col not in columns_to_sum}
    agg_dict.update({col: 'sum' for col in columns_to_sum})
    agg_dict['Weight'] = 'mean'

    # Update the rule for 'projects' to concatenate lists instead of summing
    agg_dict['projects'] = lambda x: sum(x, [])
    
    nodes_df = full_author_gp_df.groupby('url_id').agg(agg_dict)
    nodes_df = nodes_df.drop(columns=['id'])
    nodes_df = nodes_df.rename(columns={"url_id": "id", "name": "label"})
    nodes_df = nodes_df.dropna(subset=['label'])
    nodes_df = nodes_df.rename(columns={"Weight": "avg_degree"})

    
    #calculate group node degree
#     group_edges_df = pd.read_csv(group_edges)
#     group_degree_df = group_edges_df.groupby('Source')

#     adj_gp_degree_df = group_degree_df['Weight'].apply(lambda x: x.iloc[1:].sum())

    #adj_gp_degree_df = adj_gp_degree_df.reset_index(name='node_degree')
    # nodes_df = nodes_df.merge(adj_gp_degree_df, how='left', left_on='id', right_on='Source')
    #nodes_df = nodes_df.drop(columns=['Source'])
    #nodes_df['Weight'] = nodes_df['Weight'].fillna(0)
    # Save
    nodes_df.to_csv(output_group_nodes, index=None)
    log.info(f"Saved: {output_group_nodes}")

    # Create group level edgelist
    log.info(f"Loading: {input_edges}")
    edges_df = pd.read_csv(input_edges)

    mask = edges_df.apply(
        lambda row: row["Source"] in author_gp_lst and row["Target"] in author_gp_lst,
        axis=1,
    )
    edges_gp_df = edges_df[mask]
    edges_gp_df = edges_gp_df.merge(
        full_author_gp_df[["id", "url_id"]], how="left", left_on="Source", right_on="id"
    )
    #you've effectively changed author_gp_df to full_author_gp_df line 50
    edges_gp_df = edges_gp_df.rename(columns={"url_id": "Source_gp"})
    edges_gp_df = edges_gp_df.merge(
        full_author_gp_df[["id", "url_id"]], how="left", left_on="Target", right_on="id"
    )
    edges_gp_df = edges_gp_df.rename(columns={"url_id": "Target_gp"})
    edges_gp_df = edges_gp_df[["Source_gp", "Target_gp", "Weight"]]
    edges_gp_df.columns = ["Source", "Target", "Weight"]
    edges_gp_df = edges_gp_df.groupby(["Source", "Target"]).sum().reset_index()

    # Save
    edges_gp_df.to_csv(output_group_edges, index=None)
    log.info(f"Saved: {output_group_edges}")
    
    
def filter_group_edges(
    input_group_data, input_group_edges, output_internal_gp_edges):
    
    log.info(f"Loading: {input_group_data}")
    try:
        group_df = pd.read_csv(input_group_data)
        group_df["url_id"] = group_df["url"].str[31:]
    except (pd.errors.EmptyDataError, FileNotFoundError) as e:
        log.warning(f"No data found in groups file {input_group_data}. Skipping group edges filtering.")
        filtered_df = pd.DataFrame()
        filtered_df.to_csv(output_internal_gp_edges, index=None)
        return
    
    log.info(f"Loading: {input_group_edges}")
    gp_edges_df = pd.read_csv(input_group_edges)

    

    gp_list = group_df['url_id'].tolist()

    filtered_df = gp_edges_df[(gp_edges_df['Target'].isin(gp_list)) & (gp_edges_df['Source'].isin(gp_list))]

    filtered_df.to_csv(output_internal_gp_edges, index=None)
    log.info(f"Saved: {output_internal_gp_edges}")
