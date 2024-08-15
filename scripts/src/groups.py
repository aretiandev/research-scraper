#!/usr/bin/env python
# coding: utf-8
"""
Groups module: collection of helper functions to create nodes and edges for groups.
"""

import pandas as pd
import logging

log = logging.getLogger(__name__)


def create_group_networks(
    input_nodes, input_groups, input_edges, output_group_nodes, output_group_edges
):

    # Get Nodes
    log.info(f"Loading: {input_nodes}")
    author_df = pd.read_csv(input_nodes, converters={"groups": eval})

    # Filter nodelist for researchers with nonempty research groups
    mask = author_df["groups"].apply(len) > 0
    author_gp_df = author_df[mask]

    # Clean data
    author_gp_lst = list(author_gp_df["id"].unique())

    author_gp_df = author_gp_df.copy()
    
    full_author_gp_df = author_gp_df.explode('groups').reset_index(drop=True)

    full_author_gp_df["url_id"] = full_author_gp_df["groups"].str.slice(1)

    # Get group names
    group_df = pd.read_csv(input_groups)
    group_df["url_id"] = group_df["url"].str[31:]
    
    def safe_eval(x):
        if isinstance(x, str):  # Only evaluate if x is a string
            return eval(x)
        return []

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
    full_author_gp_df = full_author_gp_df[['url_id', 'name', 'institution', 'institution_2','institution_group', 'department', 'total_names', 'n_researchers',
            'n_publications', 'n_articles', 'n_chapters', 'n_books', 'n_other']]

    columns_to_sum = ['n_publications', 'n_articles', 'n_chapters', 'n_books', 'n_other']
    
    
    agg_dict = {col: 'first' for col in full_author_gp_df.columns if col not in columns_to_sum}
    agg_dict.update({col: 'sum' for col in columns_to_sum})
    
    nodes_df = full_author_gp_df.groupby('url_id').agg(agg_dict)
    nodes_df = nodes_df.rename(columns={"url_id": "id", "name": "label"})
    nodes_df = nodes_df.dropna(subset=['label'])
    
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
