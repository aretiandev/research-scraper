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

    author_df = pd.read_csv(input_nodes, converters={"groups": eval})

    # Filter nodelist for researchers with nonempty research groups
    mask = author_df["groups"].apply(len) > 0
    author_gp_df = author_df[mask]

    # Clean data
    author_gp_lst = list(author_gp_df["id"].unique())

    author_gp_df = author_gp_df.copy()
    author_gp_df.loc[:,"url_id"] = author_gp_df["groups"].apply(lambda x: x[0][1:])

    # Get group names
    group_df = pd.read_csv(input_groups)
    group_df["url_id"] = group_df["url"].str[31:]

    if 'researcher_names' in group_df.columns:
        group_df = group_df[["name", "principal_names", "researcher_names", "url_id"]]    
        group_df['principal_names'] = group_df['principal_names'].apply(lambda x: x if isinstance(x, list) else ([] if pd.isna(x) else [x]))
        group_df['researcher_names'] = group_df['researcher_names'].apply(lambda x: x if isinstance(x, list) else ([] if pd.isna(x) else [x]))

        # Combine the lists from 'principal_names' and 'researcher_names' into 'total_names'
        group_df['total_names'] = group_df.apply(lambda row: row['principal_names'] + row['researcher_names'], axis=1)
    else:
        group_df = group_df[["name", "principal_names", "url_id"]]    
        group_df =group_df.rename(columns={'principal_names':'total_names'})

    author_gp_df = author_gp_df.merge(group_df, how="left", on="url_id")
    author_gp_df = author_gp_df[['id', 'label', 'institution_2', 'department', 'institution', 'institution_group',
            'n_publications', 'n_articles', 'n_chapters', 'n_books', 'n_other',
            'url_id', 'name', 'total_names']]

    qual_nodes_df = author_gp_df.groupby("url_id")[["name","institution","institution_2", "institution_group", "department","total_names"]].first().reset_index()
    quant_nodes_df = author_gp_df.groupby("url_id").sum()

    nodes_df = qual_nodes_df.merge(quant_nodes_df, how='left', on='url_id')
    nodes_df = nodes_df.rename(columns={"url_id": "id", "name": "label"})

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
        author_gp_df[["id", "url_id"]], how="left", left_on="Source", right_on="id"
    )
    edges_gp_df = edges_gp_df.rename(columns={"url_id": "Source_gp"})
    edges_gp_df = edges_gp_df.merge(
        author_gp_df[["id", "url_id"]], how="left", left_on="Target", right_on="id"
    )
    edges_gp_df = edges_gp_df.rename(columns={"url_id": "Target_gp"})
    edges_gp_df = edges_gp_df[["Source_gp", "Target_gp", "Weight"]]
    edges_gp_df.columns = ["Source", "Target", "Weight"]
    edges_gp_df = edges_gp_df.groupby(["Source", "Target"]).sum().reset_index()

    # Save
    edges_gp_df.to_csv(output_group_edges, index=None)
    log.info(f"Saved: {output_group_edges}")
