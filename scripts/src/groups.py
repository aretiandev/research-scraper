#!/usr/bin/env python
# coding: utf-8
"""
Groups module: collection of helper functions to create nodes and edges for groups.
"""

import pandas as pd
import logging

log = logging.getLogger(__name__)


def create_group_networks(
    input_nodes, input_groups, input_edges, output_nodes, output_edges
):

    # Get Nodes
    log.info(f"Loading: {input_nodes}")
    author_df = pd.read_csv(input_nodes, converters={"groups": eval})

    # Filter nodelist for researchers with nonempty research groups
    mask = author_df["groups"].apply(len) > 0
    author_gp_df = author_df[mask]

    # Clean data
    author_gp_lst = list(author_gp_df["id"].unique())
    author_gp_df.loc[:, "url_id"] = author_gp_df.loc[:, "groups"].apply(
        lambda x: x[0][1:]
    )

    # Create group level nodelist

    # Get group names
    group_df = pd.read_csv(input_groups)

    group_df["url_id"] = group_df["url"].str[31:]

    group_df = group_df[["name", "url_id"]]
    author_gp_df = author_gp_df.merge(group_df, how="left", on="url_id")

    # Collapse at group level
    nodes_df = author_gp_df.groupby("url_id").first().reset_index()

    nodes_df = nodes_df[
        [
            "url_id",
            "name",
            "institution",
            "institution_2",
            "department",
            "institution_group",
            "n_publications",
        ]
    ]
    nodes_df = nodes_df.rename(columns={"url_id": "id", "name": "label"})

    # Save
    nodes_df.to_csv(output_nodes, index=None)
    log.info(f"Saved: {output_nodes}")

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
    edges_gp_df.to_csv(output_edges, index=None)
    log.info(f"Saved: {output_edges}")
