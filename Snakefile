# Snakefile
#
# Instructions:
# 
# There are 8 target items: 
#   author and author_links
#   paper and paper_links
#   project and project_links
#   group and group_links
#
# List of rules:
#   snakemake data/author_data_202203011.csv -c1
#   snakemake data/author_links_202203011.csv -c1
#   snakemake data/paper_data_202203011.csv -c1
#   snakemake data/paper_links_data_202203011.csv -c1
#   snakemake data/project_data_202203011.csv -c1
#   snakemake data/project_links_202203011.csv -c1
#   snakemake data/group_data_202203011.csv -c1
#   snakemake data/group_links_202203011.csv -c1
#   snakemake clean_links -c1
#   snakemake filter_authors_and_papers -c1
#   snakemake add_variables_to_nodelist -c1
#   snakemake create_edges -c1
# To run everthing:
#   snakemake -c1

import pandas as pd
from datetime import date
import asyncio

from src.scrape import scrape
from src.process import ( clean_authors, clean_papers, filter_authors, filter_papers,
                          add_publication_stats, get_date )
from src.edges import create_edgelist

# Setup global variables
date_today = get_date()
# date_today = '20220311'
url_root = 'https://portalrecerca.csuc.cat'
institution_list = ['IGTP+', 'UPC_CIMNE', 'UB', 'UPF', 'UVic-UCC', 'UOC']

rule all:
    input:
        expand(f'data/edges_{{institution}}_{date_today}.csv', institution=institution_list)

rule author_links:
    output:
        f"data/author_links_{date_today}.csv"
    run:
        batch_size = 20
        asyncio.run(scrape(items='author_links', batch_size=batch_size, out_file=output[0]))

rule paper_links:
    input:
        f"data/author_links_{date_today}.csv"
    output:
        f"data/paper_links_{date_today}.csv"
    run:
        batch_size = 20
        asyncio.run(scrape(items='paper_links', batch_size=batch_size, out_file=output[0]))

rule project_links:
    input:
        f"data/paper_links_{date_today}.csv"
    output:
        f"data/project_links_{date_today}.csv"
    run:
        batch_size = 20
        asyncio.run(scrape(items='project_links', batch_size=batch_size, out_file=output[0]))

rule group_links:
    input:
        f"data/project_links_{date_today}.csv"
    output:
        f"data/group_links_{date_today}.csv"
    run:
        batch_size = 20
        asyncio.run(scrape(items='group_links', batch_size=batch_size, out_file=output[0]))
        
rule author_data:
    input:
        f"data/group_links_{date_today}.csv",
        f'data/author_links_{date_today}.csv'
    output:
        f'data/author_data_{date_today}.csv'
    run:
        batch_size = 100
        item_urls = pd.read_csv(input[0])
        item_urls = list(item_urls['0'])
        urls = [url_root + url for url in item_urls]
        asyncio.run(scrape(items='author', urls=urls, batch_size=batch_size, out_file=output[0]))
        
rule paper_data:
    input:
        f"data/author_data_{date_today}.csv",
        f'data/paper_links_{date_today}.csv'
    output:
        f'data/paper_data_{date_today}.csv'
    run:
        batch_size = 100
        item_urls = pd.read_csv(input[0])
        item_urls = list(item_urls['0'])
        urls = [url_root + url + '?mode=full' for url in item_urls]
        asyncio.run(scrape(items='paper', urls=urls, batch_size=batch_size, out_file=output[0]))

rule project_data:
    input:
        f"data/paper_data_{date_today}.csv",
        f'data/project_links_{date_today}.csv'
    output:
        f'data/project_data_{date_today}.csv'
    run:
        batch_size = 20
        item_urls = pd.read_csv(input[0])
        item_urls = list(item_urls['0'])
        urls = [url_root + url for url in item_urls]
        asyncio.run(scrape(items='project', urls=urls, batch_size=batch_size, out_file=output[0]))

rule group_data:
    input:
        f"data/project_data{date_today}.csv",
        f'data/group_links_{date_today}.csv'
    output:
        f'data/group_data_{date_today}.csv'
    run:
        batch_size = 20
        item_urls = pd.read_csv(input[0])
        item_urls = list(item_urls['0'])
        urls = [url_root + url for url in item_urls]
        asyncio.run(scrape(items='group', urls=urls, batch_size=batch_size, out_file=output[0]))

rule clean_authors:
    input:
        f'data/author_data_{date_today}.csv'
    output:
        f'data/author_clean_{date_today}.csv'
    run:
        clean_authors()

rule clean_papers:
    input:
        f'data/paper_data_{date_today}.csv'
    output:
        f'data/paper_clean_{date_today}.csv'
    run:
        clean_papers()

rule filter_authors:
    input:
        f'data/author_clean_{date_today}.csv'
    output:
        f'data/nodes_{{institution}}_{date_today}.csv'
    run:
        filter_authors(wildcards.institution)

rule filter_papers:
    input:
        f'data/nodes_{{institution}}_{date_today}.csv',
        f'data/paper_clean_{date_today}.csv'
    output:
        f'data/papers_{{institution}}_{date_today}.csv',
        f'data/papers_{{institution}}_2plus_{date_today}.csv'
    run:
        filter_papers(wildcards.institution)

rule add_publication_stats:
    input: 
        f'data/nodes_{{institution}}_{date_today}.csv',
        f'data/papers_{{institution}}_{date_today}.csv'
    output:
        f'data/nodes_{{institution}}_full_{date_today}.csv'
    run:
        add_publication_stats(wildcards.institution)

rule create_edges:
    input:
        f'data/papers_{{institution}}_2plus_{date_today}.csv',
        f'data/nodes_{{institution}}_{date_today}.csv'
    output:
        f'data/edges_{{institution}}_{date_today}.csv'
    run:
        create_edgelist(wildcards.institution)