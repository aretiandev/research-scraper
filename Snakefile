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

# Load modules
import pandas as pd
from datetime import date
import asyncio

from src.scrape import scrape
from src.process import filter_authors_and_papers, add_variables_to_nodelist
from src.edges import create_edgelist

# Setup global variables
date_today = date.today().strftime("%Y%m%d")
url_root = 'https://portalrecerca.csuc.cat'
institution_list = ['IGTP+', 'UPC_CIMNE', 'UB', 'UPF', 'UVic-UCC', 'UOC']

rule all:
    input:
        f'data/paper_data_{date_today}.csv',
        f'data/author_data_{date_today}.csv',
        f'data/project_data_{date_today}.csv',
        f'data/group_data_{date_today}.csv',
        expand(f'./data/edges_{{institution}}_{date_today}.csv', institution=institution_list)

rule item_links:
    output:
        f"data/{{item_name}}_links_{date_today}.csv"
    run:
        batch_size = 10
        # item_name_links = wildcards.item_name + '_links'
        asyncio.run(scrape(items=wildcards.item_name+'_links', batch_size=batch_size, out_file=output[0]))

rule clean_links:
    input:
        f"data/{{item_name}}_links_{date_today}.csv"
    shell:
        "rm {input}"

rule item_data:
    input:
        f'data/{{item_name}}_links_{date_today}.csv'
    output:
        f'data/{{item_name}}_data_{date_today}.csv'
    run:
        if wilcards.item_name == 'author':
            batch_size = 100
        else:
            batch_size = 10
        item_urls = pd.read_csv(input[0])
        item_urls = list(item_urls['0'])
        if wildcards.item_name == 'paper':
            urls = [url_root + url + '?mode=full' for url in item_urls]
        else:
            urls = [url_root + url for url in item_urls]
        asyncio.run(scrape(items=wildcards.item_name, urls=urls, batch_size=batch_size, out_file=output[0]))

rule filter_authors_and_papers:
    input:
        "data/papers.csv",
        f'./data/nodes_{date_today}.csv'
    output:
        f'data/papers_{{institution}}_{date_today}.csv',
        f'data/papers_{{institution}}_2plus_{date_today}.csv',
        f'data/nodes_{{institution}}_{date_today}.csv'
    run:
        filter_authors_and_papers()

rule add_variables_to_nodelist:
    input: 
        f'./data/nodes_{{institution}}_{date_today}.csv',
        f'./data/papers_{{institution}}_{date_today}.csv'
    output:
        f'./data/nodes_{{institution}}_{date_today}.csv'
    run:
        add_variables_to_nodelist()

rule create_edges:
    input:
        f'data/papers_{{institution}}_2plus_{date_today}.csv',
        f'data/nodes_{{institution}}_{date_today}.csv'
    output:
        f'./data/edges_{{institution}}_{date_today}.csv'
    run:
        create_edgelist(institution)