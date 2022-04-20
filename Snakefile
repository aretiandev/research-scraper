# Snakefile
#
# There are 8 target items: 
#   author and author_urls
#   paper and paper_urls
#   project and project_urls
#   group and group_urls
#
# Run everthing:
#   snakemake --cores all --log-handler-script scripts/log_handler.py
#
# Note: this should be run from within a Docker container as root to avoid file permission errors
#
# Test if website is online:
#   snakemake ping_and_run -c1
#
# Useful rules:
#   snakemake --cores all --log-handler-script scripts/log_handler.py
#   snakemake ping_and_run -c1
#   snakemake dag -c1
#   snakemake rulegraph -c1
#
# View DAG or Rulegraph:
#   snakemake dag -c1
#   snakemake rulegraph -c1

from scripts.src.process import get_date

# Setup
# date_today = get_date()
date_today = '20220419'
institution_list = ['IGTP+', 'UPC_CIMNE', 'UB', 'UPF', 'UVic-UCC', 'UOC', 'Agrotecnio', 'CRAG', 'UdL', 'URV', 'UdG']
threads_max = 16

rule all:
    input:
        expand(f'data/{date_today}/{date_today}_edges_{{institution}}.csv', institution=institution_list),
        f'data/{date_today}/{date_today}_project_data.csv',
        f'data/{date_today}/{date_today}_group_data.csv'

rule ping_and_run:
    script:
        "scripts/ping_and_run.py"

rule urls:
    output: 
        f"data/{date_today}/{date_today}_{{item_name}}_urls.csv"
    threads: 
        threads_max
    script: 
        "scripts/scrape.py"

rule data:
    input:
        f'data/{date_today}/{date_today}_{{item_name}}_urls.csv'
    output:
        f'data/{date_today}/{date_today}_{{item_name}}_data.csv'
    threads: threads_max
    params:
        batch_size=200
    script:
        "scripts/scrape.py"

rule clean:
    input:
        f'data/{date_today}/{date_today}_{{item_name}}_data.csv'
    output:
        f'data/{date_today}/{date_today}_{{item_name}}_clean.csv'
    script:
        "scripts/clean.py"

rule filter_authors:
    input:
        f'data/{date_today}/{date_today}_author_clean.csv'
    output:
        f'data/{date_today}/{date_today}_nodes_{{institution}}.csv'
    script:
        "scripts/filter_authors.py"

rule filter_papers:
    input:
        f'data/{date_today}/{date_today}_nodes_{{institution}}.csv',
        f'data/{date_today}/{date_today}_paper_clean.csv'
    output:
        f'data/{date_today}/{date_today}_papers_{{institution}}.csv',
        f'data/{date_today}/{date_today}_papers_{{institution}}_2plus.csv'
    script:
        "scripts/filter_papers.py"

rule create_edges:
    input:
        f'data/{date_today}/{date_today}_nodes_{{institution}}.csv',
        f'data/{date_today}/{date_today}_papers_{{institution}}_2plus.csv'
    output:
        f'data/{date_today}/{date_today}_edges_{{institution}}.csv'
    script:
        "scripts/create_edges.py"

rule dag:
    output:
        f"figs/{date_today}/{date_today}_dag.png"
    shell:
        "snakemake --dag | dot -Tpng > {output}"

rule rulegraph:
    output:
        f"figs/{date_today}/{date_today}_rulegraph.png"
    shell:
        "snakemake --rulegraph | dot -Tpng > {output}"
