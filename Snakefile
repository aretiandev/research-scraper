# Snakefile
#
# There are 8 target items: 
#   author and author_links
#   paper and paper_links
#   project and project_links
#   group and group_links
#
# Run everthing:
#   snakemake --cores all --log-handler-script scripts/log_handler.py
#
# Test if website is online:
#   snakemake ping_and_run -c1
#
# Useful rules:
#   snakemake --cores all --log-handler-script scripts/log_handler.py
#   snakemake ping_and_run -c1
#   snakemake dag -c1
#   snakemake rulegraph -c1

from src.process import get_date

date_today = get_date()
# date_today = '20220314'

institution_list = ['IGTP+', 'UPC_CIMNE', 'UB', 'UPF', 'UVic-UCC', 'UOC']
threads_max = 16

rule all:
    input:
        expand(f'data/{date_today}_edges_{{institution}}.csv', institution=institution_list),
        f'data/{date_today}_project_data.csv',
        f'data/{date_today}_group_data.csv'

rule ping_and_run:
    script:
        "bot/bot.py"

rule links:
    output: 
        f"data/{date_today}_{{item_name}}_links.csv"
    threads: 
        threads_max
    script: 
        "src/scrape.py"

rule data:
    input:
        f'data/{date_today}_{{item_name}}_links.csv'
        # f"data/{date_today}_group_links.csv"
    output:
        f'data/{date_today}_{{item_name}}_data.csv'
    threads: threads_max
    params:
        batch_size=200
    script:
        "src/scrape.py"

rule clean:
    input:
        f'data/{date_today}_{{item_name}}_data.csv'
    output:
        f'data/{date_today}_{{item_name}}_clean.csv'
    script:
        "scripts/clean.py"

rule filter_authors:
    input:
        f'data/{date_today}_author_clean.csv'
    output:
        f'data/{date_today}_nodes_{{institution}}.csv'
    script:
        "scripts/filter_authors.py"

rule filter_papers:
    input:
        f'data/{date_today}_nodes_{{institution}}.csv',
        f'data/{date_today}_paper_clean.csv'
    output:
        f'data/{date_today}_papers_{{institution}}.csv',
        f'data/{date_today}_papers_{{institution}}_2plus.csv'
    script:
        "scripts/filter_papers.py"

rule create_edges:
    input:
        f'data/{date_today}_nodes_{{institution}}.csv',
        f'data/{date_today}_papers_{{institution}}_2plus.csv'
    output:
        f'data/{date_today}_edges_{{institution}}.csv'
    script:
        "scripts/create_edges.py"

rule dag:
    output:
        f"figs/{date_today}_dag.png"
    shell:
        "snakemake --dag | dot -Tpng > {output}"

rule rulegraph:
    output:
        f"figs/{date_today}_rulegraph.png"
    shell:
        "snakemake --rulegraph | dot -Tpng > {output}"
