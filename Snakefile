# Snakefile
#
# Scrapes information from Portal de la Recerca for authors, papers, projects and groups
# and creates nodes and edges for network visualization.
#
# Set date (if not set will default to today):
#   export DATE=20220601
#
# Run everthing:
#   snakemake --cores all
#
# Run with Slack notifications
#   snakemake --cores all --log-handler-script scripts/slack.py
#
# View DAG or Rulegraph:
#   snakemake dag -c1
#   snakemake rulegraph -c1
#
# Test if website is online:
#   snakemake ping_and_run -c1
#
# Tips:
#   To avoid crashing the server set a low value for params.batch_size

from config import Config

date_today = Config.DATE
institution_list = Config.INSTITUTION_LIST
threads_max = Config.THREADS_MAX
timeout = Config.TIMEOUT
SLACK_BOT_TOKEN = Config.SLACK_BOT_TOKEN
SLACK_MEMBER_ID = Config.SLACK_MEMBER_ID
database = Config.DATABASE
batch_size = Config.BATCH_SIZE

rule all:
    input:
        expand(f'data/{date_today}/{date_today}_nodes_{{institution}}.csv', institution=institution_list),
        expand(f'data/{date_today}/{date_today}_edges_{{institution}}.csv', institution=institution_list),
        expand(f'data/{date_today}/{date_today}_group_nodes_{{institution}}.csv', institution=institution_list),
        expand(f'data/{date_today}/{date_today}_group_edges_{{institution}}.csv', institution=institution_list),
        expand(f'data/{date_today}/{date_today}_project_data_{{institution}}.csv', institution=institution_list)

rule ping_and_run:
    params:
        SLACK_BOT_TOKEN = SLACK_BOT_TOKEN,
        SLACK_MEMBER_ID = SLACK_MEMBER_ID
    script:
        "scripts/ping_and_run.py"

rule urls:
    output: 
        f"data/{date_today}/{date_today}_{{item_name}}_urls_{{institution}}.csv"
    threads: 
        threads_max
    params:
        batch_size = batch_size,
        timeout = timeout,
        database = database,
        institution_list = institution_list
    script: 
        "scripts/scrape.py"

rule data:
    input:
        f'data/{date_today}/{date_today}_{{item_name}}_urls_{{institution}}.csv'
    output:
        f'data/{date_today}/{date_today}_{{item_name}}_data_{{institution}}.csv'
    threads: threads_max
    params:
        batch_size = batch_size,
        timeout = timeout,
        database = database
    script:
        "scripts/scrape.py"

rule clean:
    input:
        f'data/{date_today}/{date_today}_{{item_name}}_data_{{institution}}.csv'
    output:
        f'data/{date_today}/{date_today}_{{item_name}}_clean_{{institution}}.csv'
    script:
        "scripts/clean.py"

rule filter_authors:
    input:
        f'data/{date_today}/{date_today}_author_clean_{{institution}}.csv',
    output:
        f'data/{date_today}/{date_today}_nodestemp_{{institution}}.csv'
    script:
        "scripts/filter_authors.py"

rule filter_papers:
    input:
        f'data/{date_today}/{date_today}_nodestemp_{{institution}}.csv',
        f'data/{date_today}/{date_today}_paper_clean_{{institution}}.csv'
    output:
        f'data/{date_today}/{date_today}_papers_{{institution}}.csv',
        f'data/{date_today}/{date_today}_papers_{{institution}}_2plus.csv'
    script:
        "scripts/filter_papers.py"

rule create_edges:
    input:
        f'data/{date_today}/{date_today}_nodestemp_{{institution}}.csv',
        f'data/{date_today}/{date_today}_papers_{{institution}}_2plus.csv'
    output:
        f'data/{date_today}/{date_today}_edges_{{institution}}.csv'
    script:
        "scripts/create_edges.py"

rule add_nodes_stats:
    input:
        f'data/{date_today}/{date_today}_nodestemp_{{institution}}.csv',
        f'data/{date_today}/{date_today}_papers_{{institution}}.csv',
        f'data/{date_today}/{date_today}_group_data_{{institution}}.csv',
    output:
        f'data/{date_today}/{date_today}_nodes_{{institution}}.csv'
    script:
        "scripts/add_nodes_stats.py"

rule create_group_networks:
    input:
        f'data/{date_today}/{date_today}_nodes_{{institution}}.csv',
        f'data/{date_today}/{date_today}_group_data_{{institution}}.csv',
        f'data/{date_today}/{date_today}_edges_{{institution}}.csv'
    output:
        f'data/{date_today}/{date_today}_group_nodes_{{institution}}.csv',
        f'data/{date_today}/{date_today}_group_edges_{{institution}}.csv'
    script:
        "scripts/create_group_networks.py"

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

rule create_symlinks:
    params:
        target="20220419",
        date_today=date_today
    shell:
        """
        mkdir -p data/{date_today}
        ln -s ../{params.target}/{params.target}_author_urls.csv data/{date_today}/{date_today}_author_urls.csv
        ln -s ../{params.target}/{params.target}_author_data.csv data/{date_today}/{date_today}_author_data.csv
        ln -s ../{params.target}/{params.target}_paper_urls.csv data/{date_today}/{date_today}_paper_urls.csv
        ln -s ../{params.target}/{params.target}_paper_data.csv data/{date_today}/{date_today}_paper_data.csv
        ln -s ../{params.target}/{params.target}_group_urls.csv data/{date_today}/{date_today}_group_urls.csv
        ln -s ../{params.target}/{params.target}_group_data.csv data/{date_today}/{date_today}_group_data.csv
        ln -s ../{params.target}/{params.target}_project_urls.csv data/{date_today}/{date_today}_project_urls.csv
        ln -s ../{params.target}/{params.target}_project_data.csv data/{date_today}/{date_today}_project_data.csv
        """
