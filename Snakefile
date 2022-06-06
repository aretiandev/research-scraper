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

from datetime import datetime
from pathlib import Path
from scripts.src.utils import get_date
from dotenv import load_dotenv
load_dotenv()


# Configuration
configfile: 'config.yml'

date_today = os.environ.get('date') or config.get('date') or get_date()
config['date'] = date_today
date_now = datetime.now().strftime("%Y%m%d.%H%M%S")
Path(config["logfile"]).parent.mkdir(exist_ok=True)
logfile = Path(config["logfile"]).parent / f"{date_now}_{Path(config['logfile']).name}"
config["logfile"] = logfile
institution_list = config['institutions']
threads_max = config['threads_max']
timeout = config['timeout']

slack_bot_token = os.environ.get('SLACK_BOT_TOKEN')
slack_member_id = os.environ.get('SLACK_MEMBER_ID')
signing_secret = os.environ.get("SIGNING_SECRET")

n_catalog_urls = config['n_catalog_urls']
n_papers = config.get('n_papers') or 10000

Path('output').mkdir(exist_ok=True)


rule all: 
  input: expand("output/papers/{paper_number}.html", paper_number=range(n_papers))

rule fetch_catalog:
  output: expand("output/catalog/{url_number}.html", url_number=range(n_catalog_urls))
  params: target = "catalog"
  script: "scripts/fetch_html.py"

rule parse_catalog:
  input: "output/catalog/{url_number}.html"
  output: touch("output/catalog/{url_number}.done")
  script: "scripts/parse_html.py"

rule fetch_papers:
  # input: expand("output/catalog/{url_number}.done", url_number=range(n_catalog_urls))
  output: expand("output/papers/{paper_number}.html", paper_number=range(n_papers))
  script: "scripts/fetch_html.py"

rule prune:
  script: "scripts/prune.py"

########

rule ping_and_run:
    params:
        slack_bot_token = slack_bot_token,
        slack_member_id = slack_member_id
    script:
        "scripts/ping_and_run.py"

  
# rule urls:
#     output: 
#         f"data/{date_today}/{date_today}_{{item_name}}_urls.csv"
#     threads: 
#         threads_max
#     params:
#         batch_size = 50,
#         timeout = timeout
#     script: 
#         "scripts/scrape.py"

rule data:
    input:
        f'data/{date_today}/{date_today}_{{item_name}}_urls.csv'
    output:
        f'data/{date_today}/{date_today}_{{item_name}}_data.csv'
    threads: threads_max
    params:
        batch_size = 50,
        timeout = timeout
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
        f'data/{date_today}/{date_today}_author_clean.csv',
    output:
        f'data/{date_today}/{date_today}_nodestemp_{{institution}}.csv'
    script:
        "scripts/filter_authors.py"

rule filter_papers:
    input:
        f'data/{date_today}/{date_today}_nodestemp_{{institution}}.csv',
        f'data/{date_today}/{date_today}_paper_clean.csv'
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
        f'data/{date_today}/{date_today}_group_data.csv',
    output:
        f'data/{date_today}/{date_today}_nodes_{{institution}}.csv'
    script:
        "scripts/add_nodes_stats.py"

rule create_group_networks:
    input:
        f'data/{date_today}/{date_today}_nodes_{{institution}}.csv',
        f'data/{date_today}/{date_today}_group_data.csv',
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
