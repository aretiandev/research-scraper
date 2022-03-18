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
# Run everthing:
#   snakemake --cores all
#
# Test if website is online:
#   snakemake ping_and_notify -c1
#
# Full list of rules/targets:
#   snakemake -c1
#   snakemake ping_and_notify -c1
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
#   snakemake dag -c1
#   snakemake rulegraph -c1

import pandas as pd
from datetime import date
import asyncio
import pytz
from dotenv import load_dotenv

from src.scrape import scrape
from src.process import clean_authors, clean_papers, filter_authors, filter_papers, get_date
from src.edges import create_edgelist
from bot.bot import ping_and_wait, send_slack_message

# Setup Slack Bot
load_dotenv()
slack_token = os.environ["SLACK_BOT_TOKEN"]
slack_member_id = os.environ["SLACK_MEMBER_ID"]
channel = slack_member_id
error_msg = """:snake::warning: There was an error in rule '{rule_name}':
{error}"""
running_msg = ":snake::runner: Running rule {rule}."
success_msg = ":snake::ok_hand: Rule {rule} completed."

# Get date
date_today = get_date()
# date_today = '20220314'

url_root = 'https://portalrecerca.csuc.cat'
institution_list = ['IGTP+', 'UPC_CIMNE', 'UB', 'UPF', 'UVic-UCC', 'UOC']

rule all:
    input:
        expand(f'data/{date_today}_edges_{{institution}}.csv', institution=institution_list),
        f'data/{date_today}_project_data.csv',
        f'data/{date_today}_group_data.csv'

rule ping_stackoverflow:
    run:
        url = 'https://stackoverflow.com/'
        ping_and_wait(url)
        message = ":globe_with_meridians: Stackoverflow is back online."
        send_slack_message(channel, message, slack_token)

rule ping_and_notify:
    script:
        "bot/bot.py"

rule author_links:
    output:
        f"data/{date_today}_author_links.csv"
    run:
        try:
            batch_size = 20
            send_slack_message(channel,running_msg.format(rule=rule),slack_token=slack_token)
            asyncio.run(scrape(items='author_links', batch_size=batch_size, out_file=output[0]))
            send_slack_message(channel,success_msg.format(rule=rule),slack_token=slack_token)
        except Exception as e:
            send_slack_message(channel,error_msg.format(rule_name=rule,error=e),slack_token=slack_token)
            print(e)
            raise e

rule paper_links:
    input:
        f"data/{date_today}_author_links.csv"
    output:
        f"data/{date_today}_paper_links.csv"
    run:
        try:
            batch_size = 50
            send_slack_message(channel,running_msg.format(rule=rule),slack_token=slack_token)
            asyncio.run(scrape(items='paper_links', batch_size=batch_size, out_file=output[0]))
            send_slack_message(channel,success_msg.format(rule=rule),slack_token=slack_token)
        except Exception as e:
            send_slack_message(channel,error_msg.format(rule_name=rule,error=e),slack_token=slack_token)
            print(e)
            raise e

rule project_links:
    input:
        f"data/{date_today}_paper_links.csv"
    output:
        f"data/{date_today}_project_links.csv"
    run:
        try:
            batch_size = 20
            send_slack_message(channel,running_msg.format(rule=rule),slack_token=slack_token)
            asyncio.run(scrape(items='project_links', batch_size=batch_size, out_file=output[0]))
            send_slack_message(channel,success_msg.format(rule=rule),slack_token=slack_token)
        except Exception as e:
            send_slack_message(channel,error_msg.format(rule_name=rule,error=e),slack_token=slack_token)
            print(e)
            raise e

rule group_links:
    input:
        f"data/{date_today}_project_links.csv"
    output:
        f"data/{date_today}_group_links.csv"
    run:
        try:
            batch_size = 20
            send_slack_message(channel,running_msg.format(rule=rule),slack_token=slack_token)
            asyncio.run(scrape(items='group_links', batch_size=batch_size, out_file=output[0]))
            send_slack_message(channel,success_msg.format(rule=rule),slack_token=slack_token)
        except Exception as e:
            send_slack_message(channel,error_msg.format(rule_name=rule,error=e),slack_token=slack_token)
            print(e)
            raise e
        
rule author_data:
    input:
        f'data/{date_today}_author_links.csv',
        f"data/{date_today}_group_links.csv"
    output:
        f'data/{date_today}_author_data.csv'
    run:
        try:
            batch_size = 200
            item_urls = pd.read_csv(input[0])
            item_urls = list(item_urls['0'])
            urls = [url_root + url for url in item_urls]
            send_slack_message(channel,running_msg.format(rule=rule),slack_token=slack_token)
            asyncio.run(scrape(items='author', urls=urls, batch_size=batch_size, out_file=output[0]))
            send_slack_message(channel,success_msg.format(rule=rule),slack_token=slack_token)
        except Exception as e:
            send_slack_message(channel,error_msg.format(rule_name=rule,error=e),slack_token=slack_token)
            print(e)
            raise e
        
rule paper_data:
    input:
        f'data/{date_today}_paper_links.csv',
        f"data/{date_today}_author_data.csv"
    output:
        f'data/{date_today}_paper_data.csv'
    run:
        try:
            batch_size = 200
            item_urls = pd.read_csv(input[0])
            item_urls = list(item_urls['0'])
            urls = [url_root + url + '?mode=full' for url in item_urls]
            send_slack_message(channel,running_msg.format(rule=rule),slack_token=slack_token)
            asyncio.run(scrape(items='paper', urls=urls, batch_size=batch_size, out_file=output[0]))
            send_slack_message(channel,success_msg.format(rule=rule),slack_token=slack_token)
        except Exception as e:
            send_slack_message(channel,error_msg.format(rule_name=rule,error=e),slack_token=slack_token)
            print(e)
            raise e

rule project_data:
    input:
        f'data/{date_today}_project_links.csv',
        f"data/{date_today}_paper_data.csv"
    output:
        f'data/{date_today}_project_data.csv'
    run:
        try:
            batch_size = 200
            item_urls = pd.read_csv(input[0])
            item_urls = list(item_urls['0'])
            urls = [url_root + url for url in item_urls]
            send_slack_message(channel,running_msg.format(rule=rule),slack_token=slack_token)
            asyncio.run(scrape(items='project', urls=urls, batch_size=batch_size, out_file=output[0]))
            send_slack_message(channel,success_msg.format(rule=rule),slack_token=slack_token)
        except Exception as e:
            send_slack_message(channel,error_msg.format(rule_name=rule,error=e),slack_token=slack_token)
            print(e)
            raise e

rule group_data:
    input:
        f'data/{date_today}_group_links.csv',
        f"data/{date_today}_project_data.csv"
    output:
        f'data/{date_today}_group_data.csv'
    run:
        try:
            batch_size = 50
            item_urls = pd.read_csv(input[0])
            item_urls = list(item_urls['0'])
            urls = [url_root + url for url in item_urls]
            send_slack_message(channel,running_msg.format(rule=rule),slack_token=slack_token)
            asyncio.run(scrape(items='group', urls=urls, batch_size=batch_size, out_file=output[0]))
            send_slack_message(channel,success_msg.format(rule=rule),slack_token=slack_token)
        except Exception as e:
            send_slack_message(channel,error_msg.format(rule_name=rule,error=e),slack_token=slack_token)
            print(e)
            raise e

rule clean_authors:
    input:
        f'data/{date_today}_author_data.csv'
    output:
        f'data/{date_today}_author_clean.csv'
    run:
        try:
            send_slack_message(channel,running_msg.format(rule=rule),slack_token=slack_token)
            clean_authors(date_today)
            send_slack_message(channel,success_msg.format(rule=rule),slack_token=slack_token)
        except Exception as e:
            send_slack_message(channel,error_msg.format(rule_name=rule,error=e),slack_token=slack_token)
            print(e)
            raise e

rule clean_papers:
    input:
        f'data/{date_today}_paper_data.csv'
    output:
        f'data/{date_today}_paper_clean.csv'
    run:
        try:
            send_slack_message(channel,running_msg.format(rule=rule),slack_token=slack_token)
            clean_papers(date_today)
            send_slack_message(channel,success_msg.format(rule=rule),slack_token=slack_token)
        except Exception as e:
            send_slack_message(channel,error_msg.format(rule_name=rule,error=e),slack_token=slack_token)
            print(e)
            raise e

rule filter_authors:
    input:
        f'data/{date_today}_author_clean.csv'
    output:
        f'data/{date_today}_nodes_{{institution}}.csv'
    run:
        try:
            send_slack_message(channel,running_msg.format(rule=rule),slack_token=slack_token)
            filter_authors(wildcards.institution, date_today)
            send_slack_message(channel,success_msg.format(rule=rule),slack_token=slack_token)
        except Exception as e:
            rule_name = f"{rule}:({wildards.institution})"
            send_slack_message(channel,error_msg.format(rule_name=rule_name,error=e),slack_token=slack_token)
            print(e)
            raise e

rule filter_papers:
    input:
        f'data/{date_today}_nodes_{{institution}}.csv',
        f'data/{date_today}_paper_clean.csv'
    output:
        f'data/{date_today}_papers_{{institution}}.csv',
        f'data/{date_today}_papers_{{institution}}_2plus.csv'
    run:
        try:
            send_slack_message(channel,running_msg.format(rule=rule),slack_token=slack_token)
            filter_papers(wildcards.institution, date_today)
            send_slack_message(channel,success_msg.format(rule=rule),slack_token=slack_token)
        except Exception as e:
            rule_name = f"{rule}:({wildcards.institution})"
            send_slack_message(channel,error_msg.format(rule_name=rule_name,error=e),slack_token=slack_token)
            print(e)
            raise e

rule create_edges:
    input:
        f'data/{date_today}_papers_{{institution}}_2plus.csv',
        f'data/{date_today}_nodes_{{institution}}.csv'
    output:
        f'data/{date_today}_edges_{{institution}}.csv'
    run:
        try:
            send_slack_message(channel,running_msg.format(rule=rule),slack_token=slack_token)
            create_edgelist(wildcards.institution, date_today)
            send_slack_message(channel,success_msg.format(rule=rule),slack_token=slack_token)
        except Exception as e:
            rule_name = f"{rule}:({wildcards.institution})"
            send_slack_message(channel,error_msg.format(rule_name=rule_name,error=e),slack_token=slack_token)
            print(e)
            raise e

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
