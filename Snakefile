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
#   snakemake --cores 16
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

import pandas as pd
from datetime import date
import asyncio
import pytz
from dotenv import load_dotenv

from src.scrape import scrape
from src.process import ( clean_authors, clean_papers, filter_authors, filter_papers,
                          add_publication_stats, get_date )
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

# date_today = '20220311'
url_root = 'https://portalrecerca.csuc.cat'
institution_list = ['IGTP+', 'UPC_CIMNE', 'UB', 'UPF', 'UVic-UCC', 'UOC']

rule all:
    input:
        expand(f'data/edges_{{institution}}_{date_today}.csv', institution=institution_list),
        f'data/project_data_{date_today}.csv',
        f'data/group_data_{date_today}.csv'

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
        f"data/author_links_{date_today}.csv"
    run:
        try:
            batch_size = 20
            asyncio.run(scrape(items='author_links', batch_size=batch_size, out_file=output[0]))
        except Exception as e:
            send_slack_message(channel,error_msg.format(rule_name=rule,error=e),slack_token=slack_token)
            print(e)
            raise e

rule paper_links:
    input:
        f"data/author_links_{date_today}.csv"
    output:
        f"data/paper_links_{date_today}.csv"
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
        f"data/paper_links_{date_today}.csv"
    output:
        f"data/project_links_{date_today}.csv"
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
        f"data/project_links_{date_today}.csv"
    output:
        f"data/group_links_{date_today}.csv"
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
        f'data/author_links_{date_today}.csv',
        f"data/group_links_{date_today}.csv"
    output:
        f'data/author_data_{date_today}.csv'
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
        f'data/paper_links_{date_today}.csv',
        f"data/author_data_{date_today}.csv"
    output:
        f'data/paper_data_{date_today}.csv'
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
        f'data/project_links_{date_today}.csv',
        f"data/paper_data_{date_today}.csv"
    output:
        f'data/project_data_{date_today}.csv'
    run:
        try:
            batch_size = 20
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
        f'data/group_links_{date_today}.csv',
        f"data/project_data_{date_today}.csv"
    output:
        f'data/group_data_{date_today}.csv'
    run:
        try:
            batch_size = 20
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
        f'data/author_data_{date_today}.csv'
    output:
        f'data/author_clean_{date_today}.csv'
    run:
        try:
            send_slack_message(channel,running_msg.format(rule=rule),slack_token=slack_token)
            clean_authors()
            send_slack_message(channel,success_msg.format(rule=rule),slack_token=slack_token)
        except Exception as e:
            send_slack_message(channel,error_msg.format(rule_name=rule,error=e),slack_token=slack_token)
            print(e)
            raise e

rule clean_papers:
    input:
        f'data/paper_data_{date_today}.csv'
    output:
        f'data/paper_clean_{date_today}.csv'
    run:
        try:
            send_slack_message(channel,running_msg.format(rule=rule),slack_token=slack_token)
            clean_papers()
            send_slack_message(channel,success_msg.format(rule=rule),slack_token=slack_token)
        except Exception as e:
            send_slack_message(channel,error_msg.format(rule_name=rule,error=e),slack_token=slack_token)
            print(e)
            raise e

rule filter_authors:
    input:
        f'data/author_clean_{date_today}.csv'
    output:
        f'data/nodes_{{institution}}_{date_today}.csv'
    run:
        try:
            send_slack_message(channel,running_msg.format(rule=rule),slack_token=slack_token)
            filter_authors(wildcards.institution)
            send_slack_message(channel,success_msg.format(rule=rule),slack_token=slack_token)
        except Exception as e:
            rule_name = f"{rule}:({wildards.institution})"
            send_slack_message(channel,error_msg.format(rule_name=rule_name,error=e),slack_token=slack_token)
            print(e)
            raise e

rule filter_papers:
    input:
        f'data/nodes_{{institution}}_{date_today}.csv',
        f'data/paper_clean_{date_today}.csv'
    output:
        f'data/papers_{{institution}}_{date_today}.csv',
        f'data/papers_{{institution}}_2plus_{date_today}.csv'
    run:
        try:
            send_slack_message(channel,running_msg.format(rule=rule),slack_token=slack_token)
            filter_papers(wildcards.institution)
            send_slack_message(channel,success_msg.format(rule=rule),slack_token=slack_token)
        except Exception as e:
            rule_name = f"{rule}:({wildards.institution})"
            send_slack_message(channel,error_msg.format(rule_name=rule_name,error=e),slack_token=slack_token)
            print(e)
            raise e

rule add_publication_stats:
    input: 
        f'data/nodes_{{institution}}_{date_today}.csv',
        f'data/papers_{{institution}}_{date_today}.csv'
    output:
        f'data/nodes_{{institution}}_full_{date_today}.csv'
    run:
        try:
            send_slack_message(channel,running_msg.format(rule=rule),slack_token=slack_token)
            add_publication_stats(wildcards.institution)
            send_slack_message(channel,success_msg.format(rule=rule),slack_token=slack_token)
        except Exception as e:
            rule_name = f"{rule}:({wildards.institution})"
            send_slack_message(channel,error_msg.format(rule_name=rule_name,error=e),slack_token=slack_token)
            print(e)
            raise e

rule create_edges:
    input:
        f'data/papers_{{institution}}_2plus_{date_today}.csv',
        f'data/nodes_{{institution}}_{date_today}.csv'
    output:
        f'data/edges_{{institution}}_{date_today}.csv'
    run:
        try:
            send_slack_message(channel,running_msg.format(rule=rule),slack_token=slack_token)
            create_edgelist(wildcards.institution)
            send_slack_message(channel,success_msg.format(rule=rule),slack_token=slack_token)
        except Exception as e:
            rule_name = f"{rule}:({wildards.institution})"
            send_slack_message(channel,error_msg.format(rule_name=rule_name,error=e),slack_token=slack_token)
            print(e)
            raise e