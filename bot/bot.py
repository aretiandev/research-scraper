#!/usr/bin/env python
"""
Python Aretian Bot

Run with: python -u bot.py | tee log.txt
Then log out of the container with CtrlP+CtrlQ
"""
from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError
from dotenv import load_dotenv
from datetime import datetime
import pytz
import random
import requests
import time
import subprocess
import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))

from src.logging import create_logger


log = create_logger(__name__, f"log/{__name__}.log")


def ping_and_wait(url, status_code=None, wait_time=300, notifications=False):
    """Ping and wait

    Args:
        url: url to ping.
        wait_time: number of seconds to wait. Default: 5 mins.
    """
    tz_NY = pytz.timezone('America/New_York')

    log.info(f"Pinging URL: {url}")
    while True:
        try:
            r = requests.get(url)

            now = datetime.now(tz_NY)
            dt_string = now.strftime("%Y-%m-%d - %H:%M:%S")

            status_history =""
            if r.history:
                status_history = " ".join(str(s.status_code) for s in r.history) + " ->"

            log.info(f"{dt_string} - Ping! Response: {status_history} {r.status_code}: {r.url}")
        except Exception:
            log.exception("There was an exception in the response:")
            pass

        # if r.url != 'https://portalrecerca.manteniment.csuc.cat':
            # break

        # Break when response differs from status_code provided by user
        if status_code:
            if r.status_code != status_code:
                break

        # Always break when 200 OK and no redirects
        if r.status_code == 200:
            if not r.history:  # there were no redirects
                break

        if notifications:
            load_dotenv()
            slack_token = os.environ["SLACK_BOT_TOKEN"]
            slack_member_id = os.environ["SLACK_MEMBER_ID"]
            channel = slack_member_id
            message = ":hourglass: Portal de la Reserca is still down."
            send_slack_message(channel, message, slack_token)

        time.sleep(wait_time - 30 + random.random()*60)

    log.info("Success conditions satisfied, website seems to be back online. Breaking from loop.")


def send_slack_message(channel, message, slack_token):
    log.info(f"Sending Slack message: '{message}' (channel: '{channel}')")

    # Add date and time to message
    tz_NY = pytz.timezone('America/New_York') 
    now = datetime.now(tz_NY)
    dt_string = now.strftime("%m-%d %H:%M:%S")
    date_message = dt_string + ' - ' + message

    client = WebClient(token=slack_token)
    try:
        client.chat_postMessage(
            channel=channel,
            text=date_message
        )
    except SlackApiError as e:
        # You will get a SlackApiError if "ok" is False
        assert e.response["error"]    # str like 'invalid_auth', 'channel_not_found'


if __name__ == "__main__":
    load_dotenv()
    slack_token = os.environ["SLACK_BOT_TOKEN"]
    slack_member_id = os.environ["SLACK_MEMBER_ID"]

    # channel = "#proj-city_barcelona"
    channel = slack_member_id
    url = 'https://portalrecerca.csuc.cat/'

    message = ":bell: Pinging Portal de la Reserca. Will notify when the website is back online: https://portalrecerca.csuc.cat"
    send_slack_message(channel, message, slack_token)

    ping_and_wait(url)

    message = ":globe_with_meridians: Portal de la Reserca is back online: https://portalrecerca.csuc.cat"
    send_slack_message(channel, message, slack_token)

    message = ":snake: Running Snakemake from scratch with 16 cores."
    send_slack_message(channel, message, slack_token)

    smk_command = ["snakemake", "--cores", "all", "--log-handler-script", "src/log_handler_script.py"]
    subprocess.run(smk_command)
