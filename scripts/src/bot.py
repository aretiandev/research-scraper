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
import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
from src.logging import create_logger


log = create_logger(__name__, f"log/{__name__}.log")


def ping_and_wait(url=None, status_code=None, wait_time=300, notifications=None, slack_msg=None):
    """Ping and wait

    Args:
        url: url to ping.
        wait_time: number of seconds to wait. Default: 5 mins.
    """
    tz_NY = pytz.timezone('America/New_York')

    log.info(f"Pinging URL: {url}")

    if slack_msg is not None:
        if notifications is None:
            notifications = 1

    if notifications is not None:
        load_dotenv()
        slack_token = os.environ["SLACK_BOT_TOKEN"]
        slack_member_id = os.environ["SLACK_MEMBER_ID"]
        channel = slack_member_id
        message = f":warning: Website is down. Pinging URL: {url}."
        if slack_msg is not None:
            message = slack_msg
        send_slack_message(channel, message, slack_token)

    while True:
        try:
            r = requests.get(url)

            now = datetime.now(tz_NY)
            dt_string = now.strftime("%Y-%m-%d - %H:%M:%S")

            status_history = ""
            if r.history:
                status_history = " ".join(str(s.status_code) for s in r.history) + " -> "

            log.info(f"{dt_string} - Ping! Response: {status_history}{r.status_code}: {r.url}")
        except Exception:
            log.exception("There was an exception in the response:")
            pass

        if url is None:
            if r.url != 'https://portalrecerca.manteniment.csuc.cat':
                break
        
        # Break when response differs from status_code provided by user
        if status_code is not None:
            if r.status_code != status_code:
                break

        # Always break when 200 OK and no redirects
        if r.status_code == 200:
            if not r.history:  # there were no redirects
                break

        if notifications == 2:
            message = ":hourglass: Portal de la Reserca is still down."
            send_slack_message(channel, message, slack_token)

        time.sleep(max(wait_time - 30 + random.random()*60, 0))

    if notifications is not None:
        message = ":zap: Website is back online."
        send_slack_message(channel, message, slack_token)

    log.info("Website is back online. Breaking from loop.")


def send_slack_message(channel, message, slack_token, time=False):
    log.info(f"Sending Slack message: '{message}' (channel: '{channel}')")

    # Add date and time to message
    if time:
        tz_NY = pytz.timezone('America/New_York') 
        now = datetime.now(tz_NY)
        dt_string = now.strftime("%m-%d %H:%M:%S")
        message = dt_string + ' - ' + message

    client = WebClient(token=slack_token)
    try:
        client.chat_postMessage(
            channel=channel,
            text=message
        )
    except SlackApiError as e:
        # You will get a SlackApiError if "ok" is False
        assert e.response["error"]    # str like 'invalid_auth', 'channel_not_found'
