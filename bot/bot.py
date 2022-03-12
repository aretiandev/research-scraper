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


def ping_and_wait(wait_time):
    """Ping and wait
    
    Args:
        wait_time: number of seconds to wait.
    """
    # Setup
    tz_NY = pytz.timezone('America/New_York') 

    print("Pinging Portal de la Reserca...")
    while True:
        url = 'https://portalrecerca.csuc.cat/'
        try:
            r = requests.get(url)

            now = datetime.now(tz_NY)
            dt_string = now.strftime("%Y-%m-%d - %H:%M:%S")
            print(f"{dt_string} - Ping! Response status: {r.status_code}")
        except Exception:
            pass

        if r.url != 'https://portalrecerca.manteniment.csuc.cat':
            break

        if not r.history[0].is_redirect:
            break

        time.sleep(wait_time - 30 + random.random()*60)

    print("There was and exception. Stop pinging and exit.")


def send_slack_message(channel, message, slack_token):
    print("Sending message to slack.")
    print(f"channel: {channel}.")
    print(f"message: {message}.")

    client = WebClient(token=slack_token)
    try:
        client.chat_postMessage(
            channel=channel,
            text=message
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

    message = ":globe_with_meridians: Seems like Portal de la Reserca is back online: https://portalrecerca.csuc.cat"

    ping_and_wait(300)
    send_slack_message(channel, message, slack_token)
