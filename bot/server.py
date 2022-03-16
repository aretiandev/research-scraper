import os
from flask import Flask, request, Response
from slack_sdk import WebClient
# from slackeventsapi import SlackEventAdapter
from dotenv import load_dotenv
import logging
import subprocess

load_dotenv()
# signing_secret = os.environ['SIGNING_SECRET']
slack_token = os.environ["SLACK_BOT_TOKEN"]
client = WebClient(token=slack_token)
# BOT_ID = client.api_call("auth.test")["user_id"]

app = Flask(__name__)
# slack_event_adapter = SlackEventAdapter(signing_secret, '/slack/events', app)

gunicorn_error_logger = logging.getLogger('gunicorn.error')
app.logger.handlers.extend(gunicorn_error_logger.handlers)
app.logger.setLevel(logging.DEBUG)
app.logger.debug('this will show in the log')


@app.route("/checkscrape", methods=["POST"])
def checkscrape():
    data = request.form
    channel_id = data.get('channel_id')
    n_lines = 8
    if data.get('text'):
        n_lines = int(data.get('text'))
    result = subprocess.run(["/usr/bin/tmux", "capture-pane", "-pt", "SNAKEMAKE"], capture_output=True)
    result_text = result.stdout.decode("utf-8")
    result_chunk = "\n".join(result_text.splitlines()[-n_lines:])
    client.chat_postMessage(channel=channel_id, text=result_chunk)
    return Response(), 200


# Attempt to run snakemake using Slack command
# Problems:
# 1. Snakemake needs to run inside container, so I need to modify the code below
# 2. Command /checkstatus does not work if snakemake is run via Slack command, bcause it reads tmux session. I need to modify /checkstatus so that it reads a log file.
# @app.route("/snakemake", methods=["POST"])
# def snakemake():
#     data = request.form
#     channel_id = data.get('channel_id')
#     n_cores = '16'
#     if data.get('text'):
#         n_cores = data.get('text')
#     try:
#         subprocess.run(["/home/fer/miniconda3/envs/analytics/bin/snakemake", "--cores", n_cores])
#         ["docker", "exec", "barcelona", "snakemake --cores 16
#     except Exception as e:
#         error_msg = f"There was an Exception while attempting to run snakemake:\n{e}"
#         client.chat_postMessage(channel=channel_id, text=error_msg)
#         raise e
#     # result_text = result.stdout.decode("utf-8")
#     # result_chunk = "\n".join(result_text.splitlines()[-n_lines:])
#     return Response(), 200

# @slack_event_adapter.on('message')
# def message(payload):
#     event = payload.get('event', {})
#     channel_id = event.get('channel')
#     user_id = event.get('user')
#     text = event.get('text')
#     # print("Payload:")
#     # app.logger.debug(payload)
#     # print("BOT_ID:")
#     # app.logger.debug(BOT_ID)
#     if BOT_ID != user_id:
#         # client.chat_postMessage(channel=channel_id, text=text)
#         if text[:15] == "check snakemake":
#             n_lines = 8
#             try:
#                 n_lines = int(text.split(":")[1])
#             except Exception:
#                 pass
#             text = text.split(":")[0]

#             result = subprocess.run(["/usr/bin/tmux", "capture-pane", "-pt", "SNAKEMAKE"], capture_output=True)

#             result_text = result.stdout.decode("utf-8")
#             result_chunk = "\n".join(result_text.splitlines()[-n_lines:])
            
#             app.logger.debug(result_chunk)
#             # app.logger.debug(len(result_text))
#             # app.logger.debug(result_text.splitlines()[-20:])

#             client.chat_postMessage(channel=channel_id, text=result_chunk)


@app.route("/")
def hello():
    return "<h1 style='color:blue'>Hello There!</h1>"
