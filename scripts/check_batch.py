"""
Sends Slack message when given batch number is completed.
"""
from slack_sdk import WebClient
from dotenv import load_dotenv
import subprocess
import time
import os

load_dotenv()
# signing_secret = os.environ['SIGNING_SECRET']
slack_token = os.environ["SLACK_BOT_TOKEN"]
channel_id = os.environ["SLACK_MEMBER_ID"]

client = WebClient(token=slack_token)

n_lines = 1
n_batch = 0

while n_batch < 492:
    result = subprocess.run(
        ["/usr/bin/tmux", "capture-pane", "-pt", "SNAKEMAKE"], capture_output=True
    )
    result_text = result.stdout.decode("utf-8")
    result_chunk = "\n".join(result_text.splitlines()[-n_lines:])
    n_batch = result_chunk.split(")")[0].split("(")[1].split("/")[0]
    n_batch = int(n_batch)
    print(f"Batch number: {n_batch}.")
    time.sleep(300)

client.chat_postMessage(
    channel=channel_id, text="DONE! Scraping paper_data has reached batch #492."
)
