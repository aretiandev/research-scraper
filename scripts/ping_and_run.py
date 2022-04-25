import os
import subprocess
from dotenv import load_dotenv
from src.bot import send_slack_message, ping_and_wait


def main():
    load_dotenv()
    slack_token = os.environ["SLACK_BOT_TOKEN"]
    slack_member_id = os.environ["SLACK_MEMBER_ID"]

    # channel = "#proj-city_barcelona"
    channel = slack_member_id
    url = "https://portalrecerca.csuc.cat/"

    message = (
        ":bell: Pinging Portal de la Reserca. "
        + "Will notify when the website is back online: https://portalrecerca.csuc.cat"
    )
    send_slack_message(channel, message, slack_token)

    ping_and_wait(url)

    message = (
        ":globe_with_meridians: "
        + "Portal de la Reserca is back online: https://portalrecerca.csuc.cat"
    )
    send_slack_message(channel, message, slack_token)

    message = ":snake: Running Snakemake from scratch with 16 cores."
    send_slack_message(channel, message, slack_token)

    smk_command = [
        "snakemake",
        "--cores",
        "all",
        "--log-handler-script",
        "scripts/log_handler.py",
    ]
    subprocess.run(smk_command)


if __name__ == "__main__":
    main()
