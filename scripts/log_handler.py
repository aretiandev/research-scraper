from slack_sdk import WebClient
from dotenv import load_dotenv
import os
load_dotenv()

slack_token = os.environ["SLACK_BOT_TOKEN"]
slack_member_id = os.environ["SLACK_MEMBER_ID"]
client = WebClient(token=slack_token)


def log_handler(msg):
    msg_level = msg.get('level')
    if msg_level is None:
        return

    if msg_level == 'run_info':
        message = ':snake: ' + "```" + msg['msg'] + "```"
    elif msg_level == 'job_info':
        message_lst = [
            f":runner: Running rule '{msg['name']}':",
            f"    jobid: {msg['jobid']}",
            f"    input: {msg['input']}" if msg['input'] else "",
            f"    output: {msg['output']}" if msg['output'] else "",
            f"    wildcards: {msg['wildcards']}" if msg['wildcards'] else ""
            ]
        message_lst = [x for x in message_lst if x]
        message = "\n".join(message_lst)
        message = "```" + message + "```"
    elif msg_level == 'job_finished':
        message = f":white_check_mark: Finished job {msg['jobid']}."
    elif msg_level == 'progress':
        message = f":construction: {msg['done']} of {msg['total']} steps ({msg['done']/msg['total']*100:.0f}%) done"
    elif msg_level == 'job_error':
        message = f"Error in rule `{msg['name']}`: jobid {msg['jobid']}"
    elif msg_level == 'error':
        message = msg['msg']
        message = "```" + message + "```"
    else:
        return

    client.chat_postMessage(
        channel=slack_member_id,
        text=message
        )
