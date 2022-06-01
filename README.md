# Portal de la Reserca Scraping Project

This project scrapes Portal de la Reserca. 

# Instructions

## Configuration

Edit `config.py`:
- DATE: current date for versioning files and folders.
- INSTITUTIONS_LIST: list of institutions to scrape.
- THREADS_MAX: Snakemake parameter. Recommended: max number of CPU cores.
- TIMEOUT: timeout in seconds between batches avoid crashing website.
- SLACK_BOT_TOKEN: Token for notifications.
- SLACK_MEMBER_ID: Member ID for notifications.

Set current date (if not set will default to today):
```
export DATE=20220601
```

## Create Container

Create and attach to JupyterLab container with required packages and open port 8888.
```
docker compose up -d
docker exec -it --user jovyan portalrecerca zsh
```
Or `make run`.

## Run Snakemake

From within the container, check everything is installed properly:
```
snakemake --cores all --dry-run
```

Run full pipeline:
```
snakemake --cores all
```

The project is self documented in the `Snakefile`.

All output is saved into `recerca.db` sqlite database and csv files in `data/[date_today]`.

# Other features

## Slack Integration

Run:
```
snakemake --cores all --log-handler-script scripts/slack.py
```

You will need to put `SLACK_BOT_TOKEN` and `SLACK_MEMBER_ID` in an `.env` file in the root directory.

## Run tests

All tests are run from the `scripts` folder:
```
cd scripts
```

Create intermediate files for tests with:
```
python -m tests pickle
```

Run tests:
```
python -m unittest
```

Silence output while running tests with `python -m unittest -b`.

# Notes

When working from inside the container, you can git add and commit changes but not push them to the remote, since the container does not have the necessary keys.