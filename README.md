# Portal de la Reserca Scraping Project

This project scrapes Portal de la Reserca. 

## Install
Run `make install`.

## Snakemake
Check everything is installed properly:
`snakemake --cores all --dry-run`

Run full pipeline:
`snakemake --cores all`

The project is self documented in the `Snakefile`.

## Slack Integration
Run:
`snakemake --cores all --log-handler-script scripts/log_handler.py`

You will need to put your Slack tokens in an `.env` file in order for the Slack integration to work.
