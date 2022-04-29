# Portal de la Reserca Scraping Project

This project scrapes Portal de la Reserca. 

## Setup
Run Docker container:
`docker compose up -d`

Attach to container:
`docker exec -it --user jovyan portalrecerca zsh`

Install packages:
`mamba env update --name base --file environment.yml`

## Snakefile
The project is self documented in the `Snakefile`.

Run everything with the following command:
`snakemake --cores all`

## Slack Integration
Run:
`snakemake --cores all --log-handler-script scripts/log_handler.py`

You will need to put your Slack tokens in an `.env` file in order for the Slack integration to work.
