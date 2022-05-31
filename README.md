# Portal de la Reserca Scraping Project

This project scrapes Portal de la Reserca. 

## Setup
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

## Slack Integration
Run:
```
snakemake --cores all --log-handler-script scripts/log_handler.py
```

You will need to put your Slack tokens in an `.env` file in order for the Slack integration to work.

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