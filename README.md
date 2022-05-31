# Portal de la Reserca Scraping Project

This project scrapes Portal de la Reserca. 

## Setup
Create JupyterLab container with open port 8888 and install required packages.
```
make install
```
## Run Snakemake
Attach to the running container
```
docker exec -it --user jovyan portalrecerca zsh
```
or `make attach`.

Check everything is installed properly:
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

