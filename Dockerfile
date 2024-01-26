FROM yufernando/jupyterlab:snakemake

USER root

RUN apt update && apt upgrade -y && apt install -y graphviz

USER jovyan

WORKDIR /home/jovyan/work

COPY --chown=jovyan environment.yml .

RUN mamba env update --name base --file environment.yml