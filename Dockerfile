FROM yufernando/jupyterlab:snakemake

USER jovyan

WORKDIR /home/jovyan/work

COPY environment.yml .

RUN mamba env update --name base --file environment.yml