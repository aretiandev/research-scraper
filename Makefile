# Snakemake setup
# Sets up snakemake environment
#
# Run `make install` to set up.

.PHONY: help attach install install_packages snakemake help

help: ## View help
	@awk 'BEGIN {FS="^#+ ?"; header=1; body=0}; \
		  header == 1 {printf "\033[36m%s\033[0m\n", $$2} \
		  /^#\s*$$/ {header=0; body=1; next} \
		  body == 1 && /^#+ ?[^ \t]/ {print $$2} \
		  body == 1 && /^#+( {2,}| ?\t)/ {printf "\033[0;37m%s\033[0m\n", $$2} \
		  /^\s*$$/ {print "";exit}' $(MAKEFILE_LIST)
	@echo "Rules:"
	@grep -E '^[a-zA-Z_-]+:.*##[ \t]+.*$$' $(MAKEFILE_LIST) \
	| sort \
	| awk 'BEGIN {FS=":.*##[ \t]+"}; {printf "\033[36m%-20s\033[0m%s\n", $$1, $$2}'

attach: ## Attach to running container
	docker exec -it --user jovyan portalrecerca zsh

install: ## Run docker container and install required packages
	@if [ "$$(docker ps -aq -f name=portalrecerca)" ]; then \
		echo "Container portalrecerca already exists.";   \
		exit 1; \
	fi
	docker compose up -d
	docker exec --user jovyan portalrecerca make install_packages
	$(MAKE) attach

install_packages: ## Install required packages inside container
# https://github.com/mamba-org/mamba/issues/633#issuecomment-812272143
	mamba env update --name base --file environment.yml

snakemake: ## Run all Snakemake rules
	snakemake --cores all

test: ## Run tests
	cd scripts; python -m unittest -b;

