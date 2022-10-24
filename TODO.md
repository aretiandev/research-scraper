# TODO


## Update 10/12/2022

I will do one institution at a time.

## Problem

When creating the paper_id column, I am assuming the papers table is ordered by institution. That is not the case. The url_stem column is created in the parse_catalog rule which runs in parallel, meaning urls are not created by institution but in a random order. 

- Solution 1: run an ORDER BY before creating the paper_ids.
- Solution 2: Use a different id (or none) and use Snakemake checkpoints. 

The second alternative is more flexible since I no longer have to keep exact track of the names of files that are being produced, I just care about whatever was produced by the previous rule. Even better I can use the url_stem as file name.

Scraped paper urls are inserted in the same table for all institutions. therefore when I want to create the paper_ids, the id
2022-06-10 18:22:14 parse_catalog INFO: Parsed catalog HTML from EP: output/catalog/EP/8.html
2022-06-10 18:22:14 parse_catalog INFO: Parsed catalog HTML from FM: output/catalog/FM/1.html
2022-06-10 18:22:14 parse_catalog INFO: Parsed catalog HTML from FM: output/catalog/FM/3.html
2022-06-10 18:22:14 parse_catalog INFO: Parsed catalog HTML from EP: output/catalog/EP/1.html
- Use item_id as primary key

  - change variables in sqlite tables
  - use snakemake checkpoints for unexpected data
- Make sure everything is running.

- Move all outputs to a database.
- Finish all networks for all institutions. Then create networks using groups of institutions.

- Extract author affiliation from the field "USP Affiliated Authors" (example: https://repositorio.usp.br/item/003058452)
- Increase the limit of edges in output_edges.py from 10000 to the max 22224.

## Tags
- Create tag table for edges
- Create institution tag table for papers
- Create author-paper table
- Create tag table for edges


# List of institutions

Faculdade de Medicina (FM)

Escola Politécnica (EP)

