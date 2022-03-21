import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))

from src.edges import create_edgelist


def main():
    input_authors = snakemake.input[0]
    input_papers = snakemake.input[1]
    output = snakemake.output[0]
    institution = snakemake.wildcards.institution

    create_edgelist(input_authors, input_papers, output, institution)


if __name__ == "__main__":
    main()
