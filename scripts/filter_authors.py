import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))

from src.process import filter_authors


def main():
    input = snakemake.input[0]
    output = snakemake.output[0]
    institution = snakemake.wildcards.institution

    filter_authors(input, output, institution)


if __name__ == "__main__":
    main()
