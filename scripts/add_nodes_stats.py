from src.process import add_nodes_stats


def main():
    input_authors = snakemake.input[0]
    input_papers = snakemake.input[1]
    input_groups = snakemake.input[2]
    output = snakemake.output[0]
    institution = snakemake.wildcards.institution

    add_nodes_stats(input_authors, input_papers, input_groups, output, institution)


if __name__ == "__main__":
    main()
