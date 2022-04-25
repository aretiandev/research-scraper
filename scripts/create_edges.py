from src.edges import create_edges


def main():
    input_authors = snakemake.input[0]
    input_papers = snakemake.input[1]
    output = snakemake.output[0]
    institution = snakemake.wildcards.institution

    create_edges(input_authors, input_papers, output, institution)


if __name__ == "__main__":
    main()
