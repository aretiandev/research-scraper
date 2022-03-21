from src.process import filter_papers


def main():
    input_authors = snakemake.input[0]
    input_papers = snakemake.input[1]
    output_papers = snakemake.output[0]
    output_papers_2plus = snakemake.output[1]
    institution = snakemake.wildcards.institution

    filter_papers(input_authors, input_papers, output_papers, output_papers_2plus, institution)


if __name__ == "__main__":
    main()
