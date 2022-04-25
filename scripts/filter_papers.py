from src.process import filter_papers


def main():
    input_authors = snakemake.input[0]  # type: ignore # noqa
    input_papers = snakemake.input[1]  # type: ignore # noqa
    output_papers = snakemake.output[0]  # type: ignore # noqa
    output_papers_2plus = snakemake.output[1]  # type: ignore # noqa
    institution = snakemake.wildcards.institution  # type: ignore # noqa

    filter_papers(
        input_authors, input_papers, output_papers, output_papers_2plus, institution
    )


if __name__ == "__main__":
    main()
