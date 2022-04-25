from src.process import add_nodes_stats


def main():
    input_authors = snakemake.input[0]  # type: ignore # noqa
    input_papers = snakemake.input[1]  # type: ignore # noqa
    input_groups = snakemake.input[2]  # type: ignore # noqa
    output = snakemake.output[0]  # type: ignore # noqa
    institution = snakemake.wildcards.institution  # type: ignore # noqa

    add_nodes_stats(input_authors, input_papers, input_groups, output, institution)


if __name__ == "__main__":
    main()
