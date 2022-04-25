from src.process import filter_authors


def main():
    input = snakemake.input[0]  # type: ignore # noqa
    output = snakemake.output[0]  # type: ignore # noqa
    institution = snakemake.wildcards.institution  # type: ignore # noqa

    filter_authors(input, output, institution)


if __name__ == "__main__":
    main()
