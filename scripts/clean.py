from src.process import clean


def main():
    items = snakemake.wildcards.item_name   # type: ignore # noqa
    input = snakemake.input[0]              # type: ignore # noqa
    output = snakemake.output[0]            # type: ignore # noqa

    clean(items, input, output)


if __name__ == "__main__":
    main()
