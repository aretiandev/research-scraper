from src.process import clean


def main():
    items = snakemake.wildcards.item_name
    input = snakemake.input[0]
    output = snakemake.output[0]

    print()
    print("Snakemake parameters:")
    print(items)
    print(input)
    print(output)
    print()

    clean(items, input, output)


if __name__ == "__main__":
    main()
