from src.groups import filter_group_edges

def main():
    input_group_data = snakemake.input[0]
    input_group_edges = snakemake.input[1]
    output_internal_gp_edges = snakemake.output[0]

    filter_group_edges(
        input_group_data, input_group_edges, output_internal_gp_edges
    )


if __name__ == "__main__":
    main()
