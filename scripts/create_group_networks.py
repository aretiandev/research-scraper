from src.groups import create_group_networks


def main():
    input_nodes = snakemake.input[0]
    input_groups = snakemake.input[1]
    input_edges = snakemake.input[2]
    output_nodes = snakemake.output[0]
    output_edges = snakemake.output[1]

    create_group_networks(input_nodes, input_groups, input_edges,
                          output_nodes, output_edges)


if __name__ == "__main__":
    main()
