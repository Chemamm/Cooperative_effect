import argparse
import Cooperative

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Positional Consensus file")
parser.add_argument("-s", "--start", type=int, default=13, help="Relative start position of the cooperative effect region")
parser.add_argument("-e", "--end", type=int, default=39, help="Relative end position of the cooperative effect region")
parser.add_argument("-m", "--matrix", default=False, help="Expression matrix")
parser.add_argument("-o", "--output_folder", help="Output folder")
args = parser.parse_args()


def main():
    # Check if the positional consensus file is provided
    if not args.file:
        print("Error: Please provide a positional consensus file using -f or --file option.")
        return

    # Create an instance of PositionalConsensus
    positional_consensus = Cooperative.PositionalConsensus(args.file)

    # Sort the positional consensus data by gene
    positional_consensus.sort_by_gene()

    # Find adjacent region clusters based on the specified start and end parameters
    region_clusters = positional_consensus.find_adjacent_region_clusters(args.start, args.end)

    # Process region clusters to calculate valid miRNA combinations and counts
    region_clusters, valid_combos = Cooperative.process_region_clusters(region_clusters, args.start, args.end)

    if args.matrix:
        region_clusters = Cooperative.calculate_median_expression(args.matrix, region_clusters)
        valid_combos = Cooperative.calculate_median_expression(args.matrix, valid_combos)

        # Write the final TSV output
        Cooperative.write_final_tsv(region_clusters, args.output_folder)
        Cooperative.write_combo_tsv(valid_combos, args.output_folder)
    else:
        Cooperative.write_final_tsv_noexp(region_clusters, args.output_folder)
        Cooperative.write_combo_tsv_noexp(valid_combos, args.output_folder)


if __name__ == "__main__":
    main()