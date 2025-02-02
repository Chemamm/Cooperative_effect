import argparse
import MSA


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Poscon file")
parser.add_argument("-c", "--conservation", help="Conservation Filtered homolog file")
parser.add_argument("-r", "--random", help="Poscon random file")
parser.add_argument("-cr", "--conservation_random", help="Conservation Filtered homolog random file")
parser.add_argument("-t", "--taxon_file", default=False, help="JSON file")
parser.add_argument("-o", "--output", help="Output filename")
args = parser.parse_args()

MSA.plot_gene_stats(args.conservation, args.taxon_file, args.output)
MSA.plot_genes_stats_odds_ratio(args.file, args.conservation, args.random, args.conservation_random, args.output)

