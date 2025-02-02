import argparse
import Cooperative


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Combination file")
parser.add_argument("-o", "--output", help="Output")
args = parser.parse_args()

Cooperative.combinations_gene_analysis(args.file, args.output)

