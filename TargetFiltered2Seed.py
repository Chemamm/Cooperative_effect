import argparse
import MSA


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Combination file")
parser.add_argument("-s", "--seed", help="Seed conversion tsv")
parser.add_argument("-o", "--output", help="Output")
args = parser.parse_args()

MSA.target_filtered_to_seed_network(args.file, args.seed, args.output)
