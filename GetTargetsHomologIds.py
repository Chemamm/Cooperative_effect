import argparse
import Combinatory

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Positional Consensus file")
parser.add_argument("-c", "--conservation", help="merged orthology file")
parser.add_argument("-o", "--output", help="output")
args = parser.parse_args()

Combinatory.get_targets_with_orthologs(args.file, args.conservation, args.output)

