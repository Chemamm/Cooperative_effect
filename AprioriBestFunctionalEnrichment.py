
import argparse
import Combinatory

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Positional Consensus file")
parser.add_argument("-o", "--out_label", help="output label")
args = parser.parse_args()

Combinatory.lift_best_functional_enrichment(args.file, args.out_label)


