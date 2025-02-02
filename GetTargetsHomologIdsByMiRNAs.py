import argparse
import Combinatory

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="list of Positional Consensus files separated by comma")
parser.add_argument("-s", "--species", help="list of species separated by comma")
parser.add_argument("-o", "--output", help="output_label")
args = parser.parse_args()


lof = args.file.split(",")
los = args.species.split(",")

Combinatory.get_targets_with_orthologs_by_mirna(lof, los, args.output)