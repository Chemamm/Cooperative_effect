
import argparse
import Combinatory

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--conservation", help="merged orthology file")
parser.add_argument("-n", "--number", type=int, help="number of species to be considered conserved file")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()

con = Combinatory.OrthologyDataframe(args.conservation)
con_filtered = con.remove_non_conserved(args.number)
con_filtered.to_csv(args.output, sep="\t")