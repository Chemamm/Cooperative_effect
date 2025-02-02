
import argparse
import Combinatory

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Positional Consensus file")
parser.add_argument("-u", "--utr", help="UTR Fasta file")
parser.add_argument("-c", "--conservation", help="merged orthology file")
parser.add_argument("-n", "--number", type=int, help="number of species to be considered conserved file")
parser.add_argument("-o", "--output", help="output label")
args = parser.parse_args()

Combinatory.randomization_background(args.file, args.utr, args.output)
Combinatory.randomization(args.file, args.output)
Combinatory.randomization_conserved_background(args.file, args.utr, args.conservation, args.number, args.output)