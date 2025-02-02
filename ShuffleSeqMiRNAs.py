import argparse
import MSA


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Mature miRNAs file")
parser.add_argument("-mg", "--mirgenedb", help="mirgenedb mature mirnas fa")
parser.add_argument("-mb", "--mirbase", help="mirbase mature mirnas fa")
parser.add_argument("-o", "--output", help="Output filename")
args = parser.parse_args()


known_seeds = MSA.get_known_seeds(args.mirgenedb, args.mirbase)
MSA.filter_and_shuffle_miRNAs(args.file, args.output, known_seeds)

