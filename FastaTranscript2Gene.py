
import argparse
import Combinatory

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Fasta file")
parser.add_argument("-m", "--mapping", help="Mapping tsv file")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()

Combinatory.transcript_to_gene(args.file, args.mapping, args.output)