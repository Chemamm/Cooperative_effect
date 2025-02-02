import argparse
import Combinatory

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Fasta")
args = parser.parse_args()

fa = Combinatory.Fasta(args.file)
fa.get_mean_length()
