
import argparse
import Combinatory

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Positional Consensus file")
parser.add_argument("-o", "--output", help="output label")
args = parser.parse_args()

poscon = Combinatory.PositionalConsensus(args.file)
poscon.relative_intersection(args.output)


