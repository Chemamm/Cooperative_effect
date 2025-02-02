
import argparse
import Combinatory

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Positional consensus file")
parser.add_argument("-o", "--output", help="output filename")
parser.add_argument("-m", "--metric", help="metric: lift")
parser.add_argument("-p", "--poscon", help="Positional Consensus")
args = parser.parse_args()

Combinatory.apriori_best_groups(args.file, args.metric, args.output, args.poscon)