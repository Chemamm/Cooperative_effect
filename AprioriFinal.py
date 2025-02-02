
import argparse
import Combinatory

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Positional consensus file")
parser.add_argument("-o", "--out_label", help="output filename")
parser.add_argument("-ms", "--min_supp", help="minimum support")
parser.add_argument("-lm", "--low_memory",default=True, help="minimum support")
parser.add_argument("-ml", "--memory_limit", default=300, help="memory limit usage")
args = parser.parse_args()

Combinatory.apriori_final(args.file, args.min_supp, args.out_label, args.low_memory, args.memory_limit)