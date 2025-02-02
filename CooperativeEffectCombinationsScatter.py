import argparse
import Cooperative


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Combination file")
parser.add_argument("-e", "--ematrix", default=False, help="Expression matrix")
parser.add_argument("-o", "--output", help="Output")
args = parser.parse_args()

Cooperative.combinations_scatter_expression(args.file, args.ematrix, args.output)

