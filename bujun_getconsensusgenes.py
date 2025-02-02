#!/usr/bin/env/python3

import Combinatory
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="consensus file")
parser.add_argument("-i", "--input_list", help="input list")
parser.add_argument("-n", "--nmethods", type=int, help="number of minimum methods to continue the analysis")
parser.add_argument("-m", "--mature_file", help="mature file to complete the input list")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()

c_file = Combinatory.Consensus(args.file)
c_file.get_consensus_genes(args.input_list, args.mature_file, args.nmethods, args.output)
