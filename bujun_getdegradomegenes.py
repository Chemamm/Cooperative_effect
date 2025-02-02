#!/usr/bin/env/python3

import Combinatory
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="degradome file")
parser.add_argument("-i", "--input_list", help="input list")
parser.add_argument("-m", "--mature_file", help="mature file to complete the input list")
parser.add_argument("-n1", "--n1", type=int, help="cleavage site position 1")
parser.add_argument("-n2", "--n2", type=int, help="cleavage site position 2")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()

c_file = Combinatory.Degradome(args.file)
c_file.get_genes(args.input_list, args.mature_file, args.n1,args.n2,args.output)