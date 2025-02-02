#!/usr/bin/env/python3

import argparse
import Combinatory


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="positional consensus file")
parser.add_argument("-l", "--list", help="mirna list file")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()


pc = Combinatory.PositionalConsensus(args.file)
pc.get_num_targets(args.list, args.output)