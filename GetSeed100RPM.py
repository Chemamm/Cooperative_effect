#!/usr/bin/env/python3

import argparse
import Combinatory

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="mirna data file")
parser.add_argument("-s", "--seed", help="seed annot data file")
parser.add_argument("-o", "--output", help="output label")
args = parser.parse_args()


Combinatory.get_seed_100rpm(args.file, args.seed, args.output)
