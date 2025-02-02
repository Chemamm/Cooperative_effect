#!/usr/bin/env/python3

import Combinatory
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="enrichment file")
parser.add_argument("-go", "--GOobo", default=False, help="GO obo descriptions file")
parser.add_argument("-k", "--KEGobo", default=False, help="KEG descriptions file")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()

e_file = Combinatory.Enrichment(args.file)
if args.GOobo:
    e_file.depurate(args.GOobo, args.output)
elif args.KEGobo:
    e_file.depurateKEG(args.KEGobo, args.output)

