#!/usr/bin/env/python3

import argparse
import Combinatory

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--folder", help="Proteome data file")
args = parser.parse_args()

Combinatory.orthology_merge(args.folder)





