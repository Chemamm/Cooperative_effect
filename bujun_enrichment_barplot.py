#!/usr/bin/env/python3

import Combinatory
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="enrichment file")
parser.add_argument("-o", "--out_label", help="output filename label")
parser.add_argument("-n", "--number", type=int, help="number of bars to plot")
args = parser.parse_args()

en_file = Combinatory.Enrichment(args.file)
en_file.barplot_bujun(args.out_label, args.number)
