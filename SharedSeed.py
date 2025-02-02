#!/usr/bin/env/python3

import argparse
import Combinatory


parser = argparse.ArgumentParser()
parser.add_argument("-f","--file", help="miRNA host file")
parser.add_argument("-m","--mirna", help="miRNA sequence")
args = parser.parse_args()

Combinatory.shared_seed(args.mirna, args.file)