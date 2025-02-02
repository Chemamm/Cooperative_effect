import argparse
import Combinatory

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--folder", help="Folder output of intersecter")
parser.add_argument("-s", "--stats", help="stats file")
args = parser.parse_args()

Combinatory.add_zscore_to_stats(args.folder, args.stats)