import argparse
import Combinatory

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--folder", help="Folder with the files")
parser.add_argument("-s", "--species", help="Species identifiers to analyze sepparated by comma")
args = parser.parse_args()

Combinatory.merge_UTR_flanking_CDS(args.folder, args.species)

