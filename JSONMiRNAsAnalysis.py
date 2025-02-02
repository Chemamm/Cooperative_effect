import argparse
import MSA


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="JSON file")
parser.add_argument("-t", "--taxon_file", default=False, help="JSON file")
parser.add_argument("-o", "--output", help="Output filename")
args = parser.parse_args()


json_miRNAs = MSA.JSONmiRNA(args.file)
json_miRNAs.miRNAs_analysis(args.taxon_file, args.output)

