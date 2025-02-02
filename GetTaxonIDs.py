import argparse
import MSA


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Taxon species")
parser.add_argument("-n", "--names_file", help="Names file")
args = parser.parse_args()


taxonfile = MSA.SpeciesTaxon(args.file)
taxonfile.get_taxon_ids(args.names_file)

