import argparse
import MSA


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Positional file")
parser.add_argument("-s", "--structure", help="Positional file")
parser.add_argument("-c", "--chrAlias", help="Positional file")
parser.add_argument("-u", "--utr_only",action="store_true", default=False, help="Positional file")
parser.add_argument("-o", "--output", help="Output")
args = parser.parse_args()


poscon = MSA.PositionalConsensus(args.file)
poscon.add_genomic_coordinates(args.structure, args.chrAlias, args.utr_only, args.output)

