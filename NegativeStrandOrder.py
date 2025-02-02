import argparse
import MSA


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Positional file")
parser.add_argument("-o", "--output", help="Output")
args = parser.parse_args()


genomic_structure = MSA.GenomicStructure(args.file)
df = genomic_structure.parse_strand_negtive_coords()
df.to_csv(args.output, sep="\t", index=False)

