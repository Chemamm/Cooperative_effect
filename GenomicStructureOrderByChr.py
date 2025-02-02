import argparse
import MSA


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Genomic file")
parser.add_argument("-o", "--output", help="Output")
args = parser.parse_args()


genomic_structure = MSA.GenomicStructure(args.file)
genomic_structure.order_by_chr()
genomic_structure.df.to_csv(args.output, sep="\t", index=False)

