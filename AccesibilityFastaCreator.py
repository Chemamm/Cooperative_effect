import argparse
import miRNAtools

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Positional Consensus file")
parser.add_argument("-fa", "--fasta", help="UTR fasta file")
parser.add_argument("-mi", "--miRNAs", help="miRNAs fasta file")
parser.add_argument("-o", "--output", help="Output")
args = parser.parse_args()


def main():
    miRNAtools.process_files(args.fasta, args.miRNA, args.file, agrs.output )

if __name__ == "__main__":
    main()