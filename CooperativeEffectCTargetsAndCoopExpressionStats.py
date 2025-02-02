import argparse
import Cooperative

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Conserved targets")
parser.add_argument("-c", "--coop",  help="Valid cooperations")
parser.add_argument("-e", "--ematrix", help="Expression matrix")
parser.add_argument("-m", "--miRNAs", help="Saliva-abundant miRNAs")
parser.add_argument("-o", "--output_label", help="Output label")
args = parser.parse_args()

def main():
    Cooperative.expression_ctargets_stats(args.file, args.miRNAs, args.ematrix, args.output_label)
    Cooperative.expression_cooperations_stats(args.coop, args.miRNAs, args.ematrix, args.output_label)

if __name__=="__main__":
    main()