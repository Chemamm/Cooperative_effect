import argparse
import Cooperative

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Conserved poscon")
parser.add_argument("-c", "--coop", help="Valid cooperations file")
parser.add_argument("-r", "--random", help="Valid cooperations randomized file")
parser.add_argument("-m", "--min_coop", type=int, help="Minimum number of cooperations to enter the analysis.")
parser.add_argument("-o", "--output", help="Output")
args = parser.parse_args()


def main():
    Cooperative.significance_numcoop(args.file, args.coop, args.random, args.output, args.min_coop)

if __name__ == "__main__":
    main()