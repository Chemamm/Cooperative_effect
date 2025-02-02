import argparse
import Cooperative

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Problem file")
parser.add_argument("-r", "--random", help="Randomized file")
parser.add_argument("-o", "--output", help="Output")
args = parser.parse_args()


def main():
    Cooperative.spearman_rank_analysis(args.file, args.random, args.output)
if __name__ == "__main__":
    main()