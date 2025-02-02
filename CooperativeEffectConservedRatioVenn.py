import argparse
import Cooperative


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file_list", help="comma separated list of cooperative_targets.tsv")
parser.add_argument("-s", "--species", help="comma separated list of species")
parser.add_argument("-o", "--output", help="Output")
args = parser.parse_args()


def main():
    Cooperative.upset_top_ratio_mirnas(args.file_list, args.species, args.output)

if __name__ == "__main__":
    main()

