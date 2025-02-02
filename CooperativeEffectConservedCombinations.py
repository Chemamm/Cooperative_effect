import argparse
import Cooperative

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--combinations", help="Valid combinations file")
parser.add_argument("-s", "--species", help="Valid combinations files for other species sepparated by comma" )
parser.add_argument("-o", "--output_folder", help="Output folder")
parser.add_argument("-t", "--threshold", type=int,  default=3, help="Number of species needed to be considered as conserved")
args = parser.parse_args()



def main():
    Cooperative.get_conserved_combinations(args.combinations, args.species, args.threshold, args.output_folder)

if __name__ == "__main__":
    main()








