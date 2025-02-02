import argparse
import Cooperative

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Valid cooperations file")
parser.add_argument("-o", "--output", help="Output")
args = parser.parse_args()


def main():
    Cooperative.hsa_plus_iri_cooperations_analysis(args.file, args.output)
    Cooperative.hsa_plus_iri_cooperations_analysis_iri_activated(args.file, args.output)

if __name__ == "__main__":
    main()