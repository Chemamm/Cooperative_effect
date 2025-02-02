import argparse
import MSA

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Problem file")
parser.add_argument("-o", "--output", help="Output")
args = parser.parse_args()


def main():
    mirnaannot = MSA.miRNAannot(args.file)
    mirnaannot.write_df(args.output)

if __name__ == "__main__":
    main()