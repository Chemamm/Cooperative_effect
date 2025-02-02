import argparse
import Combinatory

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Expression matrix")
parser.add_argument("-o", "--output", help="output")
args = parser.parse_args()

em = Combinatory.Ematrix(args.file)
em.add_mean_median_std()
em.df.to_csv(args.output, sep="\t", index=False)
