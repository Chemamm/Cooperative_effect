import argparse
import MSA


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="JSON file")
parser.add_argument("-c", "--clusters", type=int, help="Number of clusters for Kmeans")
parser.add_argument("-o", "--output", help="Output filename")
args = parser.parse_args()


json_miRNAs = MSA.JSONmiRNA(args.file)
json_miRNAs.unsupervised_analysis(args.clusters, args.output)

