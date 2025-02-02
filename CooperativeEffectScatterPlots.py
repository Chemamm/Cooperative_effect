
import argparse
import pandas as pd
import Cooperative

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Positional Consensus file")
parser.add_argument("-o", "--output_folder", help="Output folder")
parser.add_argument("-n", "--number", type=int, default=20, help="Number of miRNAs for the comparison")
args = parser.parse_args()



def main():

    # Read the cooperative_targets.tsv file into a DataFrame
    df = pd.read_csv(args.file, sep="\t")

    # Plot scatterplots with regression lines for each combination
    Cooperative.plot_scatter_regression(df, args.output_folder)

    # Plot comparison of top and bottom miRNAs based on median expression
    Cooperative.plot_comparison_top_bottom(df, args.output_folder, args.number)

    # Plot scatterplots with regression lines for each combination
    Cooperative.plot_scatter_regression_raw(df, args.output_folder)

    # Plot comparison of top and bottom miRNAs based on median expression
    Cooperative.plot_comparison_top_bottom_raw(df, args.output_folder, args.number)

if __name__ == "__main__":
    main()
