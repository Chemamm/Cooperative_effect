import argparse
import pandas as pd
import os
import Cooperative

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Positional Consensus file")
parser.add_argument("-c", "--combinations", help="Valid combinations file")
parser.add_argument("-s", "--species", help="Positional consensus files from other species")
parser.add_argument("-m", "--matrix", help="Expression matrix")
parser.add_argument("-o", "--output_folder", help="Output folder")
parser.add_argument("-t", "--threshold", type=int, help="Threshold of conservation")
args = parser.parse_args()

class PositionalConsensus:
    def __init__(self, file):
        self.df = pd.read_csv(file, sep="\t")

    def sort_by_gene(self):
        self.df.sort_values(by=['transcript', 'end'], inplace=True)

    def find_adjacent_region_clusters(self, start, end):
        region_clusters = []
        current_gene = None
        current_cluster = []

        for index, row in self.df.iterrows():
            # If current gene changes, finalize the current cluster
            if row['transcript'] != current_gene:
                if current_cluster:
                    region_clusters.append({'transcript': current_gene, 'regions': current_cluster})
                current_gene = row['transcript']
                current_cluster = []

            # Check for adjacent regions
            if current_cluster:
                last_region = current_cluster[-1]
                if (row['end'] - last_region['end']) <= end:
                    current_cluster.append({'start': row['end'] - 6, 'end': row['end'], 'miRNA': row['miRNA']})
                else:
                    # Finalize current cluster
                    if current_cluster:
                        region_clusters.append({'transcript': current_gene, 'regions': current_cluster})
                    current_cluster = [{'start': row['end'] - 6, 'end': row['end'], 'miRNA': row['miRNA']}]
            else:
                current_cluster.append({'start': row['end'] - 6, 'end': row['end'], 'miRNA': row['miRNA']})

        # Finalize the last cluster for the last gene
        if current_cluster:
            region_clusters.append({'transcript': current_gene, 'regions': current_cluster})

        return region_clusters

    def get_dict(self):
        poscon_dict = {}
        count_dict = {}
        for index,row in self.df.iterrows():
            if "|" in row['transcript']:
                row['transcript'] = row['transcript'].split("|")[0]
            if row['miRNA'] not in poscon_dict:
                poscon_dict[row['miRNA']] = [row['transcript']]
                count_dict[row['miRNA'] + ";" + row['transcript']] = 1
            else:
                if row['transcript'] not in poscon_dict[row['miRNA']]:
                    poscon_dict[row['miRNA']].append(row['transcript'])
                    count_dict[row['miRNA'] + ";" + row['transcript']] = 1
                else:
                    poscon_dict[row['miRNA']].append(row['transcript'] + str(count_dict[row['miRNA'] + ";" + row['transcript']]))
                    count_dict[row['miRNA'] + ";" + row['transcript']] += 1
        return poscon_dict


class Ematrix:
    def __init__(self,file):
        self.df = pd.read_csv(file, sep="\t")

    def add_mean_median_std(self):
        self.df['mean'] = self.df.iloc[:, 1:].mean(axis=1)
        self.df['median'] = self.df.iloc[:, 1:].median(axis=1)
        self.df['std_dev'] = self.df.iloc[:, 1:].std(axis=1)
        # Add rank columns based on mean and median expression
        self.df['mean_rank'] = self.df['mean'].rank(ascending=False)
        self.df['median_rank'] = self.df['median'].rank(ascending=False)

        # Sort the DataFrame by median expression in descending order
        self.df = self.df.sort_values(by='median', ascending=False)


def write_tsv(df, output_folder):
    # Create the output folder if it doesn't exist
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)

    tsv_path = os.path.join(output_folder, "cooperative_targets.tsv")
    df.to_csv(tsv_path, index=False, sep="\t")

def main():
    # Check if the positional consensus file is provided
    if not args.file:
        print("Error: Please provide a positional consensus file using -f or --file option.")
        return
    df = Cooperative.get_cooperative_targets(args.file, args.combinations, args.species, args.threshold)
    df = Cooperative.add_expression_exosome(df, args.matrix)
    write_tsv(df, args.output_folder)

if __name__ == "__main__":
    main()
