import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Positional Consensus file")
parser.add_argument("-c", "--combinations", help="Valid combinations file")
parser.add_argument("-m", "--matrix", help="Expression matrix")
parser.add_argument("-o", "--output_folder", help="Output folder")
parser.add_argument("-t", "--tmp", default=False, help="Tmp file")
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

class ValidCombination:
    def __init__(self, file):
        self.df = pd.read_csv(file, sep="\t")

    def get_cooperative_targets(self):
        cooperative_targets_2 = []
        cooperative_targets_3 = []

        for _, row in self.df.iterrows():
            if row['Number of targets'] == 2:
                # Split the miRNAs and targets
                miRNAs = row['miRNAs'].split(';')
                targets = row['targets'].split(';')

                # Generate IDs for each miRNA-target pair
                for miRNA, target in zip(miRNAs, targets):
                    gene_id = row['Cooperation ID'].split('-')[0]
                    ID = f"{gene_id};{miRNA};{target}"
                    cooperative_targets_2.append(ID)

            elif row['Number of targets'] == 3:
                # Split the miRNAs and targets
                miRNAs = row['miRNAs'].split(';')
                targets = row['targets'].split(';')

                # Generate IDs for each miRNA-target pair
                for miRNA, target in zip(miRNAs, targets):
                    gene_id = row['Cooperation ID'].split('-')[0]
                    ID = f"{gene_id};{miRNA};{target}"
                    cooperative_targets_3.append(ID)
        return cooperative_targets_2, cooperative_targets_3

def get_cooperative_targets(file, combo):
    poscon = PositionalConsensus(file)
    combo = ValidCombination(combo)

    valid_2, valid_3 = combo.get_cooperative_targets()

    total_targets_count = {}
    targets_valid_2_count = {}
    targets_valid_3_count = {}
    ratio_2 = {}
    ratio_3 = {}

    # Iterate over dataframe rows
    for _, row in poscon.df.iterrows():
        # Generate ID
        gene_id = row['transcript']
        ID = f"{gene_id};{row['miRNA']};{row['start']}:{row['end']}"

        # Check if ID is in valid_combinations_2 or valid_combinations_3
        if ID in valid_2:
            # Increment count in targets_in_valid_2_count dictionary
            targets_valid_2_count[row['miRNA']] = targets_valid_2_count.get(row['miRNA'], 0) + 1
        if ID in valid_3:
            # Increment count in targets_in_valid_3_count dictionary
            targets_valid_3_count[row['miRNA']] = targets_valid_3_count.get(row['miRNA'], 0) + 1

        # Increment count in total_targets_count dictionary
        total_targets_count[row['miRNA']] = total_targets_count.get(row['miRNA'], 0) + 1
    for key in total_targets_count:
        total_targets = total_targets_count[key]
        if key in targets_valid_2_count:
            targets_valid_2 = targets_valid_2_count[key]
        else:
            targets_valid_2_count[key] = 0
            targets_valid_2 = 0
        if key in targets_valid_3_count:
            targets_valid_3 = targets_valid_3_count[key]
        else:
            targets_valid_3_count[key] = 0
            targets_valid_3 = 0

        ratio_2[key] = targets_valid_2 / total_targets
        ratio_3[key] = targets_valid_3 / total_targets
    # Convert dictionaries to DataFrames
    df_total_targets = pd.DataFrame(total_targets_count.items(), columns=['miRNA', 'total_targets_count'])
    df_valid_2_targets = pd.DataFrame(targets_valid_2_count.items(), columns=['miRNA', 'targets_valid_2_count'])
    df_valid_3_targets = pd.DataFrame(targets_valid_3_count.items(), columns=['miRNA', 'targets_valid_3_count'])
    df_ratio_2 = pd.DataFrame(ratio_2.items(), columns=['miRNA', 'ratio_2'])
    df_ratio_3 = pd.DataFrame(ratio_3.items(), columns=['miRNA', 'ratio_3'])

    # Merge DataFrames on 'miRNA'
    df_merged = df_total_targets.merge(df_valid_2_targets, on='miRNA', how='outer')
    df_merged = df_merged.merge(df_valid_3_targets, on='miRNA', how='outer')
    df_merged = df_merged.merge(df_ratio_2, on='miRNA', how='outer')
    df_merged = df_merged.merge(df_ratio_3, on='miRNA', how='outer')
    return df_merged


def add_expression_exosome(df, ematrix):
    ematrix = Ematrix(ematrix)

    # Calculate median of each row in ematrix.df
    ematrix.df['median_expression'] = ematrix.df.iloc[:, 1:].median(axis=1)

    # Calculate medians for different categories and add as new columns
    ematrix.df['median_10K'] = ematrix.df.filter(like='10K').median(axis=1)
    ematrix.df['median_150K'] = ematrix.df.filter(like='150K').median(axis=1)
    ematrix.df['median_2K'] = ematrix.df.filter(like='2K').median(axis=1)
    ematrix.df['median_SN'] = ematrix.df.filter(like='SN').median(axis=1)

    # Step 2: Merge medians from ematrix.df into df based on 'miRNA'

    # Add median expressions to df based on 'miRNA' mapping
    df['median_expression'] = df['miRNA'].map(ematrix.df.set_index('name')['median_expression'])
    df['median_10K'] = df['miRNA'].map(ematrix.df.set_index('name')['median_10K'])
    df['median_150K'] = df['miRNA'].map(ematrix.df.set_index('name')['median_150K'])
    df['median_2K'] = df['miRNA'].map(ematrix.df.set_index('name')['median_2K'])
    df['median_SN'] = df['miRNA'].map(ematrix.df.set_index('name')['median_SN'])

    return df

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
    if not args.tmp:
        os.makedirs(args.output_folder + "/tmp", exist_ok=True)
        df = get_cooperative_targets(args.file, args.combinations)
        write_tsv(df, args.output_folder + "/tmp")
    else:
        df = pd.read_csv(args.tmp, sep="\t")
    df = add_expression_exosome(df, args.matrix)
    write_tsv(df, args.output_folder)

if __name__ == "__main__":
    main()
