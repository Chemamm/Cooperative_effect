import argparse
import pandas as pd
import os
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns
from upsetplot import UpSet
from upsetplot import from_contents
from itertools import combinations
from collections import Counter
from sklearn.linear_model import LinearRegression
import numpy as np
import scipy.stats as stats
from scipy.optimize import curve_fit
from scipy.stats import spearmanr



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
        for index, row in self.df.iterrows():
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
                    poscon_dict[row['miRNA']].append(
                        row['transcript'] + str(count_dict[row['miRNA'] + ";" + row['transcript']]))
                    count_dict[row['miRNA'] + ";" + row['transcript']] += 1
        return poscon_dict

    def get_gene_dict(self):
        gene_dict = {}
        for index, row in self.df.iterrows():
            if "|" in row['transcript']:
                row['transcript'] = row['transcript'].split("|")[0]

            code = row['miRNA'] + ';' + str(row['end'] - 6) + ':' + str(row['end'])
            if row['transcript'] not in gene_dict:
                gene_dict[row['transcript']] = []
            gene_dict[row['transcript']].append(code)
        return gene_dict




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

    def calculate_medians_groups(self):
        # Extract column groups
        groups = {
            '2K': [col for col in self.df.columns if '2K' in col],
            '10K': [col for col in self.df.columns if '10K' in col],
            '150K': [col for col in self.df.columns if '150K' in col],
            'SN': [col for col in self.df.columns if 'SN' in col]
        }

        # Calculate medians for each group
        for group_name, columns in groups.items():
            self.df[f'Median_{group_name}'] = self.df[columns].median(axis=1)

        # Calculate the median across all samples
        self.df['Median_All'] = self.df[self.df.columns[1:]].median(axis=1)


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

    def filter_by_number(self, n):
        filtered_df = self.df[self.df['Number of targets'] == n]
        return filtered_df





class InteractionExtractor:
    def __init__(self, filename):
        self.interactions = self.extract_interactions(filename)

    def extract_interactions(self, filename):
        interactions = {}
        with open(filename, 'r') as file:
            lines = list(file)
            for k in range(1,len(lines)):
                if lines[k].strip():  # Skip empty lines
                    parts = lines[k].split()
                    gene = parts[0].split('-')[0]  # Extract gene ID
                    if "|" in gene:
                        gene = gene.split("|")[0]
                    miRNAs = parts[4].split(';')    # Extract miRNAs
                if gene in interactions:
                    interactions[gene].append(miRNAs)
                else:
                    interactions[gene] = [miRNAs]
        for gene in interactions:
            interactions[gene] = interactions[gene]
        return interactions

class Fasta:
    def __init__(self, file):
        with open(file) as ft:
            self.lines = list(ft)
        self.dict = self.fasta2dict()

    def fasta2dict(self):
        fasta_dict = dict()
        seqid = ""
        for i in range(0, len(self.lines)):
            if self.lines[i].startswith('>'):
                if seqid != "":
                    fasta_dict[seqid] = seq
                seqid = self.lines[i].split(' ')[0].replace(">", "").replace("\n", "")
                if "|" in seqid:
                    seqid = seqid.split("|")[0]
                if "." in seqid:
                    seqid = seqid.split(".")[0]
                seq = ""
            else:
                seq += self.lines[i].replace("\n", "")
        if seqid != "":
            fasta_dict[seqid] = seq
        return fasta_dict

class ConservedTargets:
    def __init__(self, file_path):
        # Read the file into a DataFrame
        self.df = pd.read_csv(file_path, sep="\t", header=None)

    def count_miRNA_occurrences(self):
        # Count occurrences of each miRNA in the first column
        miRNA_counts = self.df[0].value_counts().to_dict()
        return miRNA_counts





def get_conserved_interactions(poscon_dict, mirna, threshold, species_dicts):
    conserved_interactions = []

    for gene in poscon_dict[mirna]:
        conservation_level = 1
        for species in species_dicts:
            if gene in species[mirna]:
                conservation_level += 1
        if conservation_level >= threshold:
            conserved_interactions.append(gene)
    return conserved_interactions



def get_cooperative_targets(file, combo, species_files, threshold):
    poscon = PositionalConsensus(file)
    poscon_dict = poscon.get_dict()
    combo = ValidCombination(combo)
    species_dicts = []
    for species in species_files.split(","):
        species = PositionalConsensus(species)
        species_dicts.append(species.get_dict())

    valid_2, valid_3 = combo.get_cooperative_targets()

    total_targets_count = {}
    conserved_targets_count = {}
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
            conserved_targets = get_conserved_interactions(poscon_dict, key, threshold, species_dicts)
            conserved_targets_count[key] = len(conserved_targets)
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
            if len(conserved_targets) == 0:
                ratio_2[key] = "NaN"
                ratio_3[key] = "NaN"
            else:
                ratio_2[key] = targets_valid_2 / len(conserved_targets)
                ratio_3[key] = targets_valid_3 / len(conserved_targets)
    # Convert dictionaries to DataFrames
    df_total_targets = pd.DataFrame(total_targets_count.items(), columns=['miRNA', 'total_targets_count'])
    df_conserved_targets = pd.DataFrame(conserved_targets_count.items(), columns=['miRNA', 'conserved_targets_count'])
    df_valid_2_targets = pd.DataFrame(targets_valid_2_count.items(), columns=['miRNA', 'targets_valid_2_count'])
    df_valid_3_targets = pd.DataFrame(targets_valid_3_count.items(), columns=['miRNA', 'targets_valid_3_count'])
    df_ratio_2 = pd.DataFrame(ratio_2.items(), columns=['miRNA', 'ratio_2'])
    df_ratio_3 = pd.DataFrame(ratio_3.items(), columns=['miRNA', 'ratio_3'])

    # Merge DataFrames on 'miRNA'
    df_merged = df_total_targets.merge(df_conserved_targets, on='miRNA', how='outer')
    df_merged = df_merged.merge(df_valid_2_targets, on='miRNA', how='outer')
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



def get_gene_count(file):
    gene_count_dict_2 = dict()
    gene_count_dict_3 = dict()
    gene_mirnas_dict = dict()
    with open(file) as fl:
        for line in fl.readlines()[1::]:
            fields = line.split("\t")
            gene = fields[0].split("-")[0]
            mirnas = fields[4]
            n_mirnas = fields[6]
            if n_mirnas == "2":
                if gene not in gene_count_dict_2:
                    gene_count_dict_2[gene] = 1
                    gene_mirnas_dict[gene] = mirnas
                else:
                    gene_count_dict_2[gene] += 1
                    gene_mirnas_dict[gene] =  gene_mirnas_dict[gene] + ";" + mirnas
            elif n_mirnas == "3":
                if gene not in gene_count_dict_3:
                    gene_count_dict_3[gene] = 1
                    gene_mirnas_dict[gene] = mirnas
                else:
                    gene_count_dict_3[gene] += 1
                    gene_mirnas_dict[gene] = gene_mirnas_dict[gene] + ";" + mirnas
    sorted_tuples = sorted(gene_count_dict_2.items(), key=lambda x: x[1], reverse=True)
    # Convert the sorted list of tuples back to a dictionary
    gene_count_dict_2 = dict(sorted_tuples)
    return gene_count_dict_2, gene_count_dict_3, gene_mirnas_dict

def write_tsv(gene_count_dict_2, gene_count_dict_3, output):
    print(" ".join(list(gene_count_dict_2.keys())[0:50]))
    output = output + "_genes.tsv"
    with open(output, "wt") as out:
        out.write("gene\tnumber of combinations of 2 miRNAs\tnumber of combinations of 3 miRNAs\n")
        for gene in gene_count_dict_2:
            n_3 = 0
            if gene in gene_count_dict_3:
                n_3 = gene_count_dict_3[gene]
            with open(output, "wt") as out:
                out.write("%s\t%i\t%i\n" %(gene, gene_count_dict_2[gene], n_3))

def write_top_mirnas_and_genes(gene_count_dict_2, gene_mirnas_dict, output, n):
    output_mirnas = output + "_top%i_mirnas.txt" %n
    output_genes = output + "_top%i_genes.txt" %n
    top_genes = list(gene_count_dict_2.keys())[0:n]
    with open(output_genes, "wt") as out:
        out.write("\n".join(top_genes))
    mirnas = []
    for gene in top_genes:
        mirnas = mirnas + gene_mirnas_dict[gene].split(";")
    mirnas = list(set(mirnas))
    with open(output_mirnas, "wt") as out:
        out.write("\n".join(mirnas))
    return mirnas


def expression_analysis(mirnas, ematrix, mirnas_library, output, n):
    output = output + "_top%i_boxplot.png" %n
    mirnas_library = Fasta(mirnas_library)
    all_mirnas = list(mirnas_library.dict.keys())
    expression_df = pd.read_csv(ematrix, sep="\t")
    other_mirnas = [x for x in all_mirnas if x not in mirnas]

    list1_df = expression_df[expression_df['name'].isin(mirnas)]
    list2_df = expression_df[expression_df['name'].isin(other_mirnas)]

    list1_mean = list1_df.mean(axis=1)
    list2_mean = list2_df.mean(axis=1)

    # Step 4: Combine mean expression values into a single DataFrame
    mean_df = pd.DataFrame({'miRNAs in conserved cooperations': list1_mean, 'Other miRNAs': list2_mean})

    # Step 5: Create a boxplot
    sns.boxplot(data=mean_df)
    plt.yscale('log')
    plt.ylabel('Log Mean Expression (RPM)')
    plt.savefig(output, dpi=300)

def get_interactions_other_species(species):
    species_interactions = []
    for file in species.split(","):
        inex = InteractionExtractor(file)
        species_interactions.append(inex.interactions)
    return species_interactions


def get_conservation_level(line, interactions_list):
    conservation_level = 1
    parts = line.split()
    gene = parts[0].split('-')[0]  # Extract gene ID
    if "|" in gene:
        gene= gene.split("|")[0]
    miRNAs = parts[4].split(';')


    for interaction in interactions_list:
        if gene in interaction:
            if isinstance(interaction[gene][0], list):
                for lst in interaction[gene]:
                    if set(miRNAs) == set(lst):
                        conservation_level +=1
                        break
            else:
                if set(miRNAs) == set(interaction[gene]):
                    conservation_level += 1

    return conservation_level

def get_conserved_combinations(filename, species, threshold, output_folder):
    species_interactions = get_interactions_other_species(species)
    # Create the output folder if it doesn't exist
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)

    tsv_path = os.path.join(output_folder, "conserved_combinations_%i.tsv" %threshold)

    with open(filename, 'r') as file:
        lines = list(file)
        with open(tsv_path, 'wt') as out:
            out.write(lines[0].strip("\n") + "\tconservation_level\n")
            for k in range(1, len(lines)):
                if lines[k].strip():  # Skip empty lines
                    clevel = get_conservation_level(lines[k], species_interactions)
                    if clevel >= threshold:
                        out.write(lines[k].strip("\n") + "\t%i\n" %clevel)



def calculate_valid_miRNA_combinations(cluster, start_var, end_var):
    valid_combinations_2 = []
    valid_combinations_3 = []
    regions = cluster['regions']
    miRNAs = [region['miRNA'] for region in cluster['regions']]

    # Generate all combinations of miRNAs within the cluster
    for r in range(2, 4):  # Generate combinations of length 2 to len(miRNAs)
        for combo in combinations(miRNAs, r):
            combo_regions = [regions[i] for i in range(len(regions)) if regions[i]['miRNA'] in combo]
            # Check if the combination of miRNAs is valid based on the specified rules
            if is_valid_miRNA_combination(combo_regions, combo, start_var, end_var):
                if r == 2:
                    valid_combinations_2.append({'transcript': cluster['transcript'], 'miRNAs': combo, 'regions': combo_regions})
                elif r == 3:
                    valid_combinations_3.append({'transcript': cluster['transcript'], 'miRNAs': combo, 'regions': combo_regions})
    return valid_combinations_2, valid_combinations_3


def is_valid_miRNA_combination(regions, miRNAs, start_var, end_var):
    # Check if the combination of miRNAs satisfies the specified rules
    is_valid = True
    for i in range(len(miRNAs) - 1):
        first_region = regions[i]
        second_region = regions[i + 1]

        # Check the rule for each pair of consecutive regions in the combination
        if (second_region['end'] - first_region['end']) > end_var or \
                (second_region['start'] - first_region['end']) < start_var:
            is_valid = False

    return is_valid

def add_info_to_combo(region_clusters):
    region_id_count = {}
    for cluster in region_clusters:
        targets = list(region['miRNA'] for region in cluster['regions'])
        num_targets = len(targets)
        cluster['Number of Targets'] = num_targets
        start = cluster['regions'][0]['start']
        end = cluster['regions'][-1]['end']
        region_id = f"{cluster['transcript']}-{start}-{end}"
        if region_id in region_id_count:
            region_id_count[region_id] += 1
            region_id = f"{cluster['transcript']}-{start}-{end}.{region_id_count[region_id]}"
        else:
            region_id_count[region_id] = 1

        cluster['Parent ID'] = region_id

    return region_clusters

def process_region_clusters(region_clusters, start_var, end_var):
    valid_combinations = []

    for cluster in region_clusters:
        # Calculate valid combinations of 2 and 3 miRNAs for the current cluster
        valid_combinations_2, valid_combinations_3 = calculate_valid_miRNA_combinations(cluster, start_var, end_var)

        # Count the number of valid combinations of 2 miRNAs
        count_2_miRNAs = sum(1 for combo in valid_combinations_2 if combo['transcript'] == cluster['transcript'])

        # Count the number of valid combinations of 3 miRNAs
        count_3_miRNAs = sum(1 for combo in valid_combinations_3 if combo['transcript'] == cluster['transcript'])

        # Update the cluster dictionary with counts
        cluster['Number of valid combinations of 2 miRNAs'] = count_2_miRNAs
        cluster['Number of valid combinations of 3 miRNAs'] = count_3_miRNAs

        # Extend the list of valid combinations
        valid_combinations.extend(valid_combinations_2)
        valid_combinations.extend(valid_combinations_3)

    valid_combinations = add_info_to_combo(valid_combinations)
    # Return the updated region_clusters (although this is typically not necessary)
    return region_clusters, valid_combinations


def generate_combination_id(transcript, combo):
    start = combo[0]['start']
    end = combo[-1]['end']
    return f"{transcript}-{start}-{end}"


def calculate_median_expression(ematrix, region_clusters):
    ematrix = Ematrix(ematrix)
    ematrix.add_mean_median_std()

    # Create a mapping of miRNA to its median expression
    miRNA_median_dict = ematrix.df.set_index('name')['median'].to_dict()

    # Calculate median expressions of miRNAs for each cluster
    cluster_median_expressions = []

    # Add 'Median expression' field to each region cluster
    for cluster in region_clusters:
        miRNAs = [region['miRNA'] for region in cluster['regions']]
        median_expressions = [miRNA_median_dict.get(miRNA, None) for miRNA in miRNAs if miRNA in miRNA_median_dict]
        if median_expressions:
            cluster_median_expression = pd.Series(median_expressions).median()
        else:
            cluster_median_expression = None

        cluster_median_expressions.append(cluster_median_expression)

        # Add 'Median expression' field to each region cluster
    for i, cluster in enumerate(region_clusters):
        cluster['Median expression'] = cluster_median_expressions[i]

    return region_clusters

def write_final_tsv(region_clusters, output_folder):
    # Create the output folder if it doesn't exist
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)

    # Write region clusters to a TSV file in the output folder
    tsv_path = os.path.join(output_folder, "final_output.tsv")
    with open(tsv_path, 'wt') as outfile:
        outfile.write("Region ID\tStart\tEnd\tmiRNAs\ttargets\tNumber of valid combinations of 2 miRNAs\tNumber of valid combinations of 3 miRNAs\tMedian expression\n")
        region_id_count = {}
        for idx, cluster in enumerate(region_clusters):
            start = cluster['regions'][0]['start']
            end = cluster['regions'][-1]['end']
            miRNAs = ';'.join(region['miRNA'] for region in cluster['regions'])
            targets = ';'.join(f"{region['start']}:{region['end']}" for region in cluster['regions'])
            region_id = f"{cluster['transcript']}-{start}-{end}"
            if region_id in region_id_count:
                region_id_count[region_id] += 1
                region_id = f"{cluster['transcript']}-{start}-{end}.{region_id_count[region_id]}"
            else:
                region_id_count[region_id] = 1

            # Write the cluster information to the TSV file
            outfile.write(f"{region_id}\t{start}\t{end}\t{miRNAs}\t{targets}\t"
                          f"{cluster['Number of valid combinations of 2 miRNAs']}\t"
                          f"{cluster['Number of valid combinations of 3 miRNAs']}\t{cluster['Median expression']}\n")


    print("Final output saved as TSV to", tsv_path)

def write_final_tsv_noexp(region_clusters, output_folder):
    # Create the output folder if it doesn't exist
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)

    # Write region clusters to a TSV file in the output folder
    tsv_path = os.path.join(output_folder, "final_output.tsv")
    with open(tsv_path, 'wt') as outfile:
        outfile.write("Region ID\tStart\tEnd\tmiRNAs\ttargets\tNumber of valid combinations of 2 miRNAs\tNumber of valid combinations of 3 miRNAs\n")
        region_id_count = {}
        for idx, cluster in enumerate(region_clusters):
            start = cluster['regions'][0]['start']
            end = cluster['regions'][-1]['end']
            miRNAs = ';'.join(region['miRNA'] for region in cluster['regions'])
            targets = ';'.join(f"{region['start']}:{region['end']}" for region in cluster['regions'])
            region_id = f"{cluster['transcript']}-{start}-{end}"
            if region_id in region_id_count:
                region_id_count[region_id] += 1
                region_id = f"{cluster['transcript']}-{start}-{end}.{region_id_count[region_id]}"
            else:
                region_id_count[region_id] = 1

            # Write the cluster information to the TSV file
            outfile.write(f"{region_id}\t{start}\t{end}\t{miRNAs}\t{targets}\t"
                          f"{cluster['Number of valid combinations of 2 miRNAs']}\t"
                          f"{cluster['Number of valid combinations of 3 miRNAs']}\n")

def write_combo_tsv(region_clusters, output_folder):
    # Create the output folder if it doesn't exist
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)

    # Write region clusters to a TSV file in the output folder
    tsv_path = os.path.join(output_folder, "valid_cooperations.tsv")
    with open(tsv_path, 'wt') as outfile:
        outfile.write("Cooperation ID\tStart\tEnd\tRegion ID\tmiRNAs\ttargets\tNumber of targets\tMedian expression\n")
        region_id_count = {}
        for idx, cluster in enumerate(region_clusters):
            start = cluster['regions'][0]['start']
            end = cluster['regions'][-1]['end']
            miRNAs = ';'.join(region['miRNA'] for region in cluster['regions'])
            targets = ';'.join(f"{region['start']}:{region['end']}" for region in cluster['regions'])
            region_id = f"{cluster['transcript']}-{start}-{end}"
            if region_id in region_id_count:
                region_id_count[region_id] += 1
                region_id = f"{cluster['transcript']}-{start}-{end}.{region_id_count[region_id]}"
            else:
                region_id_count[region_id] = 1

            # Write the cluster information to the TSV file
            outfile.write(f"{region_id}\t{start}\t{end}\t{cluster['Parent ID']}\t{miRNAs}\t{targets}\t"
                          f"{cluster['Number of Targets']}\t{cluster['Median expression']}\n")

    print("Final output saved as TSV to", tsv_path)

def write_combo_tsv_noexp(region_clusters, output_folder):
    # Create the output folder if it doesn't exist
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)

    # Write region clusters to a TSV file in the output folder
    tsv_path = os.path.join(output_folder, "valid_cooperations.tsv")
    with open(tsv_path, 'wt') as outfile:
        outfile.write("Cooperation ID\tStart\tEnd\tRegion ID\tmiRNAs\ttargets\tNumber of targets\n")
        region_id_count = {}
        for idx, cluster in enumerate(region_clusters):
            start = cluster['regions'][0]['start']
            end = cluster['regions'][-1]['end']
            miRNAs = ';'.join(region['miRNA'] for region in cluster['regions'])
            targets = ';'.join(f"{region['start']}:{region['end']}" for region in cluster['regions'])
            region_id = f"{cluster['transcript']}-{start}-{end}"
            if region_id in region_id_count:
                region_id_count[region_id] += 1
                region_id = f"{cluster['transcript']}-{start}-{end}.{region_id_count[region_id]}"
            else:
                region_id_count[region_id] = 1

            # Write the cluster information to the TSV file
            outfile.write(f"{region_id}\t{start}\t{end}\t{cluster['Parent ID']}\t{miRNAs}\t{targets}\t"
                          f"{cluster['Number of Targets']}\n")

    print("Final output saved as TSV to", tsv_path)
def plot_scatter_regression(df, output_folder):
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Define combinations of median expressions and cooperative ratios
    median_cols = ['median_expression', 'median_10K', 'median_150K', 'median_2K', 'median_SN']
    ratio_cols = ['ratio_2', 'ratio_3']

    # Iterate over each combination and create scatterplot with regression line
    for median_col in median_cols:
        for ratio_col in ratio_cols:
            # Filter out rows with NaN values in selected columns
            filtered_df = df[[median_col, ratio_col]].dropna()

            # Plot scatterplot with regression line
            sns.set(style="whitegrid")
            plt.figure(figsize=(8, 6))
            sns.regplot(x=ratio_col, y=median_col, data=filtered_df, scatter_kws={"s": 50})

            # Set Y-axis to logarithmic scale
            plt.yscale('log')

            # Set plot title and labels
            #plt.title(f"{median_col} vs {ratio_col}")
            plt.xlabel(f"{ratio_col}")
            plt.ylabel(f"{median_col}")

            # Save plot as PNG file
            plot_filename = f"{median_col}_vs_{ratio_col}.png"
            plot_path = os.path.join(output_folder, plot_filename)
            plt.savefig(plot_path)
            plt.close()


def plot_comparison_top_bottom(df, output_folder, n):
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Sort DataFrame by median_expression
    df_sorted = df.sort_values(by='median_expression', ascending=False)

    # Select top 20 and bottom 20 miRNAs based on median_expression
    top_miRNAs = df_sorted.head(n)['miRNA']
    bottom_miRNAs = df_sorted.tail(n)['miRNA']

    # Filter DataFrame for top and bottom miRNAs
    top_df = df[df['miRNA'].isin(top_miRNAs)]
    bottom_df = df[df['miRNA'].isin(bottom_miRNAs)]

    # Combine top and bottom DataFrames
    combined_df = pd.concat([top_df.assign(Group='Top expressed'), bottom_df.assign(Group='Low expressed')])

    # Create catplots for ratio_2 and ratio_3
    for ratio_col in ['ratio_2', 'ratio_3']:
        sns.set(style="whitegrid")
        plt.figure(figsize=(10, 8))
        sns.catplot(x='Group', y=ratio_col, data=combined_df, kind='box', height=6, aspect=1.5)
        # Set plot title and labels
        #plt.title(f"Comparison of {ratio_col} between Top and Low expressed miRNAs")
        plt.xlabel("Expression Group")
        plt.ylabel(f"{ratio_col}")

        # Save plot as PNG file
        plot_filename = f"comparison_{ratio_col}_top_bottom.png"
        plot_path = os.path.join(output_folder, plot_filename)
        plt.savefig(plot_path)
        plt.close()



def plot_scatter_regression_raw(df, output_folder):
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Define combinations of median expressions and cooperative raws
    median_cols = ['median_expression', 'median_10K', 'median_150K', 'median_2K', 'median_SN']
    raw_cols = ['targets_valid_2_count', 'targets_valid_3_count']

    # Iterate over each combination and create scatterplot with regression line
    for median_col in median_cols:
        for raw_col in raw_cols:
            # Filter out rows with NaN values in selected columns
            filtered_df = df[[median_col, raw_col]].dropna()

            # Plot scatterplot with regression line
            sns.set(style="whitegrid")
            plt.figure(figsize=(8, 6))
            sns.regplot(x=raw_col, y=median_col, data=filtered_df, scatter_kws={"s": 50})

            # Set Y-axis to logarithmic scale
            plt.yscale('log')

            # Set plot title and labels
            #plt.title(f"{median_col} vs {raw_col}")
            plt.xlabel(f"{raw_col}")
            plt.ylabel(f"{median_col}")

            # Save plot as PNG file
            plot_filename = f"{median_col}_vs_{raw_col}.png"
            plot_path = os.path.join(output_folder, plot_filename)
            plt.savefig(plot_path)
            plt.close()


def plot_comparison_top_bottom_raw(df, output_folder, n):
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Sort DataFrame by median_expression
    df_sorted = df.sort_values(by='median_expression', ascending=False)

    # Select top 20 and bottom 20 miRNAs based on median_expression
    top_miRNAs = df_sorted.head(n)['miRNA']
    bottom_miRNAs = df_sorted.tail(n)['miRNA']

    # Filter DataFrame for top and bottom miRNAs
    top_df = df[df['miRNA'].isin(top_miRNAs)]
    bottom_df = df[df['miRNA'].isin(bottom_miRNAs)]

    # Combine top and bottom DataFrames
    combined_df = pd.concat([top_df.assign(Group='Top expressed'), bottom_df.assign(Group='Low expressed')])

    # Create catplots for raw_2 and raw_3
    for raw_col in ['targets_valid_2_count', 'targets_valid_3_count']:
        sns.set(style="whitegrid")
        plt.figure(figsize=(10, 8))
        sns.catplot(x='Group', y=raw_col, data=combined_df, kind='box', height=6, aspect=1.5)

        # Set plot title and labels
        #plt.title(f"Comparison of {raw_col} between Top and Low expressed miRNAs")
        plt.xlabel("Expression Group")
        plt.ylabel(f"{raw_col}")

        # Save plot as PNG file
        plot_filename = f"comparison_{raw_col}_top_bottom.png"
        plot_path = os.path.join(output_folder, plot_filename)
        plt.savefig(plot_path)
        plt.close()



def upset_top_ratio_mirnas(file_list, species_list, output):
    file_list = file_list.split(",")
    species_list = species_list.split(",")
    # Read each dataframe into a dictionary with keys as species names
    dataframes = {
        species_list[0]: pd.read_csv(file_list[0], sep="\t"),
        species_list[1]: pd.read_csv(file_list[1], sep="\t"),
        species_list[2]: pd.read_csv(file_list[2], sep="\t"),
        species_list[3]: pd.read_csv(file_list[3], sep="\t"),
        species_list[4]: pd.read_csv(file_list[4], sep="\t"),
        species_list[5]: pd.read_csv(file_list[5], sep="\t")
    }

    # Initialize empty dictionaries to store the top 20 miRNAs for each criterion and species
    top_20_ratio_2 = {}
    top_20_targets_valid_2_count = {}

    # Loop through each species' dataframe
    for species, df in dataframes.items():
        # Sort the dataframe by 'ratio_2' and 'targets_valid_2_count' in descending order
        sorted_df_ratio = df.sort_values(by=['ratio_2'], ascending=False)
        sorted_df_raw = df.sort_values(by=['targets_valid_2_count'], ascending=False)

        # Extract the top 20 miRNAs based on 'ratio_2' and 'targets_valid_2_count'
        top_20_ratio_2[species] = sorted_df_ratio.head(20)['miRNA'].tolist()
        top_20_targets_valid_2_count[species] = sorted_df_raw.head(20)['miRNA'].tolist()


    upset_ratio = from_contents(top_20_ratio_2)
    # Create an UpSet plot for top 20 miRNAs by ratio_2
    UpSet(upset_ratio, subset_size="count").plot()
    plt.title("UpSet plot of Top 20 miRNAs by ratio_2")
    plt.savefig(output + "_ratio.png")
    plt.close()

    upset_raw = from_contents(top_20_targets_valid_2_count)
    # Create an UpSet plot for top 20 miRNAs by targets_valid_2_count
    UpSet(upset_raw, subset_size="count").plot()
    plt.title("UpSet plot of Top 20 miRNAs by targets_valid_2_count")
    plt.savefig(output + "_raw.png")
    plt.close()

    sets_ratio = [{value for value in sublist} for sublist in top_20_ratio_2.values()]
    sets_raw = [{value for value in sublist} for sublist in top_20_targets_valid_2_count.values()]

    output_intersections_ratio = output + "_intersections_ratio.txt"
    output_intersections_raw = output + "_intersections_raw.txt"

    with open(output_intersections_ratio, "wt") as out:
        for comboSize in range(1, len(sets_ratio) + 1):
            for combo in combinations(range(len(sets_ratio)), comboSize):
                intersection = sets_ratio[combo[0]]
                for i in combo[1:]:
                    intersection = intersection & sets_ratio[i]
                out.write(", ".join([species_list[i] for i in combo]) + " = " + str(intersection) + "\n")

    with open(output_intersections_raw, "wt") as out:
        for comboSize in range(1, len(sets_raw)+1):
            for combo in combinations(range(len(sets_raw)), comboSize):
                intersection = sets_raw[combo[0]]
                for i in combo[1:]:
                    intersection = intersection & sets_raw[i]
                out.write(", ".join([species_list[i] for i in combo]) + " = " + str(intersection) + "\n")

def sort_mirnas(miRNA_str):
    miRNA_list = miRNA_str.split(';')
    miRNA_list.sort()
    return ';'.join(miRNA_list)

def combinations_scatter_expression(combinations, ematrix, output):
    # Function to calculate the mean expression of the miRNAs in each combination
    def calculate_mean_expression(miRNA_combination, expression_df):
        miRNAs = miRNA_combination.split(';')
        expressions = []
        for miRNA in miRNAs:
            expression = expression_df.loc[expression_df['name'] == miRNA, 'median'].values

            if expression:
                expressions.append(expression[0])
            else:
                expressions.append(0)
        if expressions:
            return sum(expressions) / len(expressions)
        else:
            return None

    df = pd.read_csv(combinations, sep="\t")
    df = df[df["Number of targets"]==2]
    # Apply the function to the miRNAs column
    df['sorted_miRNAs'] = df['miRNAs'].apply(sort_mirnas)
    print(df)

    # Count the occurrences of each sorted combination
    miRNA_counts = Counter(df['sorted_miRNAs'])


    # Convert to DataFrame for better readability
    miRNA_counts_df = pd.DataFrame(miRNA_counts.items(), columns=['miRNA_combination', 'count'])
    if ematrix:
        ematrix = Ematrix(ematrix)
        ematrix.add_mean_median_std()


        # Calculate the mean expression for each miRNA combination
        miRNA_counts_df['mean_expression'] = miRNA_counts_df['miRNA_combination'].apply(
            lambda x: calculate_mean_expression(x, ematrix.df))

    miRNA_counts_df = miRNA_counts_df.sort_values(by='count', ascending=False)

    if ematrix:
        # Create a categorical plot
        sns.set(style="whitegrid")
        cat_plot = sns.catplot(x='count', y='mean_expression', data=miRNA_counts_df, kind='strip', height=6, aspect=1.5)

        # Set plot labels and title
        cat_plot.set_axis_labels('Count', 'Mean Expression')
        cat_plot.fig.suptitle('Mean Expression vs. Count of miRNA Combinations', y=1.02)
        plt.tight_layout()

        # Show the plot
        plt.savefig(output + "_count_vs_exp.png", dpi=300)
        plt.close()

    miRNA_list = miRNA_counts_df['miRNA_combination'].str.split(';').explode()
    miRNA_counts = miRNA_list.value_counts()

    # Select the top 15 miRNAs
    top_miRNAs = miRNA_counts.head(15)

    # Create a count plot
    sns.set(style="whitegrid")
    plt.figure(figsize=(9, 6))
    sns.barplot(x=top_miRNAs.values, y=top_miRNAs.index, color=sns.color_palette()[0])
    plt.xlabel('Number of Cooperative targets')
    plt.ylabel('miRNA')
    plt.title('Top 15 miRNAs by Number of Cooperative Targets')
    plt.tight_layout()
    plt.savefig(output + "_count_miRNAs.png", dpi=300)
    plt.close()


    # Select the top 15 most frequent combinations based on the 'count' column
    top_combinations = miRNA_counts_df.nlargest(15, 'count')

    # Create the second plot: Top 15 Most Frequent miRNA Combinations
    plt.figure(figsize=(9, 6))
    sns.barplot(x=top_combinations['count'], y=top_combinations['miRNA_combination'], color=sns.color_palette()[0])
    plt.xlabel('Frequency')
    plt.ylabel('miRNA Combination')
    plt.title('Top 15 Most Frequent miRNA Combinations')
    plt.tight_layout()
    plt.savefig(output + "_count_combinations.png", dpi=300)
    plt.close()

def filter_combinations_by_exp(combinations, ematrix, output):
    def check_miRNAs(row, expression_names):
        miRNAs = row['miRNAs'].split(';')
        return all(miRNA in expression_names for miRNA in miRNAs)

    df = pd.read_csv(combinations, sep="\t")

    ematrix = Ematrix(ematrix)

    expression_names = set(ematrix.df['name'])
    df = df[df.apply(check_miRNAs, axis=1, expression_names=expression_names)]

    df.to_csv(output, sep="\t", index=False)




def combinations_gene_analysis(combinations, output):
    df = pd.read_csv(combinations, sep="\t")
    df = df[df["Number of targets"] == 2]
    # Remove rows where 'Cooperation ID' contains a dot ('.')
    df = df[~df['Cooperation ID'].str.contains('\.')]

    # Create a new column 'gene'
    df['gene'] = df['Cooperation ID'].str.split('-').str[0]

    # Generate a count plot of the top 20 genes

    top_genes = df['gene'].value_counts().head(20)
    sns.set(style="whitegrid")
    plt.figure(figsize=(9, 6))
    sns.barplot(x=top_genes.values, y=top_genes.index, color=sns.color_palette()[0])
    plt.xlabel('Number of Cooperative targets')
    plt.ylabel('Gene')
    plt.title('Top 20 Genes by Number of Cooperative Targets')
    plt.tight_layout()
    plt.savefig(output + "_top_20_genes.png", dpi=300)

    # Save the DataFrame to CSV files
    genes = df['gene'].value_counts().index
    with open(output + '_all_genes.txt', "wt") as out:
        out.write("\n".join(genes))

    # Save the top 500 frequent genes
    top_500_genes = df['gene'].value_counts().head(500).index
    with open(output + '_top_500_genes.txt', "wt") as out:
        out.write("\n".join(top_500_genes))

    # Save the top 100 frequent genes
    top_100_genes = df['gene'].value_counts().head(100).index
    with open(output + '_top_100_genes.txt', "wt") as out:
        out.write("\n".join(top_100_genes))

    # Filter genes with count >= 3
    filtered_genes = df['gene'].value_counts()[df['gene'].value_counts() >= 3].index
    with open(output + '_genes_3combinations.txt', "wt") as out:
        out.write("\n".join(filtered_genes))



def expression_ctargets_stats(conserved_targets, iri_100RPM, ematrix, output):
    iri_100RPM = Fasta(iri_100RPM)
    iri_100RPM = iri_100RPM.dict.keys()
    conserved_targets = ConservedTargets(conserved_targets)
    conserved_targets_count = conserved_targets.count_miRNA_occurrences()

    # Classify miRNAs based on abundance
    saliva_abundant_counts = []
    non_saliva_abundant_counts = []
    miRNA_classification = []

    for miRNA, count in conserved_targets_count.items():
        if miRNA in iri_100RPM:
            saliva_abundant_counts.append(count)
            miRNA_classification.append({'miRNA': miRNA, 'Count': count, 'Class': 'Saliva abundant'})
        else:
            non_saliva_abundant_counts.append(count)
            miRNA_classification.append({'miRNA': miRNA, 'Count': count, 'Class': 'Non-saliva abundant'})

    # Prepare data for plotting
    plot_data = pd.DataFrame(miRNA_classification)

    class_order = ['Saliva abundant', 'Non-saliva abundant']

    # Create boxplot
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='Class', y='Count', data=plot_data, order=class_order)
    plt.title('Distribution of the number of conserved targets by miRNA')
    plt.savefig(output + "_ctargets_boxplot.png", dpi=300)

    # Create catplot
    sns.catplot(x='Class', y='Count', kind='swarm', data=plot_data, height=6, aspect=1.5, order=class_order)
    plt.title('Distribution of the number of conserved targets by miRNA')
    plt.savefig(output + "_ctargets_catplot.png", dpi=300)

    # Create violin plot
    plt.figure(figsize=(10, 6))
    sns.violinplot(x='Class', y='Count', data=plot_data, order=class_order)
    plt.title('Distribution of the number of conserved targets by miRNA')
    plt.savefig(output + "_ctargets_violinplot.png", dpi=300)

    # Scatter with expression

    ematrix = Ematrix(ematrix)
    ematrix.calculate_medians_groups()

    all_miRNAs = set(ematrix.df['name']).union(set(conserved_targets_count.keys()))
    expression_df = ematrix.df.set_index('name').reindex(all_miRNAs, fill_value=0).reset_index()

    # Add miRNA occurrence counts to the DataFrame
    expression_df['Occurrences'] = expression_df['name'].map(conserved_targets_count).fillna(0)

    # Calculate log of median expressions
    for group in ['2K', '10K', '150K', 'SN', 'All']:
        expression_df[f'Log_Median_{group}'] = np.log1p(expression_df[f'Median_{group}'])

    # Plot scatter plots
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('Scatter Plots with Regression Lines')

    plot_scatter_with_regression(expression_df['Occurrences'], expression_df['Log_Median_2K'], '2K Group', axes[0, 0])
    plot_scatter_with_regression(expression_df['Occurrences'], expression_df['Log_Median_10K'], '10K Group', axes[0, 1])
    plot_scatter_with_regression(expression_df['Occurrences'], expression_df['Log_Median_150K'], '150K Group',
                                 axes[1, 0])
    plot_scatter_with_regression(expression_df['Occurrences'], expression_df['Log_Median_SN'], 'SN Group', axes[1, 1])

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output + "_ctargets_scatterbygroup.png", dpi=300)
    plt.close()

    # Plot scatter plot for Median_All
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_scatter_with_regression(expression_df['Occurrences'], expression_df['Log_Median_All'], 'All Samples', ax)
    plt.savefig(output + "_ctargets_scatterAll.png", dpi=300)
    plt.close()

    # Plot scatter plots no zeros
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('Scatter Plots with Regression Lines')

    plot_scatter_with_nonlinear_regression_no_zeros(expression_df['Occurrences'], expression_df['Log_Median_2K'],
                                                    '2K Group', axes[0, 0])
    plot_scatter_with_nonlinear_regression_no_zeros(expression_df['Occurrences'], expression_df['Log_Median_10K'],
                                                    '10K Group', axes[0, 1])
    plot_scatter_with_nonlinear_regression_no_zeros(expression_df['Occurrences'], expression_df['Log_Median_150K'],
                                                    '150K Group',
                                                    axes[1, 0])
    plot_scatter_with_nonlinear_regression_no_zeros(expression_df['Occurrences'], expression_df['Log_Median_SN'],
                                                    'SN Group', axes[1, 1])

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output + "_ctargets_scatterbygroup_nozeros.png", dpi=300)
    plt.close()

    # Plot scatter plot for Median_All
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_scatter_with_nonlinear_regression_no_zeros(expression_df['Occurrences'], expression_df['Log_Median_All'],
                                                    'All Samples', ax)
    plt.savefig(output + "_ctargets_scatterAll_nozeros.png", dpi=300)
    plt.close()


def expression_cooperations_stats(conserved_cooperations, iri_100RPM, ematrix, output):
    iri_100RPM = Fasta(iri_100RPM)
    iri_100RPM = iri_100RPM.dict.keys()
    conserved_coop = ValidCombination(conserved_cooperations)
    conserved_coop_2, conserved_coop_3 = conserved_coop.get_cooperative_targets()
    # Dictionary to store miRNA counts
    miRNA_counts = {}

    # Iterate over the list and count occurrences
    for item in conserved_coop_2:
        miRNA_name = item.split(';')[1]
        if miRNA_name in miRNA_counts:
            miRNA_counts[miRNA_name] += 1
        else:
            miRNA_counts[miRNA_name] = 1

    # Classify miRNAs based on abundance
    saliva_abundant_counts = []
    non_saliva_abundant_counts = []
    miRNA_classification = []

    for miRNA, count in miRNA_counts.items():
        if miRNA in iri_100RPM:
            saliva_abundant_counts.append(count)
            miRNA_classification.append({'miRNA': miRNA, 'Count': count, 'Class': 'Saliva abundant'})
        else:
            non_saliva_abundant_counts.append(count)
            miRNA_classification.append({'miRNA': miRNA, 'Count': count, 'Class': 'Non-saliva abundant'})

    # Prepare data for plotting
    plot_data = pd.DataFrame(miRNA_classification)

    class_order = ['Saliva abundant', 'Non-saliva abundant']

    # Create boxplot
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='Class', y='Count', data=plot_data, order=class_order)
    plt.title('Distribution of the number of conserved targets by miRNA')
    plt.savefig(output + "_coop_boxplot.png", dpi=300)

    # Create catplot
    sns.catplot(x='Class', y='Count', kind='swarm', data=plot_data, height=6, aspect=1.5, order=class_order)
    plt.title('Distribution of the number of conserved targets by miRNA')
    plt.savefig(output + "_coop_catplot.png", dpi=300)

    # Create violin plot
    plt.figure(figsize=(10, 6))
    sns.violinplot(x='Class', y='Count', data=plot_data, order=class_order)
    plt.title('Distribution of the number of conserved targets by miRNA')
    plt.savefig(output + "_coop_violinplot.png", dpi=300)

    # Scatter with expression

    ematrix = Ematrix(ematrix)
    ematrix.calculate_medians_groups()

    all_miRNAs = set(ematrix.df['name']).union(set(miRNA_counts.keys()))
    expression_df = ematrix.df.set_index('name').reindex(all_miRNAs, fill_value=0).reset_index()

    # Add miRNA occurrence counts to the DataFrame
    expression_df['Occurrences'] = expression_df['name'].map(miRNA_counts).fillna(0)

    # Calculate log of median expressions
    for group in ['2K', '10K', '150K', 'SN', 'All']:
        expression_df[f'Log_Median_{group}'] = np.log1p(expression_df[f'Median_{group}'])

    # Plot scatter plots
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('Scatter Plots with Regression Lines')

    plot_scatter_with_regression(expression_df['Occurrences'], expression_df['Log_Median_2K'], '2K Group', axes[0, 0])
    plot_scatter_with_regression(expression_df['Occurrences'], expression_df['Log_Median_10K'], '10K Group', axes[0, 1])
    plot_scatter_with_regression(expression_df['Occurrences'], expression_df['Log_Median_150K'], '150K Group',
                                 axes[1, 0])
    plot_scatter_with_regression(expression_df['Occurrences'], expression_df['Log_Median_SN'], 'SN Group', axes[1, 1])

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output + "_coop_scatterbygroup.png", dpi=300)
    plt.close()

    # Plot scatter plot for Median_All
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_scatter_with_regression(expression_df['Occurrences'], expression_df['Log_Median_All'], 'All Samples', ax)
    plt.savefig(output + "_coop_scatterAll.png", dpi=300)
    plt.close()


    # Plot scatter plots no zeros
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('Scatter Plots with Regression Lines')

    plot_scatter_with_nonlinear_regression_no_zeros(expression_df['Occurrences'], expression_df['Log_Median_2K'], '2K Group', axes[0, 0])
    plot_scatter_with_nonlinear_regression_no_zeros(expression_df['Occurrences'], expression_df['Log_Median_10K'], '10K Group', axes[0, 1])
    plot_scatter_with_nonlinear_regression_no_zeros(expression_df['Occurrences'], expression_df['Log_Median_150K'], '150K Group',
                                 axes[1, 0])
    plot_scatter_with_nonlinear_regression_no_zeros(expression_df['Occurrences'], expression_df['Log_Median_SN'], 'SN Group', axes[1, 1])

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output + "_coop_scatterbygroup_nozeros.png", dpi=300)
    plt.close()

    # Plot scatter plot for Median_All
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_scatter_with_nonlinear_regression_no_zeros(expression_df['Occurrences'], expression_df['Log_Median_All'], 'All Samples', ax)
    plt.savefig(output + "_coop_scatterAll_nozeros.png", dpi=300)
    plt.close()


def plot_scatter_with_regression(x, y, title, ax):
    ax.clear()
    sns.scatterplot(x=x, y=y, ax=ax)
    ax.set_title(title)


    # Fit non-linear regression
    # Fit non-linear regression
    popt, _ = curve_fit(nonlinear_func, x, y)
    x_fit = np.linspace(min(x), max(x), 100)
    y_fit = nonlinear_func(x_fit, *popt)

    # Plot regression line
    ax.plot(x_fit, y_fit, color='red', linewidth=2)
    # Calculate R-squared
    residuals = y - nonlinear_func(x, *popt)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)

    # Display R-squared and parameters
    params = f'a = {popt[0]:.2f}, b = {popt[1]:.2f}'
    r_squared_info = f'R^2 = {r_squared:.2f}\nParameters: {params}'
    ax.annotate(r_squared_info, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=12,
                ha='left', va='top', bbox=dict(boxstyle='round,pad=0.5', fc='white', ec='black', lw=1))


def plot_scatter_with_nonlinear_regression_no_zeros(x, y, title, ax):
    # Convert to numpy arrays if they are not
    x = np.array(x)
    y = np.array(y)

    # Filter out values with 0 in x
    valid_indices = y > 0
    x_valid = x[valid_indices]
    y_valid = y[valid_indices]

    # Plot with non-linear regression
    plot_scatter_with_regression(x_valid, y_valid, title, ax)


def nonlinear_func(x, a, b):
    return a * np.log(x) + b

def hsa_plus_iri_cooperations_analysis(coop, output):
    df = pd.read_csv(coop, sep="\t")

    # Step 1: Filter the DataFrame for Number of targets == 2
    filtered_df = df[df['Number of targets'] == 2]

    # Step 2: Standardize Cooperation IDs
    filtered_df['Standardized Cooperation ID'] = filtered_df['Cooperation ID'].apply(lambda x: x.split('.')[0])

    # Step 3: Perform the analysis for each Cooperation ID
    result_data = []
    cooperation_ids = filtered_df['Standardized Cooperation ID'].unique()

    for coop_id in cooperation_ids:
        coop_data = filtered_df[filtered_df['Standardized Cooperation ID'] == coop_id]
        iri_flag = False
        hsa_flag = False
        mixed_flag = False
        mimick_flag = False
        for _, row in coop_data.iterrows():
            miRNAs = row['miRNAs'].split(';')
            if all(mirna.lower().startswith('iri') for mirna in miRNAs):
                iri_flag = True
            elif all(mirna.lower().startswith('hsa') for mirna in miRNAs):
                hsa_flag = True
            elif any(mirna.lower().startswith('iri') for mirna in miRNAs) and any(
                    mirna.lower().startswith('hsa') for mirna in miRNAs):
                mixed_flag = True

        if hsa_flag and (mixed_flag or iri_flag):
            mimick_flag = True

        result_data.append({
            'Cooperation ID': coop_id,
            'Iri': iri_flag,
            'Hsa': hsa_flag,
            'Mixed': mixed_flag,
            'Mimick': mimick_flag
        })

    # Create the result DataFrame
    result_df = pd.DataFrame(result_data)
    print(result_df)

    # Count the occurrences
    count_iri = result_df['Iri'].sum()
    count_hsa = result_df['Hsa'].sum()
    count_mixed = result_df['Mixed'].sum()
    count_mimick = result_df['Mimick'].sum()

    # Print the counts
    print(f"Number of Hsa: {count_hsa}")
    print(f"Number of Iri: {count_iri}")
    print(f"Number of Mixed: {count_mixed}")
    print(f"Number of Mimick: {count_mimick}")

    # Create a barplot
    counts = {
        'Hsa': count_hsa,
        'Iri': count_iri,
        'Mixed': count_mixed,
        'Mimick': count_mimick
    }

    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.barplot(x=list(counts.keys()), y=list(counts.values()))
    plt.title('Counts of Hsa, Iri, Mixed, and Mimick')
    plt.xlabel('Category')
    plt.ylabel('Count')
    plt.savefig(output + "_barplot_stats.png", dpi=300)

    # Step 3: Remove rows with two miRNAs starting with "Hsa" or "hsa"
    def has_only_hsa_mirnas(miRNAs):
        mirna_list = miRNAs.split(';')
        return all(mirna.lower().startswith('hsa') for mirna in mirna_list)

    filtered_df = filtered_df[~filtered_df['miRNAs'].apply(has_only_hsa_mirnas)]

    # Step 4: Fix the .X annotation
    def fix_annotation(group):
        ids = group['Cooperation ID'].tolist()
        base_ids = group['Standardized Cooperation ID'].unique()
        for idx, base_id in enumerate(base_ids):
            # Reset the counter for each base ID
            count = 1
            for i in range(len(ids)):
                if ids[i].startswith(base_id):
                    if count == 1:
                        group.loc[group['Cooperation ID'] == ids[i], 'Cooperation ID'] = base_id
                    else:
                        group.loc[group['Cooperation ID'] == ids[i], 'Cooperation ID'] = f"{base_id}.{count}"
                    count += 1
        return group

    fixed_df = filtered_df.groupby('Standardized Cooperation ID').apply(fix_annotation).reset_index(drop=True)
    fixed_df.to_csv(output + "_iri_cooperations.tsv", sep="\t", index=False)


def hsa_plus_iri_cooperations_analysis_iri_activated(coop, output):

    coop = ValidCombination(coop)
    # Initialize an empty list to store hsa miRNAs and their corresponding targets
    hsa_list = []

    # Initialize a dictionary to store unique target regions for each gene
    unique_target_regions = {}

    filtered_df = coop.df[coop.df['Number of targets'] == 2]


    # First loop to populate the hsa_list
    for index, row in filtered_df.iterrows():
        miRNAs = row['miRNAs'].split(';')
        targets = row['targets'].split(';')

        # Check if both miRNAs start with 'Hsa' or 'hsa'
        if all(miRNA.lower().startswith('hsa') for miRNA in miRNAs):
            for miRNA, target in zip(miRNAs, targets):
                hsa_list.append(f"{miRNA}-{target}")

    # Initialize counters for the categories
    iri_only_count = 0
    hsa_activated_count = 0
    iri_activated_count = 0

    # Second loop to categorize the cooperations
    gene_counter = Counter()
    gene_counter_iri_activated = Counter()

    for index, row in filtered_df.iterrows():
        miRNAs = row['miRNAs'].split(';')
        targets = row['targets'].split(';')
        gene = row['Region ID'].split('-')[0]  # Extract the gene ID

        # Create a unique identifier for the target regions
        target_combination = tuple(sorted(targets))

        # Check if this combination of target regions has already been counted for this gene
        if gene not in unique_target_regions:
            unique_target_regions[gene] = set()



        # Check if both miRNAs start with 'iri'
        if all(miRNA.lower().startswith('iri') for miRNA in miRNAs):
            if target_combination not in unique_target_regions[gene]:
                unique_target_regions[gene].add(target_combination)
                iri_only_count += 1
                gene_counter[gene] += 1
                gene_counter_iri_activated[gene] += 1

        # Check if one is 'iri' and one is 'hsa'
        elif any(miRNA.lower().startswith('iri') for miRNA in miRNAs) and any(
                miRNA.lower().startswith('hsa') for miRNA in miRNAs):
            # Determine which is 'hsa' and which is 'iri'
            for miRNA, target in zip(miRNAs, targets):
                if miRNA.lower().startswith('hsa'):
                    hsa_mirna_target = f"{miRNA}-{target}"
                    if hsa_mirna_target in hsa_list:
                        if target_combination not in unique_target_regions[gene]:
                            unique_target_regions[gene].add(target_combination)
                            hsa_activated_count += 1
                    else:
                        if target_combination not in unique_target_regions[gene]:
                            unique_target_regions[gene].add(target_combination)
                            iri_activated_count += 1
                            gene_counter_iri_activated[gene] += 1
                    break  # Only need to check one 'hsa' miRNA
            gene_counter[gene] += 1



    # Print the results
    print(f"IRI Only: {iri_only_count}")
    print(f"HSA Activated: {hsa_activated_count}")
    print(f"IRI Activated: {iri_activated_count}")

    # Create DataFrames for the output files
    df_numcoop = pd.DataFrame(gene_counter.items(), columns=['Gene', 'NumCoop'])
    df_numcoop_iri_activated = pd.DataFrame(gene_counter_iri_activated.items(), columns=['Gene', 'NumCoopIriActivated'])

    # Rank the genes based on NumCoop and NumCoopIriActivated using the 'average' method for ties
    df_numcoop['Rank'] = df_numcoop['NumCoop'].rank(method='average', ascending=False)
    df_numcoop_iri_activated['Rank'] = df_numcoop_iri_activated['NumCoopIriActivated'].rank(method='average',
                                                                                            ascending=False)

    # Write the DataFrames to CSV files
    df_numcoop.to_csv(output + '_counts.tsv', sep="\t", index=False)
    df_numcoop_iri_activated.to_csv(output + '_iri_activated_counts.tsv', sep="\t", index=False)



    # Barplot for the cooperation categories
    category_counts = {
        'IRI Only': iri_only_count,
        'HSA Activated': hsa_activated_count,
        'IRI Activated': iri_activated_count
    }

    sns.set_theme()
    plt.figure(figsize=(8, 6))
    plt.bar(category_counts.keys(), category_counts.values())
    plt.title('Cooperation Categories')
    plt.xlabel('Category')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig(output + "_iri_activated_countplot.png", dpi=300)

    # Barplot for the top 20 genes with the highest number of iri_only + iri_activated cooperations
    top_genes = gene_counter.most_common(20)
    genes, counts = zip(*top_genes)
    with open(output + "_iri_activated_top20_genes.txt", "wt") as out:
        out.write("\n".join(genes))

    plt.figure(figsize=(12, 8))
    plt.barh(genes, counts, color='purple')
    plt.title('Top 20 Genes by IRI Only + IRI Activated Cooperations')
    plt.xlabel('Count')
    plt.ylabel('Gene')
    plt.gca().invert_yaxis()  # Invert y-axis to have the gene with the highest count at the top
    plt.tight_layout()
    plt.savefig(output + "_iri_activated_top20_genes.png", dpi=300)

def remove_redundant_target_regions(gene_dict):
    filtered_gene_dict = {}

    for gene, targets in gene_dict.items():
        coord_dict = {}
        for target in targets:
            miRNA, coord = target.split(';')
            if coord not in coord_dict:
                coord_dict[coord] = target

        filtered_gene_dict[gene] = list(coord_dict.values())

    return filtered_gene_dict



def significance_numcoop(hsa_targets, coop, coop_random, output_prefix, min_positives=3):

    hsa_tar = PositionalConsensus(hsa_targets).get_gene_dict()
    hsa_tar_non_redundant = remove_redundant_target_regions(hsa_tar)
    coop = ValidCombination(coop)
    coop_2, coop_3 = coop.get_cooperative_targets()
    count_dict = {}

    for gene in hsa_tar_non_redundant:
        count_dict[gene] = [0,0]

        for mirna_code in hsa_tar_non_redundant[gene]:
            code = gene + ";" + mirna_code
            if code in coop_2:
                count_dict[gene][0] += 1
            else:
                count_dict[gene][1] += 1



    # Filter the dictionary based on the number of positive cases
    filtered_gene_data = {gene: (pos, neg) for gene, (pos, neg) in count_dict.items() if pos >= min_positives}


    proportion_dict = {}
    for gene in filtered_gene_data:
        proportion_dict[gene] = filtered_gene_data[gene][0]/(filtered_gene_data[gene][0] + filtered_gene_data[gene][1])



    coop_random = ValidCombination(coop_random)
    coop_random_2, coop_random_3 = coop_random.get_cooperative_targets()
    count_random_dict = {}

    for gene in hsa_tar_non_redundant:
        count_random_dict[gene] = [0, 0]

        for mirna_code in hsa_tar_non_redundant[gene]:
            code = gene + ";" + mirna_code
            if code in coop_random_2:
                count_random_dict[gene][0] += 1
            else:
                count_random_dict[gene][1] += 1


    results = []

    for gene in filtered_gene_data:
        positive, negative = filtered_gene_data[gene]
        positive_control, negative_control = count_random_dict[gene]

        # Contingency table for chi-square test
        contingency_table = [[positive, negative], [positive_control, negative_control]]
        # Perform chi-square test
        #chi2, p_value, _, _ = stats.chi2_contingency(contingency_table)
        # Perform Fisher's exact test
        oddsratio, p_value = stats.fisher_exact(contingency_table)

        # Store the results
        results.append({'Gene': gene, 'Positive': positive, 'Negative': negative, 'Positive control':positive_control,
                        'Negative control': negative_control, 'OddsRatio': oddsratio, 'P-Value': p_value})

    # Convert results to a DataFrame for easy viewing
    results_df = pd.DataFrame(results)

    # Sort the results by p-value
    results_df_sorted = results_df.sort_values(by='P-Value')

    p_values = results_df_sorted['P-Value'].values
    m = len(p_values)

    if m > 0:
        sorted_indices = np.argsort(p_values)
        sorted_p_values = p_values[sorted_indices]
        fdr_bh = np.zeros(m)
        for i in range(m):
            fdr_bh[i] = sorted_p_values[i] * m / (i + 1)
        fdr_bh = np.minimum.accumulate(fdr_bh[::-1])[::-1]  # Adjust to ensure monotonicity
        results_df_sorted['FDR_BH'] = np.take(fdr_bh, np.argsort(sorted_indices))

    # Calculate FDR using Benjamini-Yekutieli procedure
    harmonic_sum = np.sum(1 / np.arange(1, m + 1))
    results_df_sorted['FDR_BY'] = np.nan
    if m > 0:
        fdr_by = np.zeros(m)
        for i in range(m):
            fdr_by[i] = sorted_p_values[i] * m / (i + 1) * (1 / harmonic_sum)
        fdr_by = np.minimum.accumulate(fdr_by[::-1])[::-1]
        results_df_sorted['FDR_BY'] = np.take(fdr_by, np.argsort(sorted_indices))

    # Calculate q-values using Storey’s method
    results_df_sorted['Q-Value'] = np.nan
    if m > 0:
        pi0 = 1 - (np.sum(p_values <= 0.05) / m)  # Estimate pi0
        q_values = np.zeros(m)
        for i in range(m):
            q_values[i] = p_values[i] * m / (i + 1) / pi0
        q_values = np.minimum.accumulate(q_values[::-1])[::-1]
        results_df_sorted['Q-Value'] = np.take(q_values, np.argsort(sorted_indices))


    results_df_sorted.to_csv(output_prefix + "_data.tsv", sep="\t", index=False)


    # Plotting the p-value distribution
    plt.figure(figsize=(10, 6))
    sns.histplot(results_df['P-Value'], bins=10, kde=True)
    plt.title('P-Value Distribution')
    plt.xlabel('P-Value')
    plt.ylabel('Frequency')
    plt.savefig(f'{output_prefix}_p_value_distribution.png', dpi=300)
    plt.close()

    # Plotting the top 20 genes with the most significant p-values
    top_20_genes = results_df_sorted.head(20)

    plt.figure(figsize=(12, 8))
    sns.barplot(x='P-Value', y='Gene', data=top_20_genes, palette='viridis')
    plt.title('Top 20 Genes with Most Significant P-Values')
    plt.xlabel('P-Value')
    plt.ylabel('Gene')
    plt.savefig(f'{output_prefix}_top_20_genes.png', dpi=300)
    plt.close()

    # Scatter plot of positive cases vs. p-values
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x='Positive', y='P-Value', data=results_df)
    plt.title('Positive Cases vs. P-Value')
    plt.xlabel('Number of Positive Cases')
    plt.ylabel('P-Value')
    plt.savefig(f'{output_prefix}_positive_cases_vs_p_value.png', dpi=300)
    plt.close()


def spearman_rank_analysis(problem_file, control_file, output):
    # Read the CSV files into DataFrames
    df_problem = pd.read_csv(problem_file, sep="\t")
    df_control = pd.read_csv(control_file, sep="\t")

    # Ensure the Gene column is used as the index for easier joining
    df_problem.set_index('Gene', inplace=True)
    df_control.set_index('Gene', inplace=True)

    # Join the DataFrames on the index (Gene)
    df_joined = df_problem.join(df_control, lsuffix='_problem', rsuffix='_control')

    # Drop any rows with missing values to ensure a proper comparison
    df_joined.dropna(inplace=True)

    # Extract the ranks
    ranks_problem = df_joined['Rank_problem']
    ranks_control = df_joined['Rank_control']

    # Calculate Spearman's rank correlation coefficient
    corr_coefficient, p_value = spearmanr(ranks_problem, ranks_control)
    print(f"Spearman Rank Correlation Coefficient: {corr_coefficient}")
    print(f"P-Value: {p_value}")

    # Scatter plot of the ranks
    sns.set_theme
    plt.figure(figsize=(10, 6))
    plt.scatter(ranks_problem, ranks_control, alpha=0.5)
    plt.title('Spearman Rank Correlation Scatter Plot')
    plt.xlabel('Problem Set Rank')
    plt.ylabel('Control Set Rank')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output + "_scatter.png", dpi=300)

    # Histogram of rank differences
    rank_diff = ranks_problem - ranks_control
    plt.figure(figsize=(10, 6))
    plt.hist(rank_diff, bins=20, alpha=0.7, color='orange')
    plt.title('Histogram of Rank Differences (Problem - Control)')
    plt.xlabel('Rank Difference')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output + "_historank_differences.png", dpi=300)

    # Density plot of rank differences
    plt.figure(figsize=(10, 6))
    rank_diff.plot(kind='density', color='green')
    plt.title('Density Plot of Rank Differences (Problem - Control)')
    plt.xlabel('Rank Difference')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output + "_densirank_difference.png", dpi=300)


