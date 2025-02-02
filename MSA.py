
import pandas as pd
import os
import json
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, DBSCAN
import numpy as np
import seaborn as sns
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import random
import joypy
from scipy.stats import fisher_exact

class PositionalConsensus:
    def __init__(self, file):
        self.df = pd.read_csv(file, sep="\t")
        with open(file) as fl:
            self.lines = list(fl)

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

    def add_genomic_coordinates(self, structure_file, chrAlias, utr_only, output):
        structure = TranscriptStructure(structure_file, utr_only)
        if chrAlias:
            chrAlias = ChromAlias(chrAlias)
        with open(output, "wt") as out:
            out.write("miRNA\ttranscript\tstart\tend\tstructure\tenergy\tframe\ttarget_seq\tmiRNA_seq\t"
                      "chr_ensembl\tchr_ucsc\tstrand\tgenomic_start\tgenomic_end\n")
            for line in self.lines[1::]:
                fields = line.strip("\n").split("\t")
                mirna = fields[0]
                gene = fields[1]
                if "." in gene:
                    gene = gene.split(".")[0]
                if "|" in gene:
                    gene = gene.split("|")[0]
                start = fields[2]
                end = fields[3]
                genomic_start, genomic_end, chr, strand = structure.get_genomic_coordinate_regions(gene, start, end)
                chr_ucsc = chrAlias.dict[chr]
                if strand == -1:
                    strand = "-"
                elif strand == 1:
                    strand = "+"
                out.write(line.strip("\n") + "\t" + chr + "\t" + chr_ucsc + "\t" +  strand + "\t" + genomic_start
                          + "\t" + genomic_end + "\n")



    def get_full_target_region(self, miranda, ts, output):
        ts_dict = {}
        miranda_dict = {}
        with open(ts) as ts:
            for line in ts.readlines()[1::]:
                fields = line.split("\t")
                mirna = fields[0]
                gene = fields[1]
                start = fields[3]
                end = fields[4]
                identifier = "%s:%s:%i" %(mirna,gene, int(end) - 1)
                ts_dict[identifier] = (start,end)

        with open(miranda) as miranda:
            for line in miranda.readlines()[1::]:
                fields = line.split("\t")
                mirna = fields[0]
                gene = fields[1]
                start = fields[3]
                end = fields[4]
                identifier = "%s:%s:%i" %(mirna,gene, int(end) - 1)
                miranda_dict[identifier] = (start,end)

        with open(output, "wt") as out:
            out.write(self.lines[0])
            for line in self.lines[1::]:
                fields = line.split("\t")
                mirna = fields[0]
                gene = fields[1]
                end = fields[3]
                identifier = "%s:%s:%i" %(mirna,gene, int(end))
                if identifier in ts_dict:
                    out.write("%s\t%s\t%s\t%s\n" %(mirna, gene, ts_dict[identifier][0], ts_dict[identifier][1]))
                elif identifier in miranda_dict:
                    out.write("%s\t%s\t%s\t%s\n" % (mirna, gene, miranda_dict[identifier][0], miranda_dict[identifier][1]))
                else:
                    out.write(line)
                    print("%s not in miranda or ts" %identifier)

class ChromAlias:
    def __init__(self,file):
        self.filename = file
        self.dict = self.get_dict()

    def get_dict(self):
        chrom_dict = {}
        with open(self.filename) as fl:
            for line in fl.readlines()[1::]:
                fields = line.strip("\n").split("\t")
                chrom_dict[fields[0]]=fields[1]
        return chrom_dict

class TranscriptStructure:
    def __init__(self, file, utr_only=False):
        self.df = pd.read_csv(file, sep="\t")
        self.gene_exon_dict = self.get_gene_exon_dict(utr_only)

    def get_gene_exon_dict(self, utr_only):
        # Group by "Gene stable ID" and apply the function
        if utr_only:
            gene_exon_dict = self.df.groupby("Gene stable ID").apply(generate_utr_coordinates).to_dict()
        else:
            gene_exon_dict = self.df.groupby("Gene stable ID").apply(generate_exon_coordinates).to_dict()

        return gene_exon_dict


    def get_genomic_coordinate_regions(self, gene, start_position, end_position):

        ## Getting strand
        strand = int(self.df[self.df["Gene stable ID"] == gene].iloc[0]['Strand'])

        chr =  str(self.df[self.df["Gene stable ID"] == gene].iloc[0]['Chromosome/scaffold name'])

        genomic_start, genomic_end = translate_position_to_genomic(self.gene_exon_dict[gene], int(start_position), int(end_position), strand)


        return genomic_start, genomic_end, chr, strand

class GenomicStructure:
    def __init__(self,file):
        self.df = pd.read_csv(file, sep="\t")

    def parse_strand_negtive_coords(self):
        def process_row(row):
            if row['strand'] == '-':
                # Swap genomic_start and genomic_end
                row['genomic_start'], row['genomic_end'] = row['genomic_end'], row['genomic_start']

                # If the values contain ",", split by ",", reverse the order, and join by ","
                if ',' in row['genomic_start']:
                    row['genomic_start'] = ','.join(reversed(row['genomic_start'].split(',')))
                if ',' in row['genomic_end']:
                    row['genomic_end'] = ','.join(reversed(row['genomic_end'].split(',')))
            return row

        df_new = self.df.apply(process_row, axis=1)
        df_new = df_new.sort_values(by='chr_ucsc')
        return df_new

    def order_by_chr(self):
        # Create a new column with the first value from 'genomic_start' split by commas if present
        self.df['genomic_start_temp'] = self.df['genomic_start'].apply(
            lambda x: int(str(x).split(',')[0]) if ',' in str(x) else int(x))

        # Sort the dataframe by 'chr_ucsc' and the new 'genomic_start_temp' column
        self.df = self.df.sort_values(by=['chr_ucsc', 'genomic_start_temp'])

        # Drop the temporary column
        self.df = self.df.drop(columns=['genomic_start_temp'])

class SpeciesTaxon:
    def __init__(self,file):
        with open(file) as fl:
            self.lines = list(fl)
        self.taxon_dict = self.get_taxon_dict()

    def get_taxon_dict(self):
        taxon_dict = {}
        for line in self.lines:
            fields = line.strip("\n").split("\t")
            taxon = fields[1]
            name = fields[2].replace(" ", "_")
            taxon_dict[name] = taxon

        return taxon_dict

    def get_taxon_ids(self,names_file):
        taxon_ids = []
        with open(names_file) as fl:
            for line in fl.readlines():
                name = line.strip("\n").replace(" ","_")
                taxon_ids.append(self.taxon_dict[name])

        print(",".join(taxon_ids))

class JSONmiRNA:
    def __init__(self,file):
        with open(file, 'r') as file:
            self.data = json.load(file)

    def preprocess_vectors(self):
        # Initialize lists to store the binary and conservation vectors
        binary_vectors = []
        conservation_vectors = []
        conservation_dict = {}
        max_length = max([max([int(x) for x in self.data[miRNA][k]['bindings']]) for miRNA in self.data for k in range(0, len(self.data[miRNA]))])

        for miRNA in self.data:
            for k in range(0, len(self.data[miRNA])):
                bindings = self.data[miRNA][k]['bindings']
                bind_cons = self.data[miRNA][k]['bind_cons']
                gene = self.data[miRNA][k]['line'].split("\t")[2].split("|")[0]
                line = self.data[miRNA][k]['line']


                # Transform bindings to binary vector
                binary_vector = [int(bindings[str(i)]) if str(i) in bindings else 0 for i in range(max_length)]

                # Initialize conservation vector with 0 or NaN
                conservation_vector = [0] * max_length
                conservation_species = [0] * max_length

                # Fill conservation vector based on bind_cons data
                for position in bind_cons:
                    score = float(bind_cons[position][1]) / float(bind_cons[position][2])
                    conservation_species[int(position) - 1] = bind_cons[position][4]
                    conservation_vector[int(position) - 1] = score


                # Append the vectors to the lists
                binary_vectors.append(binary_vector)
                conservation_vectors.append(conservation_vector)

                if miRNA in conservation_dict:
                    conservation_dict[miRNA].append((gene, conservation_vector, conservation_species, line))
                else:
                    conservation_dict[miRNA] = [(gene, conservation_vector, conservation_species, line)]


        # Convert lists to numpy arrays
        binary_vectors = np.array(binary_vectors)
        conservation_vectors = np.array(conservation_vectors)


        return conservation_vectors, binary_vectors, conservation_dict

    def unsupervised_analysis(self, n_cl, output):
        conservation_vectors, binary_vectors, conservation_dict = self.preprocess_vectors()

        feature_matrix = np.hstack((binary_vectors, conservation_vectors))

        pca = PCA(n_components=2)
        reduced_features = pca.fit_transform(feature_matrix)

        # Clustering using K-means
        kmeans = KMeans(n_clusters=n_cl)
        clusters = kmeans.fit_predict(reduced_features)

        # Visualization
        plt.scatter(reduced_features[:, 0], reduced_features[:, 1], c=clusters, s=0.1)
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        plt.title('miRNA-Target Clustering Kmeans')
        plt.colorbar(label='Cluster')
        plt.savefig(output + "_Kmeans.png", dpi=300)
        plt.close()

        dbscan = DBSCAN(eps=0.5, min_samples=500)
        clusters = dbscan.fit_predict(reduced_features)

        # Visualization
        plt.scatter(reduced_features[:, 0], reduced_features[:, 1], c=clusters, s=0.1)
        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        plt.title('miRNA-Target Clustering DBSCAN')
        plt.colorbar(label='Cluster')
        plt.savefig(output + "_DBSCAN.png", dpi=300)
        plt.close()

        out_boxplot = output + "_conservation_boxplot.png"
        print(conservation_vectors)

        plot_conservation_boxplot(conservation_vectors, out_boxplot)

        out_boxplot = output + "_conservation_violinplot.png"

        plot_conservation_violinplot(conservation_vectors, out_boxplot)

        out_boxplot = output + "_conservation_heatmap.png"

        plot_conservation_lineplot(conservation_vectors, out_boxplot)

        out_boxplot = output + "_conservation_lineplot.png"

        plot_conservation_heatmap(conservation_vectors, out_boxplot)

        out_boxplot = output + "_conservation_swarmplot.png"

        plot_conservation_swarmplot(conservation_vectors, out_boxplot)

        out_boxplot = output + "_conservation_ridgeplot.png"

        plot_conservation_ridgeline(conservation_vectors, out_boxplot)

        out_barplot = output + "_pairing_barplot.png"
        print(binary_vectors)

        plot_pairing_barplot(binary_vectors, out_barplot)

    def miRNAs_analysis(self, taxon_file, output):
        conservation_vectors, binary_vectors, conservation_dict = self.preprocess_vectors()

        filter_and_plot_conservation(conservation_dict, taxon_file, output)

class miRNAannot:
    def __init__(self, file):
        # Create a new list to store the modified rows
        cleaned_rows = []

        # Open the file and process it line by line
        with open(file, 'r') as file:
            for line in file:
                cleaned_rows.append(",".join(line.split(",")[0:6]))  # Keep only the part without the last field

        cleaned_rows = list(set(cleaned_rows))
        cleaned_data = "\n".join(cleaned_rows)
        with open("tmp_annot.txt", "wt") as out:
            out.write(cleaned_data)
        # Define the column names (since the file doesn't have a header)
        column_names = ['Seed', 'miRNA', 'Arm', 'Type', 'Proposed name', 'Taxon level']

        # Read the TSV file into a DataFrame, specifying no header and providing column names
        self.df = pd.read_csv("tmp_annot.txt", header=None, names=column_names)
        os.system("rm tmp_annot.txt")

    def write_df(self, output):
        self.df.to_csv(output, sep="\t", index=False)



def plot_conservation_boxplot(conservation_vectors, output):
    max_length = conservation_vectors.shape[1]
    conservation_data = []

    for i in range(max_length):
        position_data = conservation_vectors[:, i]
        # Filter out zeros since they represent unpaired positions
        position_data = position_data[position_data != 0]
        conservation_data.append(position_data)

    plt.figure(figsize=(15, 6))
    plt.boxplot(conservation_data)
    plt.xlabel('Position')
    plt.ylabel('Conservation Ratio')
    plt.title('Distribution of Conservation Ratios per Position')
    plt.savefig(output, dpi=300)
    plt.close()

def plot_conservation_violinplot(conservation_vectors, output):
    max_length = conservation_vectors.shape[1]
    conservation_data = []
    positions = []

    for i in range(max_length):
        position_data = conservation_vectors[:, i]
        position_data = position_data[position_data != 0]
        conservation_data.extend(position_data)
        positions.extend([i + 1] * len(position_data))  # 1-based index for positions

    # Create a DataFrame for easier plotting with seaborn
    df = pd.DataFrame({'Position': positions, 'Conservation Ratio': conservation_data})

    plt.figure(figsize=(15, 6))
    sns.violinplot(x='Position', y='Conservation Ratio', data=df, inner='box')
    plt.xlabel('Position')
    plt.ylabel('Conservation Ratio')
    plt.title('Distribution of Conservation Ratios per Position')
    plt.savefig(output, dpi=300)
    plt.close()

def plot_conservation_heatmap(conservation_vectors, output):
    max_length = conservation_vectors.shape[1]
    mean_conservation = []

    for i in range(max_length):
        position_data = conservation_vectors[:, i]
        position_data = position_data[position_data != 0]
        mean_conservation.append(position_data.mean())

    # Create a DataFrame for the heatmap
    df = pd.DataFrame({'Position': range(1, max_length + 1), 'Mean Conservation Ratio': mean_conservation})
    df = df.set_index('Position')

    plt.figure(figsize=(15, 6))
    sns.heatmap(df.T, cmap='viridis', annot=True, cbar=True)
    plt.xlabel('Position')
    plt.ylabel('')
    plt.title('Heatmap of Mean Conservation Ratios per Position')
    plt.savefig(output, dpi=300)
    plt.close()

def plot_conservation_lineplot(conservation_vectors, output):
    max_length = conservation_vectors.shape[1]
    mean_conservation = []
    std_conservation = []

    for i in range(max_length):
        position_data = conservation_vectors[:, i]
        position_data = position_data[position_data != 0]
        mean_conservation.append(position_data.mean())
        std_conservation.append(position_data.std())

    positions = range(1, max_length + 1)

    plt.figure(figsize=(15, 6))
    plt.errorbar(positions, mean_conservation, yerr=std_conservation, fmt='-o')
    plt.xlabel('Position')
    plt.ylabel('Mean Conservation Ratio')
    plt.title('Mean Conservation Ratios per Position with Error Bars')
    plt.savefig(output, dpi=300)
    plt.close()

def plot_conservation_swarmplot(conservation_vectors, output):
    max_length = conservation_vectors.shape[1]
    conservation_data = []
    positions = []

    for i in range(max_length):
        position_data = conservation_vectors[:, i]
        position_data = position_data[position_data != 0]
        conservation_data.extend(position_data)
        positions.extend([i + 1] * len(position_data))

    # Create a DataFrame for easier plotting with seaborn
    df = pd.DataFrame({'Position': positions, 'Conservation Ratio': conservation_data})

    plt.figure(figsize=(15, 6))
    sns.swarmplot(x='Position', y='Conservation Ratio', data=df)
    plt.xlabel('Position')
    plt.ylabel('Conservation Ratio')
    plt.title('Swarm Plot of Conservation Ratios per Position')
    plt.savefig(output, dpi=300)
    plt.close()

def plot_conservation_ridgeline(conservation_vectors, output):
    max_length = conservation_vectors.shape[1]
    conservation_data = []
    positions = []

    for i in range(max_length):
        position_data = conservation_vectors[:, i]
        position_data = position_data[position_data != 0]
        conservation_data.extend(position_data)
        positions.extend([i + 1] * len(position_data))

    # Create a DataFrame for easier plotting with joypy
    df = pd.DataFrame({'Position': positions, 'Conservation Ratio': conservation_data})

    plt.figure(figsize=(15, 8))
    joypy.joyplot(data=df, by="Position", column="Conservation Ratio", ylim='own', figsize=(15, 8))
    plt.xlabel('Position')
    plt.ylabel('Conservation Ratio')
    plt.title('Ridgeline Plot of Conservation Ratios per Position')
    plt.savefig(output, dpi=300)
    plt.close()

def plot_pairing_barplot(binary_vectors, output):
    max_length = binary_vectors.shape[1]
    pairing_counts = np.sum(binary_vectors, axis=0)

    plt.figure(figsize=(15, 6))
    plt.bar(range(max_length), pairing_counts)
    plt.xlabel('Position')
    plt.ylabel('Number of Pairings')
    plt.title('Number of Pairings per Position')
    plt.savefig(output, dpi=300)
    plt.close()


def filter_and_plot_conservation(conservation_dict, taxon_file, output, theshold=0.9):
    filtered_dict = {}
    primate_list = ["Primates","Haplorrhini","Simiiformes","Catarrhini","Hominoidea","Hominidae","Homininae","Homo"]


    for miRNA, gene_cons_tuples in conservation_dict.items():
        filtered_tuples = [t for t in gene_cons_tuples if all(pos >= theshold for pos in t[1][1:8])]
        filtered_tuples = [t for t in filtered_tuples if all(pos not in primate_list for pos in t[2][1:8])]
        filtered_dict[miRNA] = filtered_tuples

    miRNA_fractions = {miRNA: len(filtered_dict[miRNA]) / len(conservation_dict[miRNA]) for miRNA in conservation_dict
                       if len(conservation_dict[miRNA]) > 0}
    miRNA_counts = {miRNA: len(filtered_dict[miRNA]) for miRNA in conservation_dict}

    sorted_by_fraction = sorted(miRNA_fractions.items(), key=lambda x: x[1], reverse=True)
    sorted_by_count = sorted(miRNA_counts.items(), key=lambda x: x[1], reverse=True)

    top_fraction_miRNAs, top_fractions = zip(*sorted_by_fraction[:20])
    top_count_miRNAs, top_counts = zip(*sorted_by_count[:20])

    fraction_miRNAs, fractions = zip(*sorted_by_fraction)
    count_miRNAs, counts = zip(*sorted_by_count)

    out_barplot = output + "_miRNAs_cons.png"
    plt.figure(figsize=(15, 6))

    # Plot for the highest fractions
    plt.subplot(1, 2, 1)
    plt.bar(top_fraction_miRNAs, top_fractions)
    plt.xlabel('miRNA')
    plt.ylabel('Fraction of Conserved targets')
    plt.xticks(rotation=45, ha='right')

    fractions_mean = 0.017513251365241837
    fractions_mean_1sd = 0.030566339171218607
    fractions_mean_2sd = 0.04361942697719537

    plt.axhline(y=fractions_mean, color='red', linestyle='-', linewidth=2, label='Random Mean')
    plt.axhline(y=fractions_mean_1sd, color='orange', linestyle='--', label='Random Mean + 1 SD')
    plt.axhline(y=fractions_mean_2sd, color='orange', linestyle=':', label='Random Mean + 2 SD')

    plt.legend()

    # Plot for the highest raw counts
    plt.subplot(1, 2, 2)
    plt.bar(top_count_miRNAs, top_counts)
    plt.xlabel('miRNA')
    plt.ylabel('Number of Conserved targets')
    plt.xticks(rotation=45, ha='right')

    counts_mean = 45.485915492957744
    counts_mean_1sd = 95.5113647790416
    counts_mean_2sd = 145.53681406512547

    plt.axhline(y=counts_mean, color='red', linestyle='-', linewidth=2, label='Random Mean')
    plt.axhline(y=counts_mean_1sd, color='orange', linestyle='--', label='Random Mean + 1 SD')
    plt.axhline(y=counts_mean_2sd, color='orange', linestyle=':', label='Random Mean + 2 SD')

    plt.legend()

    plt.tight_layout()
    plt.savefig(out_barplot, dpi=300)
    plt.close()

    # Calculate and print statistics for fractions
    fractions_mean = np.mean(fractions)
    fractions_std = np.std(fractions)
    print(f"Fractions Mean: {fractions_mean}")
    print(f"Fractions Mean + 1 SD: {fractions_mean + fractions_std}")
    print(f"Fractions Mean + 2 SD: {fractions_mean + 2 * fractions_std}")

    # Calculate and print statistics for counts
    counts_mean = np.mean(counts)
    counts_std = np.std(counts)
    print(f"Counts Mean: {counts_mean}")
    print(f"Counts Mean + 1 SD: {counts_mean + counts_std}")
    print(f"Counts Mean + 2 SD: {counts_mean + 2 * counts_std}")

    if taxon_file != False:
        miRNA_taxon_classified = classify_miRNAs_by_taxon(taxon_file,
                                                          conservation_dict)
        create_catplots(miRNA_taxon_classified,conservation_dict, filtered_dict, output)

    out_tsv = output + "_filtered.tsv"

    with open(out_tsv, "wt") as out:
        for mirna in filtered_dict:
            for tuple in filtered_dict[mirna]:
                out.write(tuple[3] + "\n")



def classify_miRNAs_by_taxon(taxon_file, conservation_dict):
    dataframe = pd.read_csv(taxon_file, sep="\t")
    miRNA_taxon = dataframe[['miRNA', 'Taxon level']].drop_duplicates().set_index('miRNA').to_dict()['Taxon level']
    miRNA_taxon_classified = {}

    for miRNA in conservation_dict:
        if miRNA in miRNA_taxon:
            taxon_level = miRNA_taxon[miRNA]
            if taxon_level in miRNA_taxon_classified:
                miRNA_taxon_classified[taxon_level].append(miRNA)
            else:
                miRNA_taxon_classified[taxon_level] = [miRNA]

    return miRNA_taxon_classified


def create_catplots(miRNA_taxon_classified, conservation_dict, filtered_dict, output):
    miRNA_counts = {miRNA: len(filtered_dict[miRNA]) for miRNA in filtered_dict}
    miRNA_fractions = {miRNA: len(filtered_dict[miRNA]) / len(conservation_dict[miRNA]) for miRNA in filtered_dict}

    count_data = []
    fraction_data = []

    for taxon, miRNAs in miRNA_taxon_classified.items():
        if taxon not in ["Bilateria", "Protostomia", "Arthropoda"]:
            taxon = "I. ricinus"
        for miRNA in miRNAs:
            count = miRNA_counts.get(miRNA, 0)
            fraction = miRNA_fractions.get(miRNA, 0)
            count_data.append((taxon, miRNA, count))
            fraction_data.append((taxon, miRNA, fraction))

    count_df = pd.DataFrame(count_data, columns=['Taxon Level', 'miRNA', 'Count'])
    fraction_df = pd.DataFrame(fraction_data, columns=['Taxon Level', 'miRNA', 'Fraction'])
    count_output = output + "_taxon_catplot_count.png"
    fraction_output = output + "_taxon_catplot_fraction.png"
    # Define the taxon order
    taxon_order = ['Bilateria', 'Protostomia', 'Arthropoda', 'I. ricinus']

    plt.figure(figsize=(15, 6))
    sns.catplot(data=count_df, x='Taxon Level', y='Count', kind='strip', height=6, aspect=2, order=taxon_order)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(count_output, dpi=300)
    plt.close()

    plt.figure(figsize=(15, 6))
    sns.catplot(data=fraction_df, x='Taxon Level', y='Fraction', kind='strip', height=6, aspect=2, order=taxon_order)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(fraction_output, dpi=300)
    plt.close()

def plot_gene_stats(filtered_homolog, taxon_file, output):
    taxon_dict = {}
    with open(taxon_file, 'r') as tx_file:
        next(tx_file)
        for line in tx_file:
            seed, mirna, arm, name, taxon = line.strip().split('\t')
            if taxon not in ["Bilateria", "Protostomia", "Arthropoda"]:
                taxon = "I. ricinus"
            taxon_dict[mirna] = taxon



    gene_targets = {}

    # Open and parse the TSV file
    with open(filtered_homolog, 'r') as file:
        next(file)  # Skip header
        for line in file:
            miRNA, gene, start, end = line.strip().split('\t')[0:4]
            gene = gene.split("|")[0]

            # Count total target regions for each gene
            if gene not in gene_targets:
                gene_targets[gene] = {'total': 0, 'Bilateria': 0, "Protostomia": 0, "Arthropoda": 0, "I. ricinus": 0}

            gene_targets[gene]['total'] += 1
            if miRNA in taxon_dict:
                gene_targets[gene][taxon_dict[miRNA]] += 1
            else:
                gene_targets[gene]["I. ricinus"] += 1

        # Convert gene_targets into a DataFrame for easy manipulation
        df = pd.DataFrame.from_dict(gene_targets, orient='index')

        # Sort genes by the total number of target regions
        df_sorted = df.sort_values(by='total', ascending=False)

        # Write to TSV
        df_sorted.to_csv(output + '_gene_stats.tsv', sep='\t', index=True)

        # Select the top 20 and top 50 genes
        top_20 = df_sorted.head(20)
        top_50 = df_sorted.head(50)

        # Plot stacked bar plot for top 20 genes
        plot_stacked_bar(top_20, '', output + '_top20_stacked_conserved.png')

        # Plot stacked bar plot for top 50 genes
        plot_stacked_bar(top_50, '', output + '_top50_stacked_conserved_.png')

        # Plot total bar plot for top 20 genes
        plot_total_bar(top_20, '', output + '_top20_total_conserved.png')

        # Plot total bar plot for top 50 genes
        plot_total_bar(top_50, '', output + '_top50_total_conserved.png')

def plot_genes_stats_odds_ratio(poscon_real, poscon_conserved_real, poscon_random, poscon_conserved_random, output):
    def apply_smoothing_and_calculate_odds(conserved_real, non_conserved_real, conserved_random, non_conserved_random):
        # Additive smoothing
        smoothed_contingency_table = [
            [conserved_real + 1, non_conserved_real + 1],
            [conserved_random + 1, non_conserved_random + 1]
        ]
        # Calculate the odds ratio using the smoothed values
        odds_ratio, p_value = fisher_exact(smoothed_contingency_table, alternative='greater')
        return odds_ratio, p_value

    # Load real data
    real_df = pd.read_csv(poscon_real, sep='\t')
    real_df['gene'] = real_df['transcript'].apply(lambda x: x.split('|')[0])  # Extract gene ID

    # Load conservation data (real conserved targets)
    conservation_df = pd.read_csv(poscon_conserved_real, sep='\t', header=None, names=[
        'miRNA', 'transcript', 'start', 'end', 'structure', 'energy', 'conservation',
        'target_seq', 'conserved_seq', 'unknown_col', 'chromosome', 'strand', 'genomic_start', 'genomic_end'
    ])
    conservation_df['gene'] = conservation_df['transcript'].apply(lambda x: x.split('|')[0])  # Extract gene ID

    # Load randomized data
    randomized_df = pd.read_csv(poscon_random, sep='\t')
    randomized_df['gene'] = randomized_df['transcript'].apply(lambda x: x.split('|')[0])  # Extract gene ID

    # Load conservation random data (randomized conserved targets)
    conservation_random_df = pd.read_csv(poscon_conserved_random, sep='\t', header=None, names=[
        'miRNA', 'transcript', 'start', 'end', 'structure', 'energy', 'conservation',
        'target_seq', 'conserved_seq', 'unknown_col', 'chromosome', 'strand', 'genomic_start', 'genomic_end'
    ])
    conservation_random_df['gene'] = conservation_random_df['transcript'].apply(
        lambda x: x.split('|')[0])  # Extract gene ID

    # Group by gene (transcript) and count unique miRNA interactions for real and random data
    all_targets_per_gene = real_df.groupby('gene')['miRNA'].count()
    conserved_targets_per_gene = conservation_df.groupby('gene')['miRNA'].count()
    all_random_targets_per_gene = randomized_df.groupby('gene')['miRNA'].count()
    conserved_random_targets_per_gene = conservation_random_df.groupby('gene')['miRNA'].count()


    # Prepare to calculate odds ratios for each gene
    odds_ratios = []
    p_values = []
    total_real_counts = []
    real_conserved_counts = []
    total_random_counts = []
    random_conserved_counts = []

    for gene in all_targets_per_gene.index:
        # Real and randomized gene interaction counts
        all_real = all_targets_per_gene.get(gene, 0)
        all_random = all_random_targets_per_gene.get(gene, 0)
        conserved_real = conserved_targets_per_gene.get(gene, 0)
        conserved_random = conserved_random_targets_per_gene.get(gene, 0)

        # Non-conserved targets
        non_conserved_real = all_real - conserved_real
        non_conserved_random = all_random - conserved_random

        # Calculate odds ratio and p-value
        odds_ratio, p_value = apply_smoothing_and_calculate_odds(conserved_real, non_conserved_real, conserved_random,
                                                                 non_conserved_random)

        # Store results
        odds_ratios.append(odds_ratio)
        p_values.append(p_value)
        total_real_counts.append(all_real)
        real_conserved_counts.append(conserved_real)
        total_random_counts.append(all_random)
        random_conserved_counts.append(conserved_random)

    # Create a DataFrame to store the results
    gene_stats_df = pd.DataFrame({
        'Gene': all_targets_per_gene.index,
        'Real Total Targets': total_real_counts,
        'Real Conserved Targets': real_conserved_counts,
        'Random Total Targets': total_random_counts,
        'Random Conserved Targets': random_conserved_counts,
        'Odds Ratio': odds_ratios,
        'P-value': p_values
    })
    gene_stats_df['-log10(P-value)'] = -np.log10(gene_stats_df['P-value'])

    # Save the results to a CSV file
    gene_stats_df.to_csv(output + "_gene_target_statistics.csv", index=False)

    # Filter out rows where the p-value is greater than 0.05
    filtered_gene_stats_df = gene_stats_df[gene_stats_df['P-value'] <= 0.05]

    # Function to plot the top N genes by a specific column
    def plot_top_genes(df, column, title, ylabel, filename, top_n=20):
        top_genes = df.sort_values(by=column, ascending=False).head(top_n)
        plt.figure(figsize=(10, 6))
        plt.bar(top_genes['Gene'], top_genes[column])
        plt.title(title)
        plt.xlabel('Gene')
        plt.ylabel(ylabel)
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig(output + filename, dpi=300)
        plt.show()

    # Plot the top 20 genes by total targets
    plot_top_genes(gene_stats_df, 'Real Total Targets', 'Top 20 Genes by Total Targets', 'Total Targets',
                   '_top20_total_targets.png')

    # Plot the top 50 genes by total targets
    plot_top_genes(gene_stats_df, 'Real Total Targets', 'Top 50 Genes by Total Targets', 'Total Targets',
                   '_top50_total_targets.png', top_n=50)

    # Plot the top 20 genes by odds ratio
    plot_top_genes(filtered_gene_stats_df, 'Odds Ratio', 'Top 20 Genes by Odds Ratio', 'Odds Ratio', '_top20_odds_ratios.png')

    # Plot the top 50 genes by odds ratio
    plot_top_genes(filtered_gene_stats_df, 'Odds Ratio', 'Top 50 Genes by Odds Ratio', 'Odds Ratio', '_top50_odds_ratios.png',
                   top_n=50)

    # Plot the top 20 genes by -log10(p-value)
    plot_top_genes(gene_stats_df, '-log10(P-value)', 'Top 20 Genes by -log10(P-value)', '-log10(P-value)',
                   '_top20_logp_values.png')

    # Plot the top 50 genes by -log10(p-value)
    plot_top_genes(gene_stats_df, '-log10(P-value)', 'Top 50 Genes by -log10(P-value)', '-log10(P-value)',
                   '_top50_logp_values.png', top_n=50)


def plot_stacked_bar(df, title, output_file):
    ax = df[['Bilateria', 'Protostomia', 'Arthropoda', 'I. ricinus']].plot(kind='bar', stacked=True,
                                                                           figsize=(10, 7))
    plt.title(title)
    plt.xlabel('Genes')
    plt.ylabel('Number of Target Regions')
    plt.xticks(rotation=90)
    plt.legend(title='Taxon')
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()

def plot_total_bar(df, title, output_file):
    ax = df['total'].plot(kind='bar', color='skyblue', figsize=(10, 7))
    plt.title(title)
    plt.xlabel('Genes')
    plt.ylabel('Total Number of Target Regions')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()




def plot_miRNA_statistics(conservation_dict, filtered_cons):
    mirnas = list(conservation_dict.keys())
    total_counts = [len(conservation_dict[mirna]) for mirna in mirnas]
    filtered_counts = [len(filtered_cons.get(mirna, [])) for mirna in mirnas]
    fractions = [filtered / total if total > 0 else 0 for filtered, total in zip(filtered_counts, total_counts)]

    x = np.arange(len(mirnas))

    plt.figure(figsize=(15, 6))

    plt.subplot(1, 2, 1)
    plt.bar(x, fractions)
    plt.xticks(x, genes, rotation='vertical')
    plt.xlabel('miRNAs')
    plt.ylabel('Fraction of Conserved targets')

    plt.subplot(1, 2, 2)
    plt.bar(x, filtered_counts)
    plt.xticks(x, genes, rotation='vertical')
    plt.xlabel('miRNAs')
    plt.ylabel('Number of Conserved targets')

    plt.tight_layout()
    plt.show()


def generate_exon_coordinates(group):
    group = group.sort_values(by="Exon region start (bp)")

    coordinates = [f"{start}:{end}" for start, end in
                   zip(group["Exon region start (bp)"], group["Exon region end (bp)"])]
    return ";".join(coordinates)

def generate_utr_coordinates(group):
    # Drop rows with NaN values
    group = group.dropna(subset=["3' UTR start", "3' UTR end"])

    group = group.sort_values(by="3' UTR start")

    coordinates = [f"{start}:{end}" for start, end in
                   zip(group["3' UTR start"], group["3' UTR end"])]
    return ";".join(coordinates)


def translate_position_to_genomic(exon_structure, start_position, end_position, strand):
    exons = exon_structure.split(';')
    cumulative_length = 0
    start_genomic = None
    end_genomic = None
    intermediate_coords = []

    if strand == -1:
        exons = exons[::-1]  # Reverse the order of exons
        exons = [f"{end}:{start}" for start, end in (exon.split(':') for exon in exons)]  # Swap start and end
        for exon in exons:
            start, end = map(float, exon.split(':'))
            start = int(round(start))
            end = int(round(end))
            exon_length = start - end + 1
            if len(intermediate_coords) != 0:
                intermediate_coords.append((start , end))

            if start_genomic is None and cumulative_length + exon_length >= start_position:
                start_genomic = start - (start_position - cumulative_length + 1)

            if cumulative_length + exon_length >= end_position:
                end_genomic = start - (end_position - cumulative_length + 1)
                break

            if start_genomic is not None and end_genomic is None:
                # Collect intermediate exon boundaries for the case where positions span multiple exons
                intermediate_coords.append((start, end))

            cumulative_length += exon_length

        if start_genomic is None or end_genomic is None:
            raise ValueError("Positions out of range of the exon structure.")

        # Format the output based on whether positions are within the same exon or span multiple exons
        if len(intermediate_coords) == 0:
            return str(end_genomic), str(start_genomic)
        else:
            start_boundaries = []
            end_boundaries = []
            for n in range(1, len(intermediate_coords)):
                start_boundaries.append(str(intermediate_coords[n][0]))
                end_boundaries.append(str(intermediate_coords[n - 1][1]))
            start_positions_genomic = f"{','.join(end_boundaries)},{end_genomic}"
            start_positions_genomic_order = ','.join(start_positions_genomic.split(",")[::-1])
            end_positions_genomic = f"{start_genomic},{','.join(start_boundaries)}"
            end_positions_genomic_order = ','.join(end_positions_genomic.split(",")[::-1])
            return start_positions_genomic_order, end_positions_genomic_order


    else:
        for exon in exons:
            start, end = map(float, exon.split(':'))
            start = int(round(start))
            end = int(round(end))
            exon_length = end - start + 1

            if len(intermediate_coords) != 0:
                intermediate_coords.append((start, end))

            if start_genomic is None and cumulative_length + exon_length >= start_position:
                start_genomic = start + (start_position - cumulative_length - 1)

            if cumulative_length + exon_length >= end_position:
                end_genomic = start + (end_position - cumulative_length - 1)
                break

            if start_genomic is not None and end_genomic is None:
                # Collect intermediate exon boundaries for the case where positions span multiple exons
                intermediate_coords.append((start, end))

            cumulative_length += exon_length

        if start_genomic is None or end_genomic is None:
            raise ValueError("Positions out of range of the exon structure.")

        # Format the output based on whether positions are within the same exon or span multiple exons
        if len(intermediate_coords) == 0:
            return str(start_genomic), str(end_genomic)
        else:
            start_boundaries = []
            end_boundaries = []
            for n in range(1, len(intermediate_coords)):
                start_boundaries.append(str(intermediate_coords[n][0]))
                end_boundaries.append(str(intermediate_coords[n-1][1]))

            return f"{start_genomic},{','.join(start_boundaries)}", f"{','.join(end_boundaries)},{end_genomic}"



def shuffle_sequence(seq):
    return ''.join(random.sample(seq, len(seq)))

def deduplicate_list(input_list):
    seen = set()
    deduplicated_list = []
    for item in input_list:
        if item not in seen:
            deduplicated_list.append(item)
            seen.add(item)
    return deduplicated_list


def get_known_seeds(mirgenedb, mirbase):
    mirgenedb_records = SeqIO.parse(mirgenedb, "fasta")
    mirbase_records = SeqIO.parse(mirbase, "fasta")
    seeds = []

    for record in mirgenedb_records:
        seeds.append(record.seq[1:8])
    for record in mirbase_records:
        seeds.append(record.seq[1:8])

    return deduplicate_list(seeds)


def filter_and_shuffle_miRNAs(input_fasta, output_fasta, known_seeds):
    records = SeqIO.parse(input_fasta, "fasta")

    valid_records = []
    index = 1

    for original_record in records:
        shuffled_seq = shuffle_sequence(str(original_record.seq))
        seed = shuffled_seq[1:8]

        # Shuffle until a valid seed is found
        while seed in known_seeds:
            shuffled_seq = shuffle_sequence(str(original_record.seq))
            seed = shuffled_seq[1:8]

        record_id = f"iri-random-{index}"
        record = SeqRecord(Seq(shuffled_seq), id=record_id, description="")
        valid_records.append(record)
        index += 1

    SeqIO.write(valid_records, output_fasta, "fasta")

def target_filtered_to_seed_network(target_filtered, seed_tsv, output):
    df1 = pd.read_csv(target_filtered, sep='\t', header=None, usecols=[0, 1])

    # Split the gene name column to get the actual gene name
    df1[1] = df1[1].apply(lambda x: x.split('|')[0])

    # Count the unique miRNA-gene combinations
    df1_counts = df1.value_counts().reset_index(name='count')

    df2 = pd.read_csv(seed_tsv, sep='\t', header=None, names=['miRNA', 'seed'])

    # Replace miRNA names with seed names
    df1_counts[0] = df1_counts[0].map(df2.set_index('miRNA')['seed'])

    # Drop rows with NaN values after mapping
    df1_counts.dropna(inplace=True)

    # Combine counts, taking the highest count among duplicate combinations
    df1_combined = df1_counts.groupby([0, 1])['count'].max().reset_index()

    df1_combined.to_csv(output, sep="\t", index=False)
