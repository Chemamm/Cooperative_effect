import argparse
import Cooperative
import seaborn as sns
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Positional Consensus file")
parser.add_argument("-sa", "--saliva", help="Saliva specific miRNAs")
parser.add_argument("-sg", "--sg", help="Saliva specific miRNAs")
parser.add_argument("-c", "--combinations", help="Valid combinations file")
parser.add_argument("-s", "--species", help="Positional consensus files from other species")
parser.add_argument("-o", "--output_label", help="Output label")
parser.add_argument("-t", "--threshold", type=int, help="Threshold of conservation")
args = parser.parse_args()

with open(args.saliva) as sal:
    saliva_mirnas = [x.strip("\n").strip(" ") for x in list(sal)]
with open(args.sg) as sg:
    sg_mirnas = [x.strip("\n").strip(" ") for x in list(sg)]

def get_enrichment(miRNA):
    if miRNA in saliva_mirnas:
        return 'Saliva'
    elif miRNA in sg_mirnas:
        return 'SG'
    else:
        return 'None'


def saliva_catplot(df, output):
    output_ratio = output + "_ratio.png"
    output_raw = output + "_raw.png"
    output_txt = output + "_stats.txt"
    output_conserved_targets = output + "_conserved_targets.png"
    df_filtered = df[df['Enrichment'] != 'None']

    df_filtered.to_csv(output_txt, sep="\t")

    sns.catplot(data=df_filtered, x="Enrichment", y="ratio_2")
    plt.savefig(output_ratio, dpi=300)
    plt.close()

    sns.catplot(data=df_filtered, x="Enrichment", y="targets_valid_2_count")
    plt.savefig(output_raw, dpi=300)
    plt.close()

    sns.catplot(data=df_filtered, x="Enrichment", y="conserved_targets_count")
    plt.savefig(output_conserved_targets, dpi=300)
    plt.close()



def main():
    df = Cooperative.get_cooperative_targets(args.file, args.combinations, args.species, args.threshold)
    df["Enrichment"] = df["miRNA"].apply(get_enrichment)
    saliva_catplot(df, args.output_label)

if __name__ == "__main__":
    main()


