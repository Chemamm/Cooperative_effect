import argparse
import Cooperative

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Conserved combinations file")
parser.add_argument("-l", "--library", help="miRNAs library")
parser.add_argument("-n", "--number", type=int, help="number of genes to extract")
parser.add_argument("-e", "--ematrix", help="Expression matrix")
parser.add_argument("-o", "--output_label", help="Output file")
args = parser.parse_args()


def main():
    gene_count_dict_2, gene_count_dict_3, gene_mirnas_dict = Cooperative.get_gene_count(args.file)
    Cooperative.write_tsv(gene_count_dict_2, gene_count_dict_3, args.output_label)
    mirnas = Cooperative.write_top_mirnas_and_genes(gene_count_dict_2, gene_mirnas_dict, args.output_label, args.number)
    Cooperative.expression_analysis(mirnas, args.ematrix, args.library, args.output_label, args.number)


if __name__ == "__main__":
    main()




