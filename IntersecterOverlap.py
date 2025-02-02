#!/usr/bin/env/python3

import argparse


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="intersecter top file")
parser.add_argument("-t", "--top", type=int, help="number of miRNA sets to study")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()


with open(args.file) as fl:
    fl = list(fl)
    prime_seed = fl[1].split("\t")[1].split(",")
    prime_mirna = fl[1].split("\t")[2].split(",")
    prime_gene = fl[1].split("\t")[3].split(",")
    with open(args.output, "wt") as out:
        out.write("seeds\tmiRNAs\tgenes\tseed overlap\tmirna overlap\tgene overlap\toverlapping seeds\toverlapping miRNAs\toverlapping genes\n")
        for n in range(1, args.top):
            seed = fl[n].split("\t")[1].split(",")
            mirna = fl[n].split("\t")[2].split(",")
            gene = fl[n].split("\t")[3].split(",")
            overlapping_seed = intersection(prime_seed, seed)
            overlapping_mirna = intersection(prime_mirna, mirna)
            overlapping_gene = intersection(prime_gene, gene)
            out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (",".join(seed), ",".join(mirna), ",".join(gene), len(overlapping_seed), len(overlapping_mirna), len(overlapping_gene), ",".join(overlapping_seed), ",".join(overlapping_mirna), ",".join(overlapping_gene)))
