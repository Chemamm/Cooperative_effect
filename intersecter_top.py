#!/usr/bin/env/python3

import argparse
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


class Set:
    def __init__(self, df):
        self.df = df
        self.list = df.to_numpy()
        self.max = np.max(self.list)
        self.p99 = np.percentile(self.list, 99)
        self.p95 = np.percentile(self.list, 95)
        self.p80 = np.percentile(self.list, 80)
        self.p50 = np.percentile(self.list, 50)
        self.n = len(self.list)


class Empty_Set:
    def __init__(self):
        self.df = "nan"
        self.list = "nan"
        self.max = "nan"
        self.p99 = "nan"
        self.p95 = "nan"
        self.p80 = "nan"
        self.p50 = "nan"
        self.n = "nan"


def seed2mirna(seeds, coding_file):
    with open(coding_file) as cod:
        mirna_list = list()
        cod = list(cod)
        for seed in seeds:
            seed = "'%s'" % seed
            for line in cod:
                mirna = line.split("\t")[0]
                seed_cod = "'%s'" % line.split("\t")[1].replace("\n", "")
                if seed.replace("''", "") == seed_cod:
                    mirna_list.append(mirna)
    return mirna_list


def getgenes(seeds, csv):
    with open(csv) as pt:
        pt = list(pt)
        gene_dict = dict()
        genes = list()
        for line in pt:
            gene = line.split(",")[0]
            seed = line.split(",")[1]
            if gene not in gene_dict.keys():
                gene_dict[gene] = seed
            else:
                gene_dict[gene] += ","
                gene_dict[gene] += seed
        for key in gene_dict:
            if all(item in gene_dict[key].split(",") for item in seeds):
                genes.append(key)
    return genes


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="intersecter outfilename")
parser.add_argument("-l", "--label", help="label of the output_files")
parser.add_argument("-t", "--top", default=False, help="top number of sets to show")
parser.add_argument("-c", "--csv", default=False, help="__tmp.csv.parsed_unique")
parser.add_argument("-m", "--mirna", default=False, help="mirna file")
parser.add_argument("-d", "--depth", default=10, help="max set size to explore")
args = parser.parse_args()

df = pd.read_csv(args.file, sep="\t")

intersection_global = Set(df["observed intersection"])

sets = dict()

for n in range(1, args.depth + 1):
    df_tmp = df[(df['set size'] == n)]
    if not df_tmp.empty:
        sets["set_%s" % (n)] = Set(df_tmp["observed intersection"])
    else:
        sets["set_%s" % (n)] = Empty_Set()

with open("%s_global_stats.txt" % (args.label), "wt") as out:
    out.write("Set_size\tnumber_of_itemsets\tmax_intersection\tp99\tp95\tp80\tp50\n")
    out.write("*\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
        intersection_global.n, intersection_global.max, intersection_global.p99, intersection_global.p95, intersection_global.p80, intersection_global.p50))
    for n in range(1, args.depth + 1):
        out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            n, sets["set_%s" % (n)].n, sets["set_%s" % (n)].max, sets["set_%s" % (n)].p99,
            sets["set_%s" % (n)].p95, sets["set_%s" % (n)].p80, sets["set_%s" % (n)].p50))


fig = plt.figure()
fig.subplots_adjust(hspace=0.4, wspace=0.4)
sns.set(color_codes=True)
for n in range(1, args.depth + 1):
    ax = fig.add_subplot(3, 4, n)
    if not type(sets["set_%s" % (n)].df) == str and len(sets["set_%s" % (n)].df.to_numpy()) != 1:
        sns.distplot(sets["set_%s" % (n)].df.to_numpy(), axlabel="Set size %s" % (n), ax=ax)
    elif not type(sets["set_%s" % (n)].df) == str and len(sets["set_%s" % (n)].df.to_numpy()) == 1:
        plt.hist(sets["set_%s" % (n)].df.to_numpy(), bins=1)
    else:
        plt.text(0.5, 0.5, "nan", fontsize=18, ha='center')
        plt.xlabel("Set size %s" % (n))
fig.savefig("%s_per_size_distribution.pdf" % (args.label))


fig = plt.figure()
label = "%s Intersection Distribution" % args.label
sns.distplot(intersection_global.df, axlabel=label)
fig.savefig("%s.intersectionplot.pdf" % (args.label))

if args.top:
    coding_file = "%s.seed.coding" % args.mirna
    for n in range(1, args.depth + 1):
        seed_intersections = dict()
        seed_mirnas = dict()
        seed_genes = dict()
        df_tmp = df[(df['set size'] == n)][["observed intersection", "set"]]
        if df_tmp.empty == False:
            df_by_intersection = df_tmp.sort_values("observed intersection", ascending=False)
            counter = 1
            print(df_by_intersection)
            for index, row in df_by_intersection.iterrows():
                intersection = row["observed intersection"]
                seeds = row["set"].replace(" ", "").split(",")
                counter += 1
                seed_intersections[",".join(seeds)] = intersection
                seed_mirnas[",".join(seeds)] = seed2mirna(seeds, coding_file)
                seed_genes[",".join(seeds)] = getgenes(seeds, args.csv)
                if counter == int(args.top):
                    break
            with open("%s_top_%s_setsize_%s" % (args.label, args.top, n), "wt") as out:
                out.write("intersection\tseeds\tmiRNAs\tgenes\n")
                for key in seed_intersections:
                    out.write("%s\t%s\t%s\t%s\n" % (seed_intersections[key], key, ",".join(seed_mirnas[key]), ",".join(seed_genes[key])))
