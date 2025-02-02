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
parser.add_argument("-r", "--rules", help="rules_lift filename")
parser.add_argument("-f", "--frequent", help="support filename")
parser.add_argument("-l", "--label", help="label of the output_files")
parser.add_argument("-t", "--top", default=False, help="top number of sets to show")
parser.add_argument("-c", "--csv", default=False, help="__tmp.csv.parsed_unique")
parser.add_argument("-m", "--mirna", default=False, help="mirna file")
args = parser.parse_args()

df_rules = pd.read_csv(args.rules, sep="\t")
df_frequents = pd.read_csv(args.frequent, sep="\t")
df_frequents['size'] = df_frequents['itemsets'].apply(lambda x: len(x.split(",")))
df_rules['size'] = df_rules["antecedents"].apply(lambda x: len(x.split(","))) + df_rules["consequents"].apply(lambda x: len(x.split(",")))

supp_global = Set(df_frequents["support"])
lift_global = Set(df_rules["lift"])

supp_sets = dict()
lift_sets = dict()
print(df_rules)
for n in range(1, 7):
    df_freq_tmp = df_frequents[(df_frequents['size'] == n)]
    df_rul_tmp = df_rules[(df_rules['size'] == n)]
    if not df_freq_tmp.empty:
        supp_sets["supp_set_%s" % (n)] = Set(df_freq_tmp["support"])
    else:
        supp_sets["supp_set_%s" % (n)] = Empty_Set()
    if not df_rul_tmp.empty:
        lift_sets["lift_set_%s" % (n)] = Set(df_rul_tmp["lift"])
    else:
        lift_sets["lift_set_%s" % (n)] = Empty_Set()

with open("%s_global_stats.txt" % (args.label), "wt") as out:
    out.write("Set_size\tnumber_of_supp_itemsets\tmax_supp\tp99_supp\tp95_supp\tp80_supp\tp50_supp\tnumber_of_lift_relations\tmax_lift\tp99_lift\tp95_lift\tp80_lift\tp50_lift\n")
    out.write("*\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
        supp_global.n, supp_global.max, supp_global.p99, supp_global.p95, supp_global.p80, supp_global.p50, lift_global.n,
        lift_global.max, lift_global.p99, lift_global.p95, lift_global.p80, lift_global.p50))
    for n in range(1, 7):
        out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            n, supp_sets["supp_set_%s" % (n)].n, supp_sets["supp_set_%s" % (n)].max, supp_sets["supp_set_%s" % (n)].p99,
            supp_sets["supp_set_%s" % (n)].p95, supp_sets["supp_set_%s" % (n)].p80, supp_sets["supp_set_%s" % (n)].p50,
            lift_sets["lift_set_%s" % (n)].n, lift_sets["lift_set_%s" % (n)].max, lift_sets["lift_set_%s" % (n)].p99,
            lift_sets["lift_set_%s" % (n)].p95, lift_sets["lift_set_%s" % (n)].p80, lift_sets["lift_set_%s" % (n)].p50))


fig = plt.figure()
fig.subplots_adjust(hspace=0.4, wspace=0.4)
sns.set(color_codes=True)
for n in range(1, 7):
    ax = fig.add_subplot(2, 3, n)
    if not type(supp_sets["supp_set_%s" % (n)].df) == str and len(supp_sets["supp_set_%s" % (n)].df.to_numpy()) != 1:
        sns.distplot(supp_sets["supp_set_%s" % (n)].df.to_numpy(), axlabel="Set size %s" % (n), ax=ax)
    elif not type(supp_sets["supp_set_%s" % (n)].df) == str and len(supp_sets["supp_set_%s" % (n)].df.to_numpy()) == 1:
        plt.hist(supp_sets["supp_set_%s" % (n)].df.to_numpy(), bins=1)
    else:
        plt.text(0.5, 0.5, "nan", fontsize=18, ha='center')
        plt.xlabel("Set size %s" % (n))
fig.savefig("%s_supp_per_size_distribution.pdf" % (args.label))

fig = plt.figure()
fig.subplots_adjust(hspace=0.4, wspace=0.4)
sns.set(color_codes=True)
ax = fig.add_subplot(2, 3, 1)
plt.text(0.5, 0.5, "nan",
         fontsize=18, ha='center')
for n in range(2, 7):
    ax = fig.add_subplot(2, 3, n)
    if not type(lift_sets["lift_set_%s" % (n)].df) == str and len(lift_sets["lift_set_%s" % (n)].df.to_numpy()) != 1:
        sns.distplot(lift_sets["lift_set_%s" % (n)].df.to_numpy(), axlabel="Set size %s" % (n), ax=ax)
    elif not type(lift_sets["lift_set_%s" % (n)].df) == str and len(lift_sets["lift_set_%s" % (n)].df.to_numpy()) == 1:
        plt.hist(lift_sets["lift_set_%s" % (n)].df.to_numpy(), bins=1)
    else:
        plt.text(0.5, 0.5, "nan", fontsize=18, ha='center')
        plt.xlabel("Set size %s" % (n))
fig.savefig("%s_lift_per_size_distribution.pdf" % (args.label))

fig = plt.figure()
label = "%s Support Distribution" % args.label
sns.distplot(supp_global.df, axlabel=label)
fig.savefig("%s.suppplot.pdf" % (args.label))

fig = plt.figure()
label = "%s Lift Distribution" % args.label
sns.distplot(lift_global.df, axlabel=label)
fig.savefig("%s.liftplot.pdf" % (args.label))

if args.top:
    coding_file = "%s.seed.coding" % args.mirna
    for n in range(1, 7):
        seed_supps = dict()
        seed_mirnas = dict()
        seed_genes = dict()
        seed_lifts = dict()
        df_freq_tmp = df_frequents[(df_frequents['size'] == n)][["support", "itemsets"]]
        df_rul_tmp = df_rules[(df_rules['size'] == n)][["lift", "antecedents", "consequents"]]
        if df_freq_tmp.empty == False:
            df_by_supp = df_freq_tmp.sort_values("support", ascending=False)
            counter = 1
            print(df_by_supp)
            for index, row in df_by_supp.iterrows():
                supp = row["support"]
                seeds = row["itemsets"].replace(" ", "").split(",")
                seeds2 = row["itemsets"].replace(" ", "").replace("'", "").split(",")
                counter += 1
                seed_supps[",".join(seeds)] = supp
                seed_mirnas[",".join(seeds)] = seed2mirna(seeds, coding_file)
                seed_genes[",".join(seeds)] = getgenes(seeds2, args.csv)
                if counter == int(args.top):
                    break
            with open("%s_top_%s_bysupp_setsize_%s" % (args.label, args.top, n), "wt") as out:
                out.write("support\tseeds\tmiRNAs\tgenes\n")
                for key in seed_supps:
                    out.write("%s\t%s\t%s\t%s\n" % (seed_supps[key], key, ",".join(seed_mirnas[key]), ",".join(seed_genes[key])))
            seed_mirnas = dict()
            seed_genes = dict()
        if df_rul_tmp.empty == False:
            df_by_supp = df_rul_tmp.sort_values("lift", ascending=False)
            counter = 1
            dup = False
            for index, row in df_by_supp.iterrows():
                lift = row["lift"]
                seeds = row["antecedents"].replace(" ", "").split(",") + row["consequents"].replace(" ", "").split(",")
                seeds2 = row["antecedents"].replace(" ", "").replace("'", "").split(",") + row["consequents"].replace(" ", "").replace("'", "").split(",")
                seeds.sort()
                seeds2.sort()
                if ",".join(seeds) in seed_lifts.keys():
                    seed_lifts[",".join(seeds)].append(float(lift))
                else:
                    counter += 1
                    seed_lifts[",".join(seeds)] = [float(lift)]
                    seed_mirnas[",".join(seeds)] = seed2mirna(seeds, coding_file)
                    seed_genes[",".join(seeds)] = getgenes(seeds2, args.csv)
                if counter == int(args.top):
                    break
            with open("%s_top_%s_bylift_setsize_%s" % (args.label, args.top, n), "wt") as out:
                out.write("medium_lift\tseeds\tmiRNAs\tgenes\n")
                for key in seed_lifts:
                    out.write("%s\t%s\t%s\t%s\n" % (
                        np.mean(seed_lifts[key]), key, ",".join(seed_mirnas[key]), ",".join(seed_genes[key])))
            seed_mirnas = dict()
            seed_genes = dict()
