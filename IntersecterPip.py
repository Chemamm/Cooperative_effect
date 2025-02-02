#!/usr/bin/env/python3

import argparse
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import random
import os


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


def pertranscript_dicter(input):
    pertranscript_dict = dict()
    for i in range(1, len(input)):
        line = input[i]
        pertranscript_dict[line.split("\t")[0]] = line.split("\t")[1].split(",")
    return pertranscript_dict


def fasta_to_dict(ft):
    fasta_dict = dict()
    seqid = ""
    for i in range(0, len(ft)):
        if ft[i].startswith('>'):
            if seqid != "":
                fasta_dict[seqid] = seq
            seqid = ft[i].split(' ')[0].replace(">", "").replace("\n", "")
            seq = ""
        else:
            seq += ft[i].replace("\n", "")
    if seqid != "":
        fasta_dict[seqid] = seq
    return fasta_dict


parser = argparse.ArgumentParser()
parser.add_argument("-p", "--pertranscript", default=False, help="pertranscript.txt file")
parser.add_argument("-mp", "--mirconstarget_path", help="path to the folder that cointains the .jar")
parser.add_argument("-t", "--target", help="target region")
parser.add_argument("-m", "--miRNA", help="miRNA file")
parser.add_argument("-th", "--threads", help="number of threads to use")
parser.add_argument("-mo", "--mirconstarget_output", help="mirconstarget output folder")
parser.add_argument("-a", "--algorythms", help="algorythms to use for mirconstarget. i.e SEED:MIRANDA")
parser.add_argument("-c", "--consensus", help="number of algorythms that must predict the target region")
parser.add_argument("-sp", "--specie", help="miRNA specie")
parser.add_argument("-rs", "--random_seed", help="Number of random seeds u want to chose")
parser.add_argument("-la", "--label", help="label for output files (includes output folder if needed)")
parser.add_argument("-ip", "--intersecter_path", help="path to the folder that contains intersecter.jar")
parser.add_argument("-d", "--depth", default=10, help="max set size to explore")
parser.add_argument("-tp", "--top", default=10, help="top number of sets to show")
parser.add_argument("-ra", "--ratioSort", default=False, action="store_true", help="top number of sets to show")

args = parser.parse_args()


# MIRCONSTARGET
if args.mirconstarget_path:
    mirconstarget_cmd = "java -jar %s/miRNAconsTargets.jar %s %s %s %s %s" % (args.mirconstarget_path, args.miRNA, args.target, args.mirconstarget_output, args.threads, args.algorythms)
    os.system(mirconstarget_cmd)

# GET CONSENSUS
if args.consensus:
    files = args.algorythms.replace(":", ",").replace("PITA", "pita.txt").replace("TS", "ts.txt").replace("MIRANDA", "miranda.txt").replace("SEED", "seed2-8.txt")
    consensus_cmd = "java -jar %s/makeTargetPositionalConsensus.jar %s %s %s %s 7" % (args.mirconstarget_path, args.mirconstarget_output, files, args.target, args.consensus)
    os.system(consensus_cmd)
    os.system("mv %s/positionalConsensus_1_7.txt %s/seed2-8_positionalConsensus_1_7.txt" % (args.mirconstarget_output, args.mirconstarget_output))
    os.system("mv %s/perTranscript_1_7.txt %s/seed2-8_perTranscript_1_7.txt" % (args.mirconstarget_output, args.mirconstarget_output))
    os.system("mv %s/multipleTargets_1_7.txt %s/seed2-8_multipleTargets_1_7.txt" % (args.mirconstarget_output, args.mirconstarget_output))
    files = args.algorythms.replace(":", ",").replace("PITA", "pita.txt").replace("TS", "ts.txt").replace("MIRANDA", "miranda.txt").replace("SEED", "seed2-8A.txt")
    consensus_cmd = "java -jar %s/makeTargetPositionalConsensus.jar %s %s %s %s 7" % (args.mirconstarget_path, args.mirconstarget_output, files, args.target, args.consensus)
    os.system(consensus_cmd)
    os.system("mv %s/positionalConsensus_1_7.txt %s/seed2-8A_positionalConsensus_1_7.txt" % (args.mirconstarget_output, args.mirconstarget_output))
    os.system("mv %s/perTranscript_1_7.txt %s/seed2-8A_perTranscript_1_7.txt" % (args.mirconstarget_output, args.mirconstarget_output))
    os.system("mv %s/multipleTargets_1_7.txt %s/seed2-8A_multipleTargets_1_7.txt" % (args.mirconstarget_output, args.mirconstarget_output))

# GETTING SEEDS
if args.pertranscript:
    pertranscript = list(open(args.pertranscript))
    pertranscript_dict = pertranscript_dicter(pertranscript)
    tmp_name = "%s_tmp.csv" % (args.label)
    with open(tmp_name, "wt") as tmp_out:
        for item in pertranscript_dict.keys():
            for value in pertranscript_dict[item]:
                tmp_out.write("%s,%s,1\n" % (item, value))
    names = [tmp_name]
elif not args.pertranscript:
    pertranscript = list(open("%s/seed2-8_perTranscript_1_7.txt" % args.mirconstarget_output))
    pertranscript_dict = pertranscript_dicter(pertranscript)
    tmp_name = "%s/%s_seed_2-8_tmp.csv" % (args.mirconstarget_output, args.label)
    with open(tmp_name, "wt") as tmp_out:
        for item in pertranscript_dict.keys():
            for value in pertranscript_dict[item]:
                tmp_out.write("%s,%s,1\n" % (item, value))
    names = [tmp_name]
    pertranscript = list(open("%s/seed2-8_perTranscript_1_7.txt" % args.mirconstarget_output))
    pertranscript_dict = pertranscript_dicter(pertranscript)
    tmp_name = "%s/%s_seed_2-8A_tmp.csv" % (args.mirconstarget_output, args.label)
    with open(tmp_name, "wt") as tmp_out:
        for item in pertranscript_dict.keys():
            for value in pertranscript_dict[item]:
                tmp_out.write("%s,%s,1\n" % (item, value))
    names.append(tmp_name)

uniq_names = []
# Modificando para seeds.
for item in names:

    if args.miRNA:
        fasta = list(open(args.miRNA))
        mirna_dict = fasta_to_dict(fasta)

        seed_dict = dict()
        coding_dict = dict()
        count = 1

        # Creando dos diccionarios, uno que contenga la lista de semillas, y otro que contenga que miRNAs se corresponden con cada semilla.

        for key in mirna_dict:
            seed = mirna_dict[key][1:8]
            if seed not in seed_dict.keys():
                seed_dict[seed] = "%s_seed_%s" % (args.specie, str(count))
                coding_dict[key] = seed_dict[seed]
                count += 1
            else:
                coding_dict[key] = seed_dict[seed]

        if os.path.isfile("%s.seed.coding" % (args.miRNA)) == False and os.path.isfile("%s.seed" % (args.miRNA)) == False:
            # Escribiendo los ficheros de salida.
            with open("%s.seed" % args.miRNA, "wt") as out1:
                for key in seed_dict:
                    out1.write(">%s\n%s\n" % (seed_dict[key], key))
            with open("%s.seed.coding" % (args.miRNA), "wt") as out2:
                for key in coding_dict:
                    out2.write("%s\t%s\n" % (key, coding_dict[key]))
            # Modificando el csv
        with open(item) as tmp_file:
            tmp_str = tmp_file.read()
            # Para tomar un numero de semillas aleatorias de la muestra
            if args.random_seed:
                random_sample = random.sample(list(coding_dict.values()), int(args.random_seed))
                for key, seed in coding_dict.items():
                    if seed in random_sample:
                        tmp_str = tmp_str.replace("%s," % key, "%s," % seed)
            else:
                for key in coding_dict:
                    tmp_str = tmp_str.replace("%s," % key, "%s," % coding_dict[key])
        tmp_name = "%s.parsed" % (item)
        with open(tmp_name, "wt") as out:
            out.write(tmp_str)
        if args.random_seed:
            delete_non_seed_cmd = "sed -i '/seed/!d' %s" % tmp_name
            os.system(delete_non_seed_cmd)
        cmd = "sort -u -o %s_unique %s" % (tmp_name, tmp_name)
        os.system(cmd)
        tmp_name = "%s_unique" % (tmp_name)

    os.system("sed -i '1i gene,miRNA,boolean' %s" % (tmp_name))
    uniq_names.append(tmp_name)

# INTERSECTER
os.system("mkdir %s/intersecter" % args.mirconstarget_output)
for item in uniq_names:
    intersecter_outpath = "%s/intersecter/%s" % (args.mirconstarget_output, item.split("/")[-1])
    if args.ratioSort:
        intersecter_cmd = "java -jar %s/intersecter.jar input=%s output=%s mode=ap percentile=0.01 keyCol=1 valueCol=0 sep=, maxLevel=%s aprioriMode=ratioSort" % (args.intersecter_path, item, intersecter_outpath, args.depth)
    else:
        intersecter_cmd = "java -jar %s/intersecter.jar input=%s output=%s mode=ap percentile=0.01 keyCol=1 valueCol=0 sep=, maxLevel=%s" % (args.intersecter_path, item, intersecter_outpath, args.depth)
    os.system(intersecter_cmd)
    intersecter_outfile = "%s/random_apriori.tsv" % intersecter_outpath
    intersecter_label = "%s/random_apriori" % intersecter_outpath
# INTERSECTER_TOPS
    if not args.ratioSort:
        df = pd.read_csv(intersecter_outfile, sep="\t")

        intersection_global = Set(df["observed intersection"])

        sets = dict()

        for n in range(1, args.depth + 1):
            df_tmp = df[(df['set size'] == n)]
            if not df_tmp.empty:
                sets["set_%s" % (n)] = Set(df_tmp["observed intersection"])
            else:
                sets["set_%s" % (n)] = Empty_Set()

        with open("%s_global_stats.txt" % (intersecter_label), "wt") as out:
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
        fig.savefig("%s_per_size_distribution.pdf" % (intersecter_label))

        fig = plt.figure()
        label = "Intersection Distribution"
        sns.distplot(intersection_global.df, axlabel=label)
        fig.savefig("%s.intersectionplot.pdf" % (intersecter_label))

        if args.top:
            coding_file = "%s.seed.coding" % args.miRNA
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
                        seed_genes[",".join(seeds)] = getgenes(seeds, item)
                        if counter == int(args.top):
                            break
                    with open("%s_top_%s_setsize_%s" % (intersecter_label, args.top, n), "wt") as out:
                        out.write("intersection\tseeds\tmiRNAs\tgenes\n")
                        for key in seed_intersections:
                            out.write("%s\t%s\t%s\t%s\n" % (seed_intersections[key], key, ",".join(seed_mirnas[key]), ",".join(seed_genes[key])))
    else:
        df = pd.read_csv(intersecter_outfile, sep="\t")

        intersection_global = Set(df["obs-exp ratio"])

        sets = dict()

        for n in range(1, args.depth + 1):
            df_tmp = df[(df['set size'] == n)]
            if not df_tmp.empty:
                sets["set_%s" % (n)] = Set(df_tmp["obs-exp ratio"])
            else:
                sets["set_%s" % (n)] = Empty_Set()

        with open("%s_global_stats.txt" % (intersecter_label), "wt") as out:
            out.write("Set_size\tnumber_of_itemsets\tmax_obs-exp_ratio\tp99\tp95\tp80\tp50\n")
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
        fig.savefig("%s_per_size_distribution.pdf" % (intersecter_label))

        fig = plt.figure()
        label = "Obs-exp ratio Distribution"
        sns.distplot(intersection_global.df, axlabel=label)
        fig.savefig("%s.intersectionplot.pdf" % (intersecter_label))

        if args.top:
            coding_file = "%s.seed.coding" % args.miRNA
            for n in range(1, args.depth + 1):
                seed_intersections = dict()
                seed_mirnas = dict()
                seed_genes = dict()
                df_tmp = df[(df['set size'] == n)][["obs-exp ratio", "set"]]
                if df_tmp.empty == False:
                    df_by_intersection = df_tmp.sort_values("obs-exp ratio", ascending=False)
                    counter = 1
                    print(df_by_intersection)
                    for index, row in df_by_intersection.iterrows():
                        intersection = row["obs-exp ratio"]
                        seeds = row["set"].replace(" ", "").split(",")
                        counter += 1
                        seed_intersections[",".join(seeds)] = intersection
                        seed_mirnas[",".join(seeds)] = seed2mirna(seeds, coding_file)
                        seed_genes[",".join(seeds)] = getgenes(seeds, item)
                        if counter == int(args.top):
                            break
                    with open("%s_top_%s_setsize_%s" % (intersecter_label, args.top, n), "wt") as out:
                        out.write("obs-exp_ratio\tseeds\tmiRNAs\tgenes\n")
                        for key in seed_intersections:
                            out.write("%s\t%s\t%s\t%s\n" % (seed_intersections[key], key, ",".join(seed_mirnas[key]), ",".join(seed_genes[key])))
