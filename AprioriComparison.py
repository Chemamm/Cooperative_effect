#!/usr/bin/env/python3

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob
import os
import collections


def pad_dict_list(dict_list, padel):
    lmax = 0
    for lname in dict_list.keys():
        lmax = max(lmax, len(dict_list[lname]))
    for lname in dict_list.keys():
        ll = len(dict_list[lname])
        if ll < lmax:
            dict_list[lname] += [padel] * (lmax - ll)
    return dict_list


def match(s1, s2):
    ok = False

    for c1, c2 in zip(s1, s2):
        if c1 != c2:
            if ok:
                return False
            else:
                ok = True
    return ok


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="tsv with folder/file\tlabel\n")
parser.add_argument("-l", "--label", help="initial part of files to open")
parser.add_argument("-o", "--output", help="output name")
parser.add_argument("-d", "--depth", default=10, help="max set size to explore")
parser.add_argument("-to", "--top", default=10, help="numer of sets in the top")
parser.add_argument("-fo", "--folder", help="mirconstarget outfolder form IntersecterPip")
args = parser.parse_args()


if args.folder:
    # seed2-8
    intersections_max = dict()
    folder = glob("%s/intersecter/*seed_2-8_*" % args.folder)
    label = folder[0].replace("'", "").split("/")[-1].split("_seed")[0]
    intersections_max[label] = []
    for n in range(2, args.depth + 1):
        if os.path.isfile("%s/random_apriori_top_%s_setsize_%s" % (folder[0].replace("'", ""), args.top, n)):
            with open("%s/random_apriori_top_%s_setsize_%s" % (folder[0].replace("'", ""), args.top, n)) as ra:
                ra = list(ra)
                max_inter = ra[1].split("\t")[0]
                intersections_max[label].append(float(max_inter))
    folder_list = glob("%s/*/" % args.folder)
    folder_list.remove('%s/intersecter/' % args.folder)
    for fd in folder_list:
        folder = glob("%s/intersecter/*seed_2-8_*" % fd)
        label = folder[0].replace("'", "").split("/")[-1].split("_seed")[0]
        intersections_max[label] = []
        for n in range(2, args.depth + 1):
            if os.path.isfile("%s/random_apriori_top_%s_setsize_%s" % (folder[0].replace("'", ""), args.top, n)):
                with open("%s/random_apriori_top_%s_setsize_%s" % (folder[0].replace("'", ""), args.top, n)) as ra:
                    ra = list(ra)
                    max_inter = ra[1].split("\t")[0]
                    intersections_max[label].append(float(max_inter))
    intersections_max = pad_dict_list(intersections_max, 0)
else:
    with open(args.file) as filelist:
        filelist = list(filelist)
        intersections_max = dict()
        intersections_grouped = dict()
        for line in filelist:
            folder = line.split("\t")[0]
            label = line.split("\t")[1].replace("\n", "")
            intersections_max[label] = []
            for n in range(2, args.depth + 1):
                with open("%s/%s%s" % (folder, args.label, n)) as ra:
                    ra = list(ra)
                    max_inter = ra[1].split("\t")[0]
                    intersections_max[label].append(float(max_inter))


index = range(2, len(intersections_max[label]) + 2)
index = [int(x) for x in index]
intersections_max = collections.OrderedDict(sorted(intersections_max.items()))
df = pd.DataFrame.from_dict(intersections_max)
df["set_size"] = index
df = df.set_index("set_size")
df = df.sort_index()
print(df)

tsv = "%s_seed2-8.tsv" % args.output
df.to_csv(path_or_buf=tsv, sep="\t")

fig = plt.figure()
colors = ["#d64b4b", "#d6934b", "#e3d46f", "#bce36f", "#75e36f", "#6fe3af", "#6fe3dd", "#6f96e3", "#3c30db", "#9964d1", "#cf64d1", "#d16495", "#736e70", "#000000", "#f58700", "#00f521"]
sns.set_palette(sns.color_palette(colors))
sns.set(style="whitegrid")
sns.lineplot(data=df, linewidth=1.5, dashes=False)
fig.savefig("%s_seed2-8.pdf" % args.output)


if args.folder:
    # seed2-8A
    intersections_max = dict()
    folder = glob("%s/intersecter/*seed_2-8A*" % args.folder)
    label = folder[0].replace("'", "").split("/")[-1].split("_seed")[0]
    intersections_max[label] = []
    for n in range(2, args.depth + 1):
        if os.path.isfile("%s/random_apriori_top_%s_setsize_%s" % (folder[0].replace("'", ""), args.top, n)):
            with open("%s/random_apriori_top_%s_setsize_%s" % (folder[0].replace("'", ""), args.top, n)) as ra:
                ra = list(ra)
                max_inter = ra[1].split("\t")[0]
                intersections_max[label].append(float(max_inter))
    folder_list = glob("%s/*/" % args.folder)
    folder_list.remove('%s/intersecter/' % args.folder)
    for fd in folder_list:
        folder = glob("%s/intersecter/*seed_2-8A*" % fd)
        label = folder[0].replace("'", "").split("/")[-1].split("_seed")[0]
        intersections_max[label] = []
        for n in range(2, args.depth + 1):
            if os.path.isfile("%s/random_apriori_top_%s_setsize_%s" % (folder[0].replace("'", ""), args.top, n)):
                with open("%s/random_apriori_top_%s_setsize_%s" % (folder[0].replace("'", ""), args.top, n)) as ra:
                    ra = list(ra)
                    max_inter = ra[1].split("\t")[0]
                    intersections_max[label].append(float(max_inter))
    intersections_max = pad_dict_list(intersections_max, 0)


index2 = range(2, len(intersections_max[label]) + 2)
index2 = [int(x) for x in index2]
intersections_max = collections.OrderedDict(sorted(intersections_max.items()))
df2 = pd.DataFrame.from_dict(intersections_max)
df2["set_size"] = index2
df2 = df2.set_index("set_size")
df2 = df2.sort_index()
print(df2)

tsv = "%s_seed2-8A.tsv" % args.output
df2.to_csv(path_or_buf=tsv, sep="\t")

fig = plt.figure()
colors = ["#d64b4b", "#d6934b", "#e3d46f", "#bce36f", "#75e36f", "#6fe3af", "#6fe3dd", "#6f96e3", "#3c30db", "#9964d1", "#cf64d1", "#d16495", "#736e70", "#000000", "#f58700", "#00f521"]
sns.set_palette(sns.color_palette(colors))
sns.set(style="whitegrid")
sns.lineplot(data=df2, linewidth=1.5, dashes=False)
fig.savefig("%s_seed2-8A.pdf" % args.output)


# Obtaining means

with open("%s_seed2-8.tsv" % args.output) as tsv:
    tsv_lines = list(tsv)
    count = -1
    mean_groups = {}
    mean_groups[0] = []
    item_before = "laksjd"
    for item in tsv_lines[0].replace("\n", "").split("\t")[1::]:
        if match(item, item_before):
            mean_groups[count].append(item)
        else:
            count += 1
            mean_groups[count] = [item]
        item_before = item
mean_dict = {}
data_dict = {}

for key in mean_groups:
    ID = mean_groups[key][0].replace("_1", "")
    data_dict[ID] = []
    mean_dict[ID] = []
    for item in mean_groups[key]:
        item_list = df[item].tolist()
        data_dict[ID].append(item_list)
    for n in range(0, len(item_list)):
        mean_list = []
        for j in range(0, len(data_dict[ID])):
            mean_list.append(data_dict[ID][j][n])
        mean_dict[ID].append(np.mean(mean_list))

df_mean = pd.DataFrame.from_dict(mean_dict)
df_mean["set_size"] = index
df_mean = df_mean.set_index("set_size")
df_mean = df_mean.sort_index()
print(df_mean)
tsv = "%s_seed2-8.mean.tsv" % args.output
df_mean.to_csv(path_or_buf=tsv, sep="\t")

fig = plt.figure()
colors = ["#d64b4b", "#d6934b", "#e3d46f", "#bce36f", "#75e36f", "#6fe3af", "#6fe3dd", "#6f96e3", "#3c30db", "#9964d1", "#cf64d1", "#d16495", "#736e70", "#000000", "#f58700", "#00f521"]
sns.set_palette(sns.color_palette(colors))
sns.set(style="whitegrid")
sns.lineplot(data=df_mean, linewidth=1.5, dashes=False)
fig.savefig("%s_seed2-8.mean.pdf" % args.output)
