#!/usr/bin/env/python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="multipleTargets.txt file")
parser.add_argument("-rs", "--random_seed", action="store_true", help="random seed type of analysis")
parser.add_argument("-af", "--apriori_file", help="unique file from ApriorimiRNA.py")
parser.add_argument("-m", "--miRNA", help="miRNA file")
args = parser.parse_args()

targets = 0
if args.random_seed:
    with open(args.apriori_file) as file:
        file = list(file)
        seeds = []
        for n in range(1, len(file)):
            seed = file[n].split(",")[1]
            if seed not in seeds:
                seeds.append(seed)
    with open("%s.seed.coding" % args.miRNA) as code:
        code = list(code)
        code_dict = {}
        for line in code:
            seed = line.split("\t")[1].replace("\n", "")
            mirna = line.split("\t")[0]
            if seed not in code_dict.keys() and seed in seeds:
                code_dict[seed] = [mirna]
            elif seed in code_dict.keys() and seed in seeds:
                code_dict[seed].append(mirna)
    with open(args.file) as ft:
        ft = list(ft)
        for n in range(1, len(ft)):
            mirna = ft[n].split("#")[0]
            if mirna in [x for v in code_dict.values() for x in v]:
                target = ft[n].split("\t")[2].replace("\n", "")
                targets += int(target)
    print(targets)


else:
    with open(args.file) as ft:
        ft = list(ft)
        for n in range(1, len(ft)):
            target = ft[n].split("\t")[2].replace("\n", "")
            targets += int(target)

    print(targets)
