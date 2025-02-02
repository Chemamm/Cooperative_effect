#!/usr/bin/env/python3

import argparse
import Combinatory

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Orthology data file")
parser.add_argument("-n", "--number", type=int, help="Number of species to be considered as conserved")
parser.add_argument("-mir", "--mirna", help="mirna file")
parser.add_argument("-hsa", "--hsa", help="Hsa PosConsensus")
parser.add_argument("-mmu", "--mmu", help="Mmu PosConsensus")
parser.add_argument("-mmo", "--mmo", help="Mmo PosConsensus")
parser.add_argument("-cya", "--cya", help="Cya PosConsensus")
parser.add_argument("-ocu", "--ocu", help="Ocu PosConsensus")
parser.add_argument("-vvu", "--vvu", help="Vvu PosConsensus")
parser.add_argument("-s", "--stat", help="Stat output filename")
parser.add_argument("-o", "--out_label", help="output label")
args = parser.parse_args()


with open(args.mirna) as seeds:
    mirna_list = [x.replace("\n","") for x in list(seeds)]
Combinatory.target_conservation_per_mirna(args.file, mirna_list, args.hsa, args.mmu, args.mmo, args.cya, args.ocu, args.vvu, args.out_label)
Combinatory.target_conservation_stats(mirna_list, args.out_label, args.number, args.stat)
Combinatory.target_conservation_species(mirna_list, args.out_label, args.number)
