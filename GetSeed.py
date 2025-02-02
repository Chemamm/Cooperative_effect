#!/usr/bin/env/python3

import argparse
import time

def fasta_to_dict(ft):
    fasta_dict = dict()
    seqid = ""
    for i in range(0, len(ft)):
        if ft[i].startswith('>'):
            if seqid != "":
                fasta_dict[seqid] = seq
            seqid = ft[i].split(' ')[0].replace(">","").replace("\n","")
            seq = ""
        else:
            seq += ft[i].replace("\n","")
    if seqid != "":
        fasta_dict[seqid] = seq
    return fasta_dict




parser = argparse.ArgumentParser()
parser.add_argument("-f","--file", help="miRNA file")
parser.add_argument("-s","--specie", help="miRNA specie")
args = parser.parse_args()

fasta = list(open(args.file))
mirna_dict = fasta_to_dict(fasta)

seed_dict = dict()
coding_dict = dict()
count = 1

## Creando dos diccionarios, uno que contenga la lista de semillas, y otro que contenga que miRNAs se corresponden con cada semilla.

for key in mirna_dict:
    seed = mirna_dict[key][1:8]
    if seed not in seed_dict.keys():
        seed_dict[seed] = "%s_seed_%s" % (args.specie, str(count))
        coding_dict[key] = seed_dict[seed]
        count += 1
    else:
        coding_dict[key] = seed_dict[seed]

## Escribiendo los ficheros de salida.
with open("%s.seed" % args.file, "wt") as out1:
    for key in seed_dict:
        out1.write(">%s\n%s\n" %(seed_dict[key], key))
with open("%s.seed.coding" % (args.file), "wt") as out2:
    for key in coding_dict:
        out2.write("%s\t%s\n" %(key, coding_dict[key]))

print(args.file,", ",count-1)
