#!/usr/bin/env/python3

import pandas as pd
from mlxtend.frequent_patterns import apriori
from mlxtend.frequent_patterns import association_rules
import os
import argparse
import time
import random
import sys
import resource

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="pertranscript.txt file")
parser.add_argument("-c", "--csv", help="csv file with gene,miRNA,boolean")
parser.add_argument("-su", "--support", help="min_support")
parser.add_argument("-fi", "--filter", action="store_true", help="chose this to filter by length")
parser.add_argument("-s", "--size", type=int, help="size of the subset")
parser.add_argument("-g", "--gene", action="store_true", help="True to obtain the list of genes")
parser.add_argument("-l", "--lift", help="minimum lift to be significative")
parser.add_argument("-se", "--seed", help="miRNA fasta that miRNAconstarget has used")
parser.add_argument("-sp", "--specie", help="miRNA specie")
parser.add_argument("-rs", "--random_seed", help="Number of random seeds u want to chose")
parser.add_argument("-lm", "--low_memory", action="store_true", default=False, help="Low memory usage")
parser.add_argument("-la", "--label", help="label for output files (includes output folder if needed)")
parser.add_argument("-a", "--apriori", default=False, action="store_true", help="do the apriori calculation")
args = parser.parse_args()

# En este programa se realiza el algoritmo a priori sobre un perTranscript.txt obtenido a partir de miRNAconsTarget.
# Primero obtiene un set_tmp.csv que sirve para leerlo como dataframe, una vez leído como dataframe generamos una matriz
# de 1 y 0 donde 1 es el miRNA targetea el gen y 0 es no lo targetea. Una vez hecho esto se pasa el algoritmo apriori de mlxtend.
# Esto proporciona un grupo de subconjuntos de miRNA que superon el support de support.
# Además, he añadido opciones para filtrar los subconjuntos por longitud y para obtener la lista inversa, es decir, que genes se ven
# targeteados por estos subconjuntos.

# Obteniendo un diccionario a partir del pertranscript.txt
def pertranscript_dicter(input):
    pertranscript_dict = dict()
    for i in range(1, len(input)):
        line = input[i]
        pertranscript_dict[line.split("\t")[0]] = line.split("\t")[1].split(",")
    return pertranscript_dict


# Proceso para tomar únicamente los valores 1 u 0 en la matriz. Solo funciona como seguro en nuestro caso.
def encode_units(x):
    if x <= 0:
        return 0
    if x >= 1:
        return 1


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




# Generando el csv a partir del diccionario.

def apriori(args):
    if args.file:
        pertranscript = list(open(args.file))
        pertranscript_dict = pertranscript_dicter(pertranscript)
        tmp_name = "%s_tmp.csv" % (args.label)
        with open(tmp_name, "wt") as tmp_out:
            for item in pertranscript_dict.keys():
                for value in pertranscript_dict[item]:
                    tmp_out.write("%s,%s,1\n" % (item, value))

        # Modificando para seeds.
        if args.seed:
            fasta = list(open(args.seed))
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

            localtime = time.asctime(time.localtime(time.time()))
            print("Local current time :", localtime)
            if os.path.isfile("%s.seed.coding" % (args.seed)) == False and os.path.isfile("%s.seed" % (args.seed)) == False:
                # Escribiendo los ficheros de salida.
                with open("%s.seed" % args.seed, "wt") as out1:
                    for key in seed_dict:
                        out1.write(">%s\n%s\n" % (seed_dict[key], key))
                with open("%s.seed.coding" % (args.seed), "wt") as out2:
                    for key in coding_dict:
                        out2.write("%s\t%s\n" % (key, coding_dict[key]))
                # Modificando el csv
            with open(tmp_name) as tmp_file:
                tmp_str = tmp_file.read()
                # Para tomar un numero de semillas aleatorias de la muestra
                if args.random_seed:
                    random_sample = random.sample(list(coding_dict.values()), int(args.random_seed))
                    for key, seed in coding_dict.items():
                        if seed in random_sample:
                            tmp_str = tmp_str.replace(key, seed)
                else:
                    for key in coding_dict:
                        tmp_str = tmp_str.replace(key, coding_dict[key])
            tmp_name = "%s.parsed" % (tmp_name)
            with open(tmp_name, "wt") as out:
                out.write(tmp_str)
            if args.random_seed:
                delete_non_seed_cmd = "sed -i '/seed/!d' %s" % tmp_name
                os.system(delete_non_seed_cmd)
            localtime = time.asctime(time.localtime(time.time()))
            print("Local current time :", localtime)
            cmd = "sort -u -o %s_unique %s" % (tmp_name, tmp_name)
            os.system(cmd)
            localtime = time.asctime(time.localtime(time.time()))
            print("Local current time :", localtime)
            tmp_name = "%s_unique" % (tmp_name)

        os.system("sed -i '1i gene,miRNA,boolean' %s" % (tmp_name))

        # Obteniendo el dataframe.
        df = pd.read_csv(tmp_name)
        print(df.head())

    if args.csv:
        df = pd.read_csv(args.csv)
        print(df.head())

    if args.apriori:

        # generando la matriz necesaria para el algoritmo a priori
        basket = df.groupby(['gene', 'miRNA'])['boolean'].sum().unstack().reset_index().fillna(0).set_index("gene")

        localtime = time.asctime(time.localtime(time.time()))
        print("Local current time :", localtime)

        basket_sets = basket.applymap(encode_units)
        print(basket_sets)

        frequent_itemsets = apriori(basket_sets, min_support=float(args.support), use_colnames=True, low_memory=args.low_memory)

        localtime = time.asctime(time.localtime(time.time()))
        print("Local current time :", localtime)

        frequent_itemsets.to_csv("%sfrequent_items_supp_%s.txt" % (args.label, args.support), sep="\t")

        rules = association_rules(frequent_itemsets, metric="lift", min_threshold=float(args.lift))
        rules_ordered = rules.sort_values("lift", ascending=False)
        rules_ordered.to_csv("%srules_supp_%s_lift_%s.txt" % (args.label, args.support, args.lift), sep="\t")

        localtime = time.asctime(time.localtime(time.time()))
        print("Local current time :", localtime)

        # Formateando los ficheros

        cmd_freq1 = "sed -i 's/frozenset({//g' %sfrequent_items_supp_%s.txt" % (args.label, args.support)
        cmd_freq2 = "sed -i 's/})//g' %sfrequent_items_supp_%s.txt" % (args.label, args.support)
        cmd_rule1 = "sed -i 's/frozenset({//g' %srules_supp_%s_lift_%s.txt" % (args.label, args.support, args.lift)
        cmd_rule2 = "sed -i 's/})//g' %srules_supp_%s_lift_%s.txt" % (args.label, args.support, args.lift)

        os.system(cmd_freq1)
        os.system(cmd_freq2)
        os.system(cmd_rule1)
        os.system(cmd_rule2)

        # filtrando por longitud si procede
        if args.filter:
            frequent_itemsets['size'] = frequent_itemsets['itemsets'].apply(lambda x: len(x))
            frequent_itemsets = frequent_itemsets[(frequent_itemsets['size'] == args.size)]
            frequent_itemsets.to_csv("frequent_items_supp_%s.filtered_by_size_%s.txt" % (args.support, str(args.size)), sep="\t")

            rules["itemset_len"] = rules["antecedents"].apply(lambda x: len(x)) + rules["consequents"].apply(lambda x: len(x))
            rules_filtered = rules[(rules['itemset_len'] == args.size)]
            rules_filtered_ordered = rules_filtered.sort_values("lift", ascending=False)
            rules_filtered_ordered.to_csv("%srules_supp_%s_lift_%s.filtered_by_size_%s.txt" % (args.label, args.support, args.lift, str(args.size)), sep="\t")

            cmd_freq1 = "sed -i 's/frozenset({//g' %sfrequent_items_supp_%s.filtered_by_size_%s.txt" % (args.label, args.support, str(args.size))
            cmd_freq2 = "sed -i 's/})//g' %sfrequent_items_supp_%s.filtered_by_size_%s.txt" % (args.label, args.support, str(args.size))
            cmd_rule1 = "sed -i 's/frozenset({//g' %srules_supp_%s_lift_%s.filtered_by_size_%s.txt" % (args.label, args.support, args.lift, str(args.size))
            cmd_rule2 = "sed -i 's/})//g' %srules_supp_%s_lift_%s.filtered_by_size_%s.txt" % (args.label, args.support, args.lift, str(args.size))

            os.system(cmd_freq1)
            os.system(cmd_freq2)
            os.system(cmd_rule1)
            os.system(cmd_rule2)

        if args.gene:

            genes = []

            # Tomando los diez antecedentes-consecuentes con mayor lift
            if args.size:
                top10 = rules_filtered_ordered[:10]
            else:
                top10 = rules_ordered[:10]

            # Pasando los antecedentes y consecuentes a una lista, y uniéndolos.
            print(top10)
            antecedentes = top10["antecedents"].tolist()
            consecuentes = top10["consequents"].tolist()

            # Uniendo las dos listas
            subconj = []

            for i in range(0, len(antecedentes)):
                merge = list(antecedentes[i])
                merge.extend(list(consecuentes[i]))
                subconj.append(",".join(merge))

            # sacar la longitud de match entre listas con set(a) & set(b) y escribir ambas listas
            if args.size:
                out = open("genes_size_%s_supp_%s.txt" % (str(args.size), args.support), "wt")
            else:
                out = open("genes_supp_%s.txt" % (args.support), "wt")
            out.write("genes\tmiRNAs\n")
            for item in subconj:
                value = item.split(",")
                for key in pertranscript_dict.keys():
                    if len(set(value) & set(pertranscript_dict[key])) == len(value):
                        genes.append(key)
                genes_formated = ",".join(genes)
                out.write("%s\t%s\n" % (genes_formated, value))
            localtime = time.asctime(time.localtime(time.time()))
            print("Local current time :", localtime)

        # Removing tmp files
        if args.seed:
            #os.system("rm %s" %tmp_name)
            os.system("rm %s" % tmp_name.replace("_unique", ""))
            os.system("rm %s" % tmp_name.replace(".parsed_unique", ""))

    def main(args):
        apriori(args)

    def memory_limit():
        soft, hard = resource.getrlimit(resource.RLIMIT_AS)
        resource.setrlimit(resource.RLIMIT_AS, (get_memory() * 1024 / 4 * 3, hard))

    def get_memory():
        with open('/proc/meminfo', 'r') as mem:
            free_memory = 0
            for i in mem:
                sline = i.split()
                if str(sline[0]) in ('MemFree:', 'Buffers:', 'Cached:'):
                    free_memory += int(sline[1])
        return free_memory

    if __name__ == '__main__':
        memory_limit()  # Limitates maximun memory usage to 3/4
        try:
            main()
        except MemoryError:
            sys.stderr.write('\n\nERROR: Memory Exception\n')
            sys.exit(1)
