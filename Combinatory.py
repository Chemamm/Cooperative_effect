#!/usr/bin/env/python3

import matplotlib.pyplot as plt
import seaborn as sns
import operator
import pandas as pd


class Consensus:
    def __init__(self, file):
        with open(file) as cns:
            self.lines = list(cns)

    class Line:
        def __init__(self, line):
            self.row = line.split("\t")
            self.mirna = self.row[0]
            self.gene = self.row[1]
            self.nmethods = self.row[2]
            self.methods = self.row[3]

    def get_consensus_genes(self, input_list, mature_file, n, output):
        genes = []

        with open(input_list) as inp:
            inp_list = remove_linebreak_from_list(list(inp))
        inp_list = [x for x in inp_list if x]

        inp_list, inp_dict = get_complete_input_list(inp_list, mature_file)
        for line in self.lines:
            c_line = Consensus.Line(line)
            mirna = c_line.mirna
            gene = c_line.gene.split(":")[0]
            nmethods = int(c_line.nmethods)
            if mirna in inp_list and gene not in genes and nmethods >= n:
                genes.append(gene)

        with open(output, "wt") as out:
            out.write("\n".join(genes))

    def get_overlap_by_mirna(self, input_list, degradome_file, mature_file, n, n1, n2, output):
        out_gene = output.replace(".txt","_genes.txt")
        gene_dict = dict()
        degradome_dict = dict()
        overlap_dict = dict()
        overlap_genes = []
        with open(input_list) as inp:
            inp_list = remove_linebreak_from_list(list(inp))
        inp_list = [x for x in inp_list if x]

        inp_list, inp_dict = get_complete_input_list(inp_list, mature_file)

        dgf = Degradome(degradome_file)
        with open(output, "wt") as out:
            out.write("miRNAs (same sequence)\tNumber of genes\tGenes\tNumber of genes confirmed\tGenes confirmed\tNumber of degradome genes\tDegradome genes\n")
            for key in inp_dict:
                gene_dict[key] = []
                overlap_dict[key] = []
                degradome_dict[key] = []
                for line in self.lines:
                    c_line = Consensus.Line(line)
                    mirna = c_line.mirna
                    gene = c_line.gene.split(":")[0]
                    nmethods = int(c_line.nmethods)
                    if mirna in inp_dict[key] and gene not in gene_dict[key] and nmethods >= n:
                        gene_dict[key].append(gene)
                for j in range(1, len(dgf.lines)):
                    d_line = Degradome.Line(dgf.lines[j])
                    gene = d_line.gene.split(":")[0]
                    mirna = d_line.mirna
                    csite = int(d_line.cleavagesite)
                    if mirna in inp_dict[key] and gene not in degradome_dict[key] and n1 <= csite <= n2:
                        degradome_dict[key].append(gene)

                for gene in gene_dict[key]:
                    if gene in degradome_dict[key]:
                        overlap_dict[key].append(gene)
                        if gene not in overlap_genes:
                            overlap_genes.append(gene)
                n_genes = len(gene_dict[key])
                n_degradome = len(degradome_dict[key])
                n_overlap = len(overlap_dict[key])
                if not gene_dict[key]:
                    gene_dict[key] = ["NaN"]
                if not degradome_dict[key]:
                    degradome_dict[key] = ["NaN"]
                if not overlap_dict[key]:
                    overlap_dict[key] = ["NaN"]
                out.write("%s\t%i\t%s\t%i\t%s\t%i\t%s\n" % (",".join(inp_dict[key]), n_genes, ",".join(gene_dict[key]),
                                                   n_overlap, ",".join(overlap_dict[key]), n_degradome, ",".join(degradome_dict[key])))
        with open(out_gene,"wt") as out:
            for gene in overlap_genes:
                out.write("%s\n" %gene)

class Degradome:
    def __init__(self, file):
        with open(file) as dgd:
            self.lines = list(dgd)

    class Line:
        def __init__(self, line):
            self.row = line.split("\t")
            self.gene = self.row[0]
            self.mirna = self.row[4]
            self.cleavagesite = self.row[6]

    def get_genes(self, input_list, mature_file, n1, n2, output):
        genes = []

        with open(input_list) as inp:
            inp_list = remove_linebreak_from_list(list(inp))
        inp_list = [x for x in inp_list if x]

        inp_list, inp_dict = get_complete_input_list(inp_list, mature_file)

        for j in range(1, len(self.lines)):
            if not self.lines[j].startswith("target"):
                d_line = Degradome.Line(self.lines[j])
                gene = d_line.gene.split(":")[0].split(".")[0]
                mirna = d_line.mirna
                csite = int(d_line.cleavagesite)
                if mirna in inp_list and gene not in genes and n1 <= csite <= n2:
                    genes.append(gene)
        print(len(genes))
        with open(output, "wt") as out:
            out.write("\n".join(genes))

    def get_genes_bymiRNA(self, input_list, mature_file, n1, n2, out_folder):
        genes = dict()

        if input_list:
            with open(input_list) as inp:
                inp_list = remove_linebreak_from_list(list(inp))
            inp_list = [x for x in inp_list if x]


            inp_list, inp_dict = get_complete_input_list(inp_list, mature_file)


            for j in range(1, len(self.lines)):
                d_line = Degradome.Line(self.lines[j])
                gene = d_line.gene.split(":")[0].split(".")[0]
                mirna = d_line.mirna
                csite = int(d_line.cleavagesite)
                if mirna in inp_list and mirna not in genes.keys() and n1 <= csite <= n2:
                    genes[mirna] = [gene]
                elif mirna in inp_list and gene not in genes[mirna] and n1 <= csite <= n2:
                    genes[mirna].append(gene)
        else:
            for j in range(1, len(self.lines)):
                if not self.lines[j].startswith("target"):
                    d_line = Degradome.Line(self.lines[j])
                    gene = d_line.gene.split(":")[0].split(".")[0]
                    mirna = d_line.mirna.replace("*","-")
                    print(self.lines[j])
                    csite = int(d_line.cleavagesite)
                    print(csite)

                    if mirna not in genes.keys() and n1 <= csite <= n2:
                        genes[mirna] = [gene]
                    elif n1 <= csite <= n2 and gene not in genes[mirna]:
                        genes[mirna].append(gene)

        for key in genes:
            outname_degradome = "%s/%s_degradome_genes.txt" % (out_folder, key)
            with open(outname_degradome, "wt") as out:
                out.write("\n".join(genes[key]))

class Background:
    def __init__(self,file):
        with open(file) as bc:
            self.lines = list(bc)
        self.dict = {}
        for line in self.lines:
            row = line.replace("\n","").split("\t")
            self.dict[row[0]] = row[1].split(",")


    def removeRedundancy(self,output):
        non_redundant_dict = dict()
        for key in self.dict:
            non_redundant_key = key.split(".")[0]
            if non_redundant_key not in non_redundant_dict.keys():
                non_redundant_dict[non_redundant_key] = self.dict[key]
            else:
                non_redundant_dict[non_redundant_key] = union(non_redundant_dict[non_redundant_key], self.dict[key])

        with open(output, "wt") as out:
            for key in non_redundant_dict:
                out.write("%s\t%s\n" %(key, ",".join(non_redundant_dict[key])))

class Enrichment:
    def __init__(self, file):
        with open(file) as enr:
            self.lines = list(enr)

    class Line:
        def __init__(self, line):
            self.row = line.split("\t")
            self.GO = self.row[0]
            self.description = self.row[1]
            self.inbackground = int(self.row[2])
            self.notinbackground = int(self.row[3])
            self.percentageinbackground = self.inbackground / (self.inbackground + self.notinbackground) * 100
            self.ingenelist = int(self.row[4])
            self.notingenelist = int(self.row[5])
            self.percentageingenelist = self.ingenelist / (self.ingenelist + self.notingenelist) * 100
            self.pvalue = float(self.row[6])
            self.fdr = float(self.row[7])
            self.oddsratio = float(self.row[8])
            self.RE = float(self.row[9])
            self.type = self.row[10]


    def depurate(self, GO_description_file, output):
        go_descriptions = GOobo(GO_description_file)

        with open(output, "wt") as out:
            out.write(self.lines[0])
            for j in range(1, len(self.lines)):
                en_line = Enrichment.Line(self.lines[j])
                if int(en_line.ingenelist) > 0:
                    en_line.row[1] = go_descriptions.dict[en_line.GO]
                    out.write("\t".join(en_line.row))

    def depurateKEG(self, KEG_description_file, output):
        go_descriptions = KEGobo(KEG_description_file)

        with open(output, "wt") as out:
            out.write(self.lines[0])
            for j in range(1, len(self.lines)):
                en_line = Enrichment.Line(self.lines[j])
                if int(en_line.ingenelist) > 0:
                    en_line.row[1] = go_descriptions.dict[en_line.GO]
                    out.write("\t".join(en_line.row))

    def barplot_rawdata(self,outname, number):
        bar_dict = dict()
        for j in range(1, len(self.lines)):
            en_line = Enrichment.Line(self.lines[j])
            if float(en_line.pvalue) < 0.05:
                bar_dict[en_line.description.split(",")[0]] = int(en_line.ingenelist)
        bar_dict = order_dictionary(bar_dict)
        barhplot_from_dict(bar_dict,outname,"Absolute value", number)

    def barplot_proportion(self,outname, number, filter_n=False):
        bar_dict = dict()
        for j in range(1, len(self.lines)):
            en_line = Enrichment.Line(self.lines[j])
            if float(en_line.pvalue) < 0.05:
                if filter_n and float(en_line.ingenelist) >= filter_n:
                    bar_dict[en_line.description.split(",")[0]]=float(en_line.ingenelist)/float(en_line.inbackground) * 100
                elif filter_n==False:
                    bar_dict[en_line.description.split(",")[0]] = float(en_line.ingenelist) / float(
                        en_line.inbackground) * 100
        bar_dict = order_dictionary(bar_dict)
        barhplot_from_dict(bar_dict, outname,"inGeneList/inBackground (%)", number)

    def barplot_oddsratio(self, outname, number):
        bar_dict = dict()
        for j in range(1, len(self.lines)):
            en_line = Enrichment.Line(self.lines[j])
            if float(en_line.pvalue) < 0.05:
                bar_dict[en_line.description.split(",")[0]] = float(en_line.oddsratio)
        bar_dict = order_dictionary(bar_dict)
        barhplot_from_dict(bar_dict, outname, "Odds Ratio", number)

    def barplot_bujun(self, out_label, number):
        bar_depleted_dict = dict()
        bar_enriched_dict = dict()
        for j in range(1, len(self.lines)):
            en_line = Enrichment.Line(self.lines[j])
            if en_line.pvalue < 0.05 and en_line.type == "depleted":
                bar_depleted_dict[en_line.description.split(",")[0]] = (en_line.percentageingenelist,
                                                                        en_line.percentageinbackground)
            elif en_line.pvalue < 0.05 and en_line.type == "enriched":
                bar_enriched_dict[en_line.description.split(",")[0]] = (en_line.percentageingenelist,
                                                                        en_line.percentageinbackground)
        if len(bar_depleted_dict) > 0:
            barhplot_bujun(bar_depleted_dict, "%s_depleted_%i.png" %(out_label, number), number)
        if len(bar_enriched_dict) > 0:
            barhplot_bujun(bar_enriched_dict, "%s_enriched_%i.png" % (out_label, number), number)
            print(bar_enriched_dict)

class Genelist:
    def __init__(self, file):
        with open(file) as gen:
            self.lines = list(gen)

        self.genes = []
        for line in self.lines:
            self.genes.append(line.replace("\n", ""))

    def getGO(self, file, out):
        with open(file) as genGO:
            genGO = list(genGO)
        with open(out, "wt") as out:
            for line in genGO:
                gene = line.split("\t")
                if gene in self.genes:
                    out.write(line)

class GOobo:
    def __init__(self, file):
        with open(file) as obo:
            self.lines = list(obo)
        self.dict = {}
        self.namespace_dict = {}
        for n in range(0, len(self.lines)):
            if "[Term]" in self.lines[n]:
                ID = self.lines[n + 1].split(": ")[1].replace("\n", "")
                name = self.lines[n + 2].split(": ")[1].replace("\n", "")
                namespace = self.lines[n + 3].split(": ")[1].replace("\n", "")
                self.dict[ID] = name
                self.namespace_dict[ID] = namespace
            elif "alt_id:" in self.lines[n]:
                ID = self.lines[n].split(": ")[1].replace("\n", "")
                self.dict[ID] = name
                self.namespace_dict[ID] = namespace

class KEGobo:
    def __init__(self,file):
        with open(file) as ft:
            self.lines = list(ft)
        self.dict = dict()
        for line in self.lines:
            ID = line.split("\t")[0].split(":")[1]
            description = line.split("\t")[1].replace("\n","")
            self.dict[ID] = description

class Fasta:
    def __init__(self, file):
        with open(file) as ft:
            self.lines = list(ft)
        self.dict = self.fasta2dict()

    def split(self, n, output):
        file_number = 1
        count = 0
        out = open("%s_%s.fa" % (output, file_number), "wt")
        for key in self.dict:
            out.write(">%s\n%s\n" % (key, self.dict[key]))
            count += 1
            if count >= int(n):
                out.close()
                file_number += 1
                count = 0
                out = open("%s_%s.fa" % (output, file_number), "wt")

    def fasta2dict(self):
        fasta_dict = dict()
        seqid = ""
        for i in range(0, len(self.lines)):
            if self.lines[i].startswith('>'):
                if seqid != "":
                    fasta_dict[seqid] = seq
                seqid = self.lines[i].split(' ')[0].replace(">", "").replace("\n", "")
                seq = ""
            else:
                seq += self.lines[i].replace("\n", "")
        if seqid != "":
            fasta_dict[seqid] = seq
        return fasta_dict

class Overlap:
    def __init__(self, file):
        with open(file) as ft:
            self.lines = list(ft)

    class Line:
        def __init__(self, line):
            line = line.replace("\n","").split("\t")
            self.transcript = line[0]
            self.n_consensus_genes = line[1]
            self.consensus_genes = line[2]
            self.n_overlap_genes = line[3]
            self.overlap_genes = line[4]
            self.n_degradome_genes = line[5]
            self.degradome_genes = line[6]



    def get_gene_lists(self,out_folder):
        for j in range(1, len(self.lines)):
            line = Overlap.Line(self.lines[j])
            outname_consensus = "%s/%s_consensus_genes.txt" %(out_folder, line.transcript.split(",")[0])
            outname_degradome = "%s/%s_degradome_genes.txt" % (out_folder, line.transcript.split(",")[0])
            outname_overlap = "%s/%s_overlap_genes.txt" % (out_folder, line.transcript.split(",")[0])
            if line.consensus_genes != "NaN":
                with open(outname_consensus,"wt")as out:
                    out.write("\n".join(line.consensus_genes.split(",")))
            if line.degradome_genes != "NaN":
                with open(outname_degradome, "wt")as out:
                    out.write("\n".join(line.degradome_genes.split(",")))
            if line.overlap_genes != "NaN":
                with open(outname_overlap, "wt")as out:
                    out.write("\n".join(line.overlap_genes.split(",")))

class Ematrix:
    def __init__(self,file):
        with open(file) as em:
            self.lines = list(em)
        self.header = self.lines[0]
        self.conditions = self.header.replace("+","_POS").replace("-","_NEG").split("\t")[1::]

    def getBujunTableSuperExactTest(self, out_label):
        dict1 = dict()
        dict2 = dict()
        for condition in self.conditions:
            dict1[condition] = []
            dict2[condition] = []
        for j in range(1, len(self.lines)):
            line = self.lines[j].split("\t")
            ID = line[0]
            for k in range(1,len(line)):
                if float(line[k]) > 0:
                    dict1[self.conditions[k-1]].append(ID)
                if float(line[k]) >= 100:
                    dict2[self.conditions[k-1]].append(ID)
        df1 = pd.DataFrame.from_dict(dict1, orient="index").T
        df2 = pd.DataFrame.from_dict(dict2, orient = "index").T
        df1.to_csv("%s_0.tsv" % out_label, index=False, sep="\t")
        df2.to_csv("%s_100.tsv" %out_label, index=False, sep="\t")

        union_1=[]
        union_2=[]
        for key in dict1:
            union_1 = union(union_1, dict1[key])
            union_2 = union(union_2,dict2[key])
        print("Background for 0: %i" %len(union_1))
        print("Background for 100: %i" % len(union_2))

class PositionalConsensus:
    def __init__(self,file):
        with open(file) as fl:
            self.lines = list(fl)

    class Line:
        def __init__(self,line):
            self.row = line.split("\t")
            self.mirna = self.row[0]
            self.gene = self.row[1]
            self.start = int(self.row[2])
            self.end = int(self.row[3])

    def get_num_targets(self, list_mirna, output):
        with open(list_mirna) as lst:
            lst = [x.replace("\n","") for x in list(lst)]
        target_dict = dict()
        num_target_dict = dict()
        for k in range(1,len(self.lines)):
            line = PositionalConsensus.Line(self.lines[k])
            if line.gene in target_dict:
                target_dict[line.gene].append(line.mirna)
            else:
                target_dict[line.gene] = [line.mirna]
        print(target_dict)
        for gene in target_dict:
            if not gene.startswith("BRAF_"):
                num_target_dict[gene] = []
                for mirna in lst:
                    if mirna in target_dict[gene]:
                        num_target_dict[gene].append(target_dict[gene].count(mirna))
                    else:
                        num_target_dict[gene].append(0)
        df = pd.DataFrame.from_dict(num_target_dict)
        df["miRNA"] = lst
        df = df.set_index("miRNA")
        df["Total targets"] = df.sum(axis=1)
        df.to_csv(output, sep="\t")










def get_complete_input_list(input_list, mature_file):
    complete_input_list = []
    complete_dict = dict()
    mf = Fasta(mature_file)
    mf_dict = mf.fasta2dict()
    for mirna in input_list:
        complete_dict[mirna]= []
        if mirna in mf_dict.keys():
            for mirna_mature in mf_dict:
                if mf_dict[mirna] == mf_dict[mirna_mature]:
                    complete_input_list.append(mirna_mature)
                    complete_dict[mirna].append(mirna_mature)
        else:
            complete_input_list.append(mirna)
            complete_dict[mirna].append(mirna)


    return complete_input_list, complete_dict

def remove_linebreak_from_list(list):
    list_clean = []
    for item in list:
        list_clean.append(item.replace("\n", ""))
    return list_clean

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def union(lst1, lst2):
    final_list = list(set(lst1) | set(lst2))
    return final_list

def overlap_plot(file, output_name):
    output_name_2 = output_name.split(".")[0] + "_percentage.png"
    with open(file) as list_of_files:
        list_of_files = list(list_of_files)
    overlap_dict = dict()
    overlap_perc_dict = dict()
    for file in list_of_files:
        filename = file.replace("\n","")
        filename2 = filename.replace("consensus", "KUK10_EXC10")
        name = filename.split("/")[-1].replace("_consensus_genes.txt","")
        outname = name + "_overlap.txt"
        with open(filename) as f1:
            f1 = list(f1)
        with open(filename2) as f2:
            f2=list(f2)
        overlap_list = intersection(f1,f2)
        overlap_dict[name] = len(overlap_list)
        overlap_perc_dict[name] = len(overlap_list) / len(f1) * 100
        with open(outname,"wt") as out:
            for item in overlap_list:
                out.write(item)
    barplot_from_dict(overlap_dict,output_name, "Number of genes overlapping")
    barplot_from_dict(overlap_perc_dict, output_name_2, "Percentage of genes overlapping")

def get_overlap_genes(file1,file2, out):
    f1 = Genelist(file1)
    f2 = Genelist(file2)
    overlap = intersection(f1.genes, f2.genes)
    for gene in overlap:
        out.write("%s\n" %gene)

def order_dictionary(D):
    D_ord = {}
    key_ord = sorted(D, key=D.get)
    for key in reversed(key_ord):
        D_ord[key] = D[key]
    return D_ord

def filter_dictionary(D, n):
    D_filt = {}
    for key in D:
        if float(D[key]) >= n:
            D_filt[key] = D[key]
    return D_filt

def piechart(labels, sizes, out):
    # Pie chart, where the slices will be ordered and plotted counter-clockwise:
    # labels = 'Frogs', 'Hogs', 'Dogs', 'Logs'
    # sizes = [15, 30, 45, 10]
    # explode = (0, 0.1, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')
    sns.set_theme(style="whitegrid")

    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
            shadow=True, startangle=90)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    fig1.savefig(out)

def barplot_from_dict(D,out,label):
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(20, 20))
    plt.bar(range(len(D)), list(D.values()), align='center')
    plt.xticks(range(len(D)), list(D.keys()), rotation='vertical')
    plt.ylabel(label)
    plt.tight_layout()
    #plt.subplots_adjust(bottom=0.3)
    plt.savefig(out)
    plt.close()

def barhplot_from_dict(D,out,label,number, filter_n=False):
    values = list(D.values())
    values.reverse()
    keys = list(D.keys())
    keys.reverse()
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(20, 20))
    plt.barh(range(len(D)), values, align='center')
    plt.yticks(range(len(D)), keys)
    plt.xlabel(label)
    plt.rc("ytick", labelsize=10)
    plt.tight_layout()
    #plt.subplots_adjust(bottom=0.3)
    plt.savefig(out)
    plt.close()
    if len(D) > number:
        values = list(D.values())[0:number]
        values.reverse()
        keys = list(D.keys())[0:number]
        keys.reverse()
        sns.set_theme(style="whitegrid")
        plt.figure(figsize=(20, 10))
        plt.barh(range(number), values, align='center')
        plt.yticks(range(number), keys)
        plt.xlabel(label)
        plt.tight_layout()
        # plt.subplots_adjust(bottom=0.3)
        plt.savefig("%s_%i.png" %(out.replace(".png",""), number))
        plt.close()
    else:
        sns.set_theme(style="whitegrid")
        plt.figure(figsize=(20, 10))
        plt.barh(range(len(D)), values, align='center')
        plt.yticks(range(len(D)), keys)
        plt.xlabel(label)
        plt.tight_layout()
        # plt.subplots_adjust(bottom=0.3)
        plt.savefig("%s_%i.png" %(out.replace(".png",""), number))
        plt.close()

def barhplot_bujun(D,out,number):
    transcript_list = []
    ingene_list = []
    background_list = []

    if len(D) > number:
        n = 0
        for key in D:
            transcript_list.append(key)
            ingene_list.append(D[key][0])
            background_list.append(D[key][1])
            n += 1
            if n == int(number):
                break
    else:
        for key in D:
            transcript_list.append(key)
            ingene_list.append(D[key][0])
            background_list.append(D[key][1])

    df = pd.DataFrame({'% in gene list': ingene_list, '% in background': background_list}, index=transcript_list)
    df = df.reindex(index=df.index[::-1])
    sns.set_theme()
    ax = df.plot.barh(fontsize=16)
    ax.set_ylabel("GO term")
    fig = ax.get_figure()
    fig.set_figheight(10)
    fig.set_figwidth(18)
    plt.tight_layout()
    fig.savefig(out)
    plt.close()

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


def shared_seed(mirna, mirna_set):
    fasta = list(open(mirna_set))
    mirna_dict = fasta_to_dict(fasta)
    problem_seed = mirna[1:8]
    sharing_seed = []
    for key in mirna_dict:
        seed = mirna_dict[key][1:8]
        if seed == problem_seed:
            sharing_seed.append(key)

    if sharing_seed:
        print("The miRNA shares its seed with %i miRNAs from the host:\n%s" %(len(sharing_seed), ",".join(sharing_seed)))
    else:
        print("Specific seed")



