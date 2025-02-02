import requests
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Fasta file with the identifiers")
parser.add_argument("-l", "--length", type=int, help="length of the flanking region" )
parser.add_argument("-o", "--output", help="Output")
args = parser.parse_args()


class Fasta:
    def __init__(self, file):
        with open(file) as ft:
            self.lines = list(ft)
        self.dict = self.fasta2dict()



    def fasta2dict(self):
        fasta_dict = dict()
        seqid = ""
        for i in range(0, len(self.lines)):
            if self.lines[i].startswith('>'):
                if seqid != "":
                    fasta_dict[seqid] = seq
                seqid = self.lines[i].split(' ')[0].replace(">", "").replace("\n", "")
                if "|" in seqid:
                    seqid = seqid.split("|")[0]
                if "." in seqid:
                    seqid = seqid.split(".")[0]
                seq = ""
            else:
                seq += self.lines[i].replace("\n", "")
        if seqid != "":
            fasta_dict[seqid] = seq
        return fasta_dict

def retrieve_flanking_sequence(ensembl_gene_id, flanking, n=1, length=1):
    url = f"https://rest.ensembl.org/sequence/id/{ensembl_gene_id}?expand_5prime=0;expand_3prime=%i" %flanking

    response = requests.get(url, headers={"Content-Type": "application/json"})

    if response.ok:
        downstream_sequence = response.json()['seq']
        print(f'Sequence of {ensembl_gene_id} retrieved {n}/{length}')
    else:
        print(f"Failed to retrieve downstream flanking sequence. Status code: {response.status_code}")

    if downstream_sequence:
        downstream_sequence = downstream_sequence[-flanking::]
    else:
        downstream_sequence = "Sequence not available"

    return downstream_sequence

def main():
    fa = Fasta(args.file)
    old_fa = False
    if os.path.exists(args.output):
        print("Previous file detected, looking for sequences already retrieved")
        old_fa = Fasta(args.output)

    with open(args.output, "a") as out:
        n = 1
        for id in fa.dict:
            if old_fa and id in old_fa.dict:
                downstream_seq = old_fa[id]
                n += 1
                out.write(f'>{id}\n{downstream_seq}\n')
            else:
                downstream_seq = retrieve_flanking_sequence(id, args.length, n, len(fa.dict))
                n +=1
                out.write(f'>{id}\n{downstream_seq}\n')


if __name__ == "__main__":
    main()

