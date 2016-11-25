#! /usr/bin/env python3

from Bio import SeqIO
import gzip

codon_dict = {}

with gzip.open("./T7_ref.fasta.gz", "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        i = 0
        while i <= (len(record.seq) - 6):
            codon_pair = str(record.seq[i:i+6])
            if codon_pair in codon_dict:
                codon_dict[codon_pair] += 1
            else:
                codon_dict[codon_pair] = 1
            i += 3

with open("T7_cpb.csv", "w") as outfile:
    outfile.write("codon_pair,count\n")
    for key, value in codon_dict.items():
        outfile.write(key + "," + str(value) + "\n")
