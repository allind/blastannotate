#!/usr/bin/env python

import sys
from Bio import SeqIO

def parse_fasta(input_file):
    for record in SeqIO.parse(input_file, "fasta"):
        seq_id = record.id
        seq_length = len(record.seq)
        print(seq_id + "\t" + str(seq_length))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_fasta>")
        sys.exit(1)

    input_file = sys.argv[1]

    parse_fasta(input_file)
