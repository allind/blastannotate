#!/usr/bin/env python

import sys
import gzip
from Bio import SeqIO

def calculate_n50(sequence_lengths):
    sequence_lengths.sort(reverse=True)
    total_length = sum(sequence_lengths)
    target_length = total_length / 2

    current_sum = 0
    for length in sequence_lengths:
        current_sum += length
        if current_sum >= target_length:
            return length

def main(filename):
    sequence_lengths = []

    # Determine file format and open the file accordingly
    if filename.endswith(".gz"):
        handle = gzip.open(filename, "rt")
    else:
        handle = open(filename, "r")

    # Parsing the FASTA/FASTQ file and extracting sequence lengths
    with handle as handle:
        for record in SeqIO.parse(handle, "fasta" if (filename.endswith(".fasta") or filename.endswith(".fa")) else "fastq"):
            sequence_lengths.append(len(record.seq))

    if not sequence_lengths:
        print("No sequences found in the file.")
        return

    n50 = calculate_n50(sequence_lengths)
    print("N50 value:", n50)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_file>")
        sys.exit(1)
    filename = sys.argv[1]
    main(filename)
