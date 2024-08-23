#!/usr/bin/env python

import sys
from Bio import SeqIO


if __name__ == "__main__":
	if len(sys.argv) != 2:
		print("Usage: python script.py <input_fasta>")
		sys.exit(1)

	input_file = sys.argv[1]
	length=0	
	for seq in SeqIO.parse(input_file, 'fasta'):
		length += len(str(seq.seq))
	print(length)
