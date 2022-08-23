#! /usr/bin/env python
from Bio import SeqIO
import sys,re
def main(argv):

	highest_iso = {line.split('\t')[0]: line.strip('\n').split('\t')[1] for line in open(sys.argv[1])}
	for seq in SeqIO.parse(sys.argv[2], 'fasta'):
		gene = seq.description.split(' ')[1]
		iso = seq.id
		if gene in highest_iso:
			if iso == highest_iso[gene]:
				print(">" + seq.description)
				print(str(seq.seq))
		else:
			print(">" + seq.description)
			print(str(seq.seq))
if __name__ == "__main__":
  main(sys.argv)
