#! /usr/bin/env python
from Bio import SeqIO
import sys
def main(argv):
	
	genes = {}
	for seq in SeqIO.parse(sys.argv[1], 'fasta'):
		gene_id = seq.description.split(' ')[1]
		
		has_stop = True

		if "*" not in str(seq.seq):
			has_stop = False
		
		if gene_id not in genes:
			genes[gene_id] = []
		genes[gene_id].append(has_stop)
	

	for gene in genes:
		if False in genes[gene]:
			count_stop = genes[gene].count(True)
			count_nostop = genes[gene].count(False)
			print(gene + '\t' + str(count_nostop) + '\t' + str(count_stop)  + '\t' + str(count_nostop / len(genes[gene])))
if __name__ == "__main__":
  main(sys.argv)
