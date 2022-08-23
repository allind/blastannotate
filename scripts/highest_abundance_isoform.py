#! /usr/bin/env python
import sys,re
def main(argv):

	genes = {}
	for line in open(sys.argv[1]):
		line = line.strip('\n')
		if "target_id" not in line:
			#gene = ".".join(line.split('\t')[0].split('.')[0:2])
			gene = re.sub("\.\d+\w*","", line.split('\t')[0])
			#gene = re.sub("\.\d+","",line.split('\t')[0])
			transcript = line.split('\t')[0]
			tpm = float(line.split('\t')[-1])
			if gene not in genes:
				genes[gene] = [[],[]]
			genes[gene][0].append(transcript)
			genes[gene][1].append(tpm)
	for gene in genes:
		#if len(genes[gene][1]) > 1:
		index = genes[gene][1].index(max(genes[gene][1]))
		maxtrans = genes[gene][0][index]
		print(gene + '\t' + maxtrans)
if __name__ == "__main__":
  main(sys.argv)
