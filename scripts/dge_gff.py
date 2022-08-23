#! /usr/bin/env python
from Bio import SeqIO
import sys
def main(argv):

	dge = 'TGTTTGTT'
	dge_rc = 'AACAAACA'
	counter = 0
	for seq in SeqIO.parse(sys.argv[1], 'fasta'):
		fullseq = str(seq.seq)
		seqid = seq.id
		for i in range(0, len(fullseq) - 9):
			curr = fullseq[i: i+8]
			if curr == dge:
				counter += 1
				print(seqid + '\tdge_finder\tDGE\t' + str(i+1) + '\t' + str(i + 8) + '\t.\t+\t.\tID=DGE_' + str(counter) + ";")
			elif curr == dge_rc:
				counter += 1
				print(seqid + '\tdge_finder\tDGE\t' + str(i+1) + '\t' + str(i + 8) + '\t.\t-\t.\tID=DGE_' + str(counter) + ";")

		
if __name__ == "__main__":
  main(sys.argv)
