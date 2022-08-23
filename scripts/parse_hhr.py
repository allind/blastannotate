#! /usr/bin/env python
import sys
def main(argv):

	toprint = ""
	print("Seq\tAnnot\tNumprot\tTax\tTaxID\tRepID\tProbab\tEval\tScore\tAligned_cols\tIdent\tSimilarity\tSum_probs\tTemplate_Neff")
	for line in open(sys.argv[1]):

		line = line.strip('\n')
		if line == "--":
			print(toprint)
			toprint = ""
		else:
			if ">" in line:
				seq = line.split(' ')[0].strip('>')
				descrip = " ".join(line.split(' ')[1:]).split('n=')[0].strip(' ')
				toprint += seq + '\t' + descrip + '\t'
				

				num=line.split('n=')[1].split('Tax=')[0].strip(' ')
				tax = line.split('Tax=')[1].split('TaxID=')[0].strip(' ')
				taxid = line.split("TaxID=")[1].split('RepID=')[0].strip(' ')
				repid = line.split('RepID=')[1].strip(' ')
				toprint += num + '\t' + tax + '\t' + taxid + '\t' + repid + '\t'
			else:
				#no extra spaces here
				vals = [e.split('=')[1].strip(' ') for e in line.split(' ') if "=" in e]

				toprint += '\t'.join(vals)

	
	print(toprint)



if __name__ == "__main__":
  main(sys.argv)
