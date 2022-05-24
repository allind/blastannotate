#! /usr/bin/env python
#usage: script.py [input] [output gff] [output log]
import sys
def main(argv):

	removed = []
	towrite = []
	ok = False
	for line in open(sys.argv[1]):
		if "#" in line[0]:
			towrite.append(line + '\n')
		else:
			feat = line.split('\t')[2]
			if feat == "transcript":
				length = int(line.split('\t')[4]) - int(line.split('\t')[3])
				if length < 100:
					ok = False
					removed.append(line.split('\t')[-1].split(';')[0].split(' ')[1].strip('"'))
				else:
					ok = True
					towrite.append(line + '\n')
			else:
				if ok:
					towrite.append(line + '\n')
	
	dest = open(sys.argv[2], 'w')
	for d in towrite:
		dest.write(d)
	dest.close()
	
	log = open(sys.argv[3], 'w')
	for r in removed:
		log.write("Removed: " + r + '\n')
	log.close()
	
		
	
if __name__ == "__main__":
  main(sys.argv)
