#! /usr/bin/env python
import sys
def main(argv):
	

	nc = [line.strip('\n') for line in open(sys.argv[1])]
	
	for line in open(sys.argv[2]):
		line = line.strip('\n')
		if line[0] == "#":
			print(line)
		else:
			fid = '.'.join(line.split('\t')[-1].split('ID=')[1].split(';')[0].split('.')[0:2])
			if fid not in nc:
				print(line)
		
if __name__ == "__main__":
  main(sys.argv)
