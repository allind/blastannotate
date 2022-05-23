#! /usr/bin/env python
import sys
def main(argv):
	#minunique = int(sys.argv[2])
	minreads = int(sys.argv[2])
	for line in open(sys.argv[1]):
		line = line.strip('\n')
		unique = int(line.split('\t')[6])
		multi = int(line.split('\t')[7])
		total = unique + multi
		#if unique >= minunique:
		if total >= minreads:
			print(line)


if __name__ == "__main__":
  main(sys.argv)
