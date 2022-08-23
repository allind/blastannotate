#! /usr/bin/env python
#input: script.py depthfile_of_zeros gff gff_output log_output
from Bio import SeqIO
import pprint
import sys

def main(argv):

	logwrite = []

	uncovered = {}


	for line in open(sys.argv[1]):
		line = line.strip('\n')
		contig = line.split('\t')[0]
		pos = line.split('\t')[1]
		if contig not in uncovered:
			uncovered[contig] = []
		uncovered[contig].append(pos)
	
	
	
	
	#parse stringtie gtf
	#dict: transcripts: {transcript: [start, stop, chrom, defline], transcript: [start,stop, chrom, defline]}
	#transcript_exons: {transcript: [exon, exon, exon], transcript: [exon, exon]}
	#exons: {exon: [start, stop, defline]}
	transcripts = {}
	exons = {}
	transcript_exons = {}

	transcripts_to_remove = []

	for line in open(sys.argv[2]):
		line = line.strip('\n')
		if "#" not in line[0]:
			start = line.split('\t')[3]
			stop = line.split('\t')[4]
			defline = line
			ident = line.split('\t')[2]
			chrom = line.split('\t')[0]

			if ident == "transcript":

				tid = line.split('\t')[-1].split("transcript_id")[1].split(';')[0].strip(' ').strip('"')
				if tid not in transcripts:
					transcripts[tid] = []
				currexons = []
				transcripts[tid] = [start, stop, chrom, defline]
			else:
				#is exon
				exonid = ';'.join(line.split('\t')[-1].split(';')[1:3])
				parent = line.split('\t')[-1].split("transcript_id")[1].split(';')[0].strip(' ').strip('"')
				exons[exonid] = [start, stop, defline]
				if parent not in transcript_exons:
					transcript_exons[parent] = []
				transcript_exons[parent].append(exonid)


	transcript_list = list(transcripts.keys())
	transcript_cuts = {} #transcript: # of cuts

	for t in transcript_list:
		transcript_cuts[t] = 0
		tstart = transcripts[t][0]
		tstop = transcripts[t][1]
		chrom = transcripts[t][2]
		oriname = t

		#iterate over exons to determine which need to be cut

		exons_to_cut = {}
		#dont examine first or last exon NEW CHANGE ACTUALLY DO
		#exons_to_examine = transcript_exons[t][1:len(transcript_exons[t])-2]
		exons_to_examine = transcript_exons[t]
		for exon in exons_to_examine:
			estart = exons[exon][0]
			estop = exons[exon][1]

			span = [str(e) for e in range(int(estart), int(estop))]


			zeros = []
			for base in span:
				if chrom in uncovered:
					if base in uncovered[chrom]:
						zeros.append(base)


			#collapse into runs. Discard runs of 1
			run = []
			groups = [run]
			expect = None 
			for b in zeros:
				if (int(b) == expect) or (expect is None):
					run.append(b)
				else:
					run = [b]
					groups.append(run)
				expect = int(b) + 1

			#remove single base drops and starts & stops w/in 5 bases of start/stop
			if groups != [[]]:
				for g in groups:
					removed = False
					if len(g) == 1:
						groups.remove(g)
						removed = True
					if (int(g[0]) - 6) < int(transcripts[t][0]) and removed == False:
						groups.remove(g)
						removed = True
					if (int(g[-1]) + 6) > int(transcripts[t][1]) and removed == False:
						groups.remove(g)
						removed = True

			#add the transcript cuts

			if len(groups) > 0 and groups != [[]]:
				#exons_to_cut[exon] = [zeros[0], zeros[-1], len(zeros)]
				#exons_to_cut[exon] =
				exons_to_cut[exon] = [] 
				#for g in groups:
				for i in range(0, len(groups)):
					transcript_cuts[t] += 1
					gstart = groups[i][0]
					gend = groups[i][-1]
					glen = len(groups[i])

					exons_to_cut[exon].append([gstart, gend, glen])
					#exons_to_cut structure: exon: [[gstart, gend, glen], [gstart, gend, glen]]

		#add log lines

		exons_to_cut_ids = list(exons_to_cut.keys())
		#sizes = [exons_to_cut[e][2] for e in exons_to_cut]
		if len(exons_to_cut_ids) > 0:
			logwrite.append("Below cutoff bases in transcript " + t + ". Exons affected total: " + str(len(exons_to_cut_ids)) + " exons. In exon numbers: " + ','.join(exons_to_cut_ids) + '\n')

		#if len(exons_to_cut_ids) > 0:
		#transcript_cuts = transcript_cuts[t]
		transcript_cut_count = transcript_cuts[t]

		if transcript_cut_count > 0:
			transcripts_to_remove.append(t)

			new_transcripts = transcript_cut_count + 1

			new_transcript_names = []

			oristart = transcripts[t][0]
			oriend = transcripts[t][1]

			ori_exons = transcript_exons[t]

			for i in range(1, new_transcripts + 1):
				newname = t.split('.')[0] + '.' + t.split('.')[1] + '-' + str(i) + "." + t.split('.')[2]
				new_transcript_names.append(newname)

			transcript_exon_list = {t: [] for t in new_transcript_names}

			new_exon_ids = {}

			transcript_count = 1

			for e in ori_exons:
				seen = []
				seene = []
				#if no cut
				if e not in exons_to_cut:
					#add it to current transcript
					curr = t.split('.')[0] + '.' + t.split('.')[1] + '-' + str(transcript_count) + "." + t.split('.')[2]
					seen.append(curr)
					start = exons[e][0]
					end = exons[e][1]
					olddefline = exons[e][2]
					origname = '.'.join(oriname.split('.')[0:2])
					currgname = ".".join(curr.split('.')[0:2])
					count = 'exon_number "' + str(len(transcript_exon_list[curr]) + 1) + '"'
					ename = curr + '-' + count
					newdefline = olddefline.replace(origname, currgname)
					exons[ename] = [start, end, newdefline]
					transcript_exon_list[curr].append(ename)
					seene.append(ename)

				elif e in exons_to_cut: #and len(exons_to_cut[e]) == 1:

					oldstart = exons[e][0]
					oldend = exons[e][1]
					olddefline = exons[e][2]
					oldcount = 'exon_number "' + olddefline.split('\t')[-1].split("exon_number")[1].split(';')[0].strip(' ').strip('"') + '"'
					

					for cut in exons_to_cut[e]:
						#seen = [] #list of transcripts that have been edited
						#while i < len(exons_to_cut[e]):
						#for i in range(0, len(exons_to_cut[e]) - 1):
						#for each cut make two new
						#cut = exons_to_cut[e][i]
						curr = t.split('.')[0] + '.' + t.split('.')[1] + '-' + str(transcript_count) + "." + t.split('.')[2]
						e1end =  cut[0]
						oriname = '.'.join(exons[e][2].split('\t')[-1].split(';')[0].split(' ')[-1].strip('"').split('.')[0:2])

						#e1start = 
						if curr in seen:
							curre = transcript_exon_list[curr][-1]
							e1start = exons[curre][0]
							oldend = exons[curre][1]
							olddefline = exons[curre][2]
							e1defline = olddefline.replace(oldend, e1end)
							e1name = curre
							#exons[e1name]

						else:
							curre = e
							e1start = exons[e][0]
							origname = '.'.join(oriname.split('.')[0:2])
							currgname = ".".join(curr.split('.')[0:2])
							#e1count = 'exon_number "' + str(len(transcript_exon_list[curr]) + 1) + '"'
							e1count = str(len(transcript_exon_list[curr]) + 1)
							#e1defline = olddefline.replace(oldcount, e1count).replace(oldstart, e1start).replace(oldend, e1end).replace(origname, currgname)
							e1defline = "\t".join(olddefline.split('\t')[0:3]) + '\t' + e1start + '\t' + e1end + '\t' + '\t'.join(olddefline.split('\t')[5:-1]) + '\tgene_id "' + currgname + '"; transcript_id "' + currgname + '.1"; exon number "' + e1count + '"; ' + ' '.join(olddefline.split(' ')[-2:])
							e1name = curr + '-' + e1count
							#t_exons_to_update[curr].append(e1name)
							if curr not in transcript_exon_list:
								transcript_exon_list[curr] = []
							if e1name not in transcript_exon_list[curr]:
								transcript_exon_list[curr].append(e1name)
							#t_exons[e1name] = [e1start, e1end, e1defline]
						exons[e1name] = [e1start, e1end, e1defline]

						seen.append(curr)
						seene.append(e1name)
							#e1end =  cut[0]
							#print("curr1:", curr)

							#origname = '.'.join(oriname.split('.')[0:2])
							#currgname = ".".join(curr.split('.')[0:2])

							#e1count = 'exon_number "' + str(len(transcript_exon_list[curr]) + 1) + '"'
							#e1defline = olddefline.replace(oldcount, e1count).replace(oldstart, e1start).replace(oldend, e1end).replace(origname, currgname)
							#e1name = curr + '-' + e1count
							#exons[e1name] = [e1start, e1end, e1defline]
						if e1name not in transcript_exon_list[curr]:
							transcript_exon_list[curr].append(e1name)
							
							
							

						transcript_count += 1
						curr = t.split('.')[0] + '.' + t.split('.')[1] + '-' + str(transcript_count) + "." + t.split('.')[2]


						if curr not in seen: #and e not in seene:

							oldstart = exons[e][0]
							e2start = cut[1]
							e2stop = oldend

							e2count = str(len(transcript_exon_list[curr]) + 1)
							origname = '.'.join(oriname.split('.')[0:2])
							currgname = ".".join(curr.split('.')[0:2])

							#e2defline = olddefline.replace(oldcount, e2count).replace(oldstart, e2start).replace(oldend, e2stop).replace(origname, currgname)
							e2defline = "\t".join(olddefline.split('\t')[0:3]) + '\t' + e2start + '\t' + e2stop + '\t' + '\t'.join(olddefline.split('\t')[5:-1]) + '\tgene_id "' + currgname + '"; transcript_id "' + currgname + '.1"; exon number "' + e2count + '"; ' + ' '.join(olddefline.split(' ')[-2:])
							e2name = curr + '-' + e2count
						else:
							e2name = e
							e2start = cut[1]
							e2stop = oldend
							e2defline = exons[e][2].replace(oldstart, e2start).replace(oldend, e2stop)

						exons[e2name] = [e2start, e2stop, e2defline]
						if curr not in transcript_exon_list:
							transcript_exon_list[curr] = []
						if e2name not in transcript_exon_list[curr]:
							transcript_exon_list[curr].append(e2name)

						seen.append(curr)
						seene.append(e2name)

						#for e in t_exons:
						#	exons[e]
						#update all exons
						#for e in t_exons_to_update:

					#add all the stuff to the transcript exon list


			#edit transcript deflines and start and stops
			first_t = new_transcript_names[0]
			last_t = new_transcript_names[-1]
			for trans in new_transcript_names:
				#get new starts and stops
				if trans != first_t:

					tstart = exons[transcript_exon_list[trans][0]][0]
				else:
					tstart = oristart
				if trans != last_t:
					tend = exons[transcript_exon_list[trans][-1]][1]
				else:
					tend = oriend

				#replace oristart, oriend, and oriname
				origname = '.'.join(t.split('.')[0:2])
				currgname = ".".join(trans.split('.')[0:2])
				old_defline = transcripts[t][3]
				new_defline = old_defline.replace(oristart, tstart).replace(oriend, tend).replace(origname, currgname)
				chrom = transcripts[t][2]
				transcripts[trans] = [tstart, tend, chrom, new_defline]
				transcript_exons[trans] = []
				for e in transcript_exon_list[trans]:
					transcript_exons[trans].append(e)




	log = open(sys.argv[4], 'w')
	for l in logwrite:
		log.write(l)
	log.close()

	gff = open(sys.argv[3], 'w')
	for t in transcripts:
		if t not in transcripts_to_remove:
			gff.write(transcripts[t][-1] + '\n')
			for e in transcript_exons[t]:
				gff.write(exons[e][-1] + '\n')
		
	gff.close()


		




if __name__ == "__main__":
  main(sys.argv)
