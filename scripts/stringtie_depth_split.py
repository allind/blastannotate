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
	for t in transcript_list:

		tstart = transcripts[t][0]
		tstop = transcripts[t][1]
		chrom = transcripts[t][2]

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
				if base in uncovered[chrom]:
					zeros.append(base)
			
			#cut it at the start and the stop of the zeros

			if len(zeros) > 1:
				exons_to_cut[exon] = [zeros[0], zeros[-1], len(zeros)]

		#add log lines


		exons_to_cut_ids = list(exons_to_cut.keys())
		sizes = [exons_to_cut[e][2] for e in exons_to_cut]
		if len(exons_to_cut_ids) > 0:
			logwrite.append("Zeros in transcript " + t + ". Exons affected total: " + str(len(exons_to_cut_ids)) + " exons. In exon numbers: " + ','.join(exons_to_cut_ids) + '\n')
		#print(t + '\t' + str(len(exons_to_cut_ids)))
		if len(exons_to_cut_ids) > 0:
			transcripts_to_remove.append(t)
			#print(t, exons_to_cut_ids, exons_to_cut)
			#make new transcripts
			oristart = transcripts[t][0]
			oriend = transcripts[t][1]
			new_transcripts = {}
			ori_exons = transcript_exons[t]
			currstart = oristart


			currexons = []
			new = 0
			for i in range(0, len(ori_exons)-1):
				e = ori_exons[i]

				if e in exons_to_cut_ids:
					new += 1
					#make a new transcript
					oldend = exons_to_cut[e][0]
					newstart = exons_to_cut[e][1]


					#make a new transcript
					new_transcript_id = t.split('.')[0] + '.' + t.split('.')[1] + '-' + str(new) + "." + t.split('.')[2]
					new_gene_id = t.split('.')[0] + '.' + t.split('.')[1] + '-' + str(new)
					new_transcript_start = currstart
					new_transcript_end = oldend
					chrom = transcripts[t][2]
					start_replace = transcripts[t][0]
					end_replace = transcripts[t][1]
					gene_replace = '.'.join(t.split('.')[0:2])

					new_transcript_defline = transcripts[t][3].replace(t, new_transcript_id).replace(gene_replace, new_gene_id).replace('\t' + start_replace + '\t', '\t' + new_transcript_start+  '\t').replace('\t' + end_replace + '\t', "\t" + new_transcript_end + '\t')
					transcripts[new_transcript_id] = [new_transcript_start, new_transcript_end, chrom, new_transcript_defline]
					exon_count = 0
					currtranscript_exons = []
					for exon in currexons:
						exon_count += 1
						defline = exons[exon][2]
						old_exon_count = 'exon_number "' + defline.split('\t')[-1].split('exon_number')[1].split(';')[0].strip(' ').strip('"') + '"'
						new_exon_count = 'exon_number "' + str(exon_count) + '"'
						start = exons[exon][0]
						end = exons[exon][1]
						newdef = defline.replace(gene_replace, new_gene_id).replace(t, new_transcript_id).replace(old_exon_count, new_exon_count)
						new_name = exon + '-' + str(new)
						exons[new_name] = [start, end, newdef]
						currtranscript_exons.append(new_name)

					exon_count += 1
					last_exon = e + '-' + str(new)
					last_exon_start = exons[e][0]
					last_exon_oriend = exons[e][1]
					last_exon_end = new_transcript_end
					old_exon_count = 'exon_number "' + defline.split('\t')[-1].split('exon_number')[1].split(';')[0].strip(' ').strip('"') + '"'
					new_exon_count = 'exon_number "' + str(exon_count) + '"'

					defline = exons[e][2]
					newdef = defline.replace(gene_replace, new_gene_id).replace(t, new_transcript_id).replace('\t' + last_exon_oriend + '\t', '\t' + last_exon_end + '\t').replace(old_exon_count, new_exon_count)
					exons[last_exon] = [last_exon_start, last_exon_end, newdef]
					currtranscript_exons.append(last_exon)

					transcript_exons[new_transcript_id] = currtranscript_exons


					currexons = []
					currtranscript_exons = []

					#make new exon and add it
					exon_count += 1
					new_exon_id = e + '-' + str(new + 1)
					new_exon_start = newstart
					new_exon_oristart = exons[e][0]
					new_exon_end = exons[e][1]
					old_exon_count = 'exon_number "' + defline.split('\t')[-1].split('exon_number')[1].split(';')[0].strip(' ').strip('"') + '"'
					new_exon_count = 'exon_number "' + str(exon_count) + '"'
					new_exon_defline = exons[e][2].replace('\t' + new_exon_oristart + '\t', '\t' + new_exon_start + '\t').replace(old_exon_count, new_exon_count)
					exons[new_exon_id] = [new_exon_start, new_exon_end, new_exon_defline]
					currexons.append(new_exon_id)

					currstart = newstart

				else:
					currexons.append(e)

			new += 1
			#make a new transcript
			oldend = transcripts[t][1]
			new_transcript_id = t.split('.')[0] + '.' + t.split('.')[1] + '-' + str(new) + '.' + t.split('.')[2]
			new_gene_id = t.split('.')[0] + '.' + t.split('.')[1] + '-' + str(new)
			new_transcript_start = currstart
			new_transcript_end = oldend

			chrom = transcripts[t][2]
			start_replace = transcripts[t][0]
			end_replace = transcripts[t][1]
			gene_replace = '.'.join(t.split('.')[0:2])

			new_transcript_defline = transcripts[t][3].replace(t, new_transcript_id).replace(gene_replace, new_gene_id).replace('\t' + start_replace + '\t', '\t' + new_transcript_start+  '\t').replace('\t' + end_replace + '\t', "\t" + new_transcript_end + '\t')
			transcripts[new_transcript_id] = [new_transcript_start, new_transcript_end, chrom, new_transcript_defline]

			currtranscript_exons = []
			#renumber exons
			exon_count = 0
			for exon in currexons:
				exon_count += 1

				defline = exons[exon][2]
				old_exon_count = 'exon_number "' + defline.split('\t')[-1].split('exon_number')[1].split(';')[0].strip(' ').strip('"') + '"'
				new_exon_count = 'exon_number "' + str(exon_count) + '"'
				start = exons[exon][0]
				end = exons[exon][1]
				newdef = defline.replace(gene_replace, new_gene_id).replace(t, new_transcript_id).replace(old_exon_count, new_exon_count)
				#new_name = exon + '-' + str(new)
				new_name = ';'.join(newdef.split('\t')[-1].split(';')[1:3])
				exons[new_name] = [start, end, newdef]
				currtranscript_exons.append(new_name)

			exon_count += 1
			e = transcript_exons[t][-1]
			defline = exons[e][2]

			old_exon_count = 'exon_number "' + defline.split('\t')[-1].split('exon_number')[1].split(';')[0].strip(' ').strip('"') + '"'
			new_exon_count = 'exon_number "' + str(exon_count) + '"'

			newdef = exons[e][2].replace(gene_replace, new_gene_id).replace(t, new_transcript_id).replace(old_exon_count, new_exon_count)

			last_exon_start = exons[e][0]

			last_exon_end = new_transcript_end
			last_exon = ';'.join(newdef.split('\t')[-1].split(';')[1:3])
			exons[last_exon] = [last_exon_start, last_exon_end, newdef]
			currtranscript_exons.append(last_exon)

			transcript_exons[new_transcript_id] = currtranscript_exons



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
