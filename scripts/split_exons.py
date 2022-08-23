#! /usr/bin/env python
#usage:script.py [original gff] [original exons depth] [original rnaseq tab file]

import sys, copy


def parse_exons(gff):
    exons = []
    for line in gff:
      feat = line.split('\t')[2]
      if feat == "exon":
        #exonname = whatever?
        #exonname = "blah"
        start = line.split('\t')[3]
        end = line.split('\t')[4]
        #exons.append([exonname, start, end])
        exons.append([start, end, 0])
      #elif feat == "CDS":
       # start = line.split('\t')[3]
        #end = line.split('\t')[4]
       # offset = int(line.split('\t')[7])
        #exons[[start, end, 0]] = [start, end, offset]
    return exons

def parse_introns(gff):
  introns = []
  for line in gff:
    feat = line.split('\t')[2]
    if feat == "intron":
      #exonname = whatever?
      #intronname = "blah"
      start = line.split('\t')[3]
      end = line.split('\t')[4]
      #introns.append([intronname, start, end])
      introns.append([start,end])
  return introns

class Transcript:

  def __init__(self, defline, chrom, exons):#, introns):
    self.defline = defline
    self.chrom = chrom
    self.exons = exons
    #self.introns = introns

def main(argv):


  #create all the transcripts
  curr = []
  original_transcripts = []

  for line in open(sys.argv[1]):
    if "#" not in line[0]:
      line = line.strip('\n')
      feat = line.split('\t')[2]
      if feat == "transcript":
        if curr != []:
          exons = parse_exons(curr)
          introns = parse_introns(curr)
          defline = curr[0]
          chrom = defline.split('\t')[0]
          trans = Transcript(defline, chrom, exons)#, introns)
          original_transcripts.append(trans)

          curr = []
        curr.append(line)
      else:
        curr.append(line)
  #add the last one
  exons = parse_exons(curr)
  introns = parse_introns(curr)
  defline = curr[0]
  chrom = defline.split('\t')[0]
  trans = Transcript(defline, chrom, exons)#, introns)
  original_transcripts.append(trans)

  uncovered = {}


  for line in open(sys.argv[2]):
    line = line.strip('\n')
    contig = line.split('\t')[0]
    pos = line.split('\t')[1]
    if contig not in uncovered:
      uncovered[contig] = []
    uncovered[contig].append(pos)

  exon_cov_transcripts = []

  for trans in original_transcripts:
    chrom = trans.chrom
    tstart = trans.defline.split('\t')[3]
    tend = trans.defline.split('\t')[4]


    all_exons_status = []

    exons_to_cut = []

    for exon in trans.exons:
      start = exon[0]
      stop = exon[1] 
      elen = int(stop) - int(start) - 1
      ename = chrom + '-' + start + '-' + stop
      span = [str(e) for e in range(int(start), int(stop))]


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
          run = [int(b)]
          groups.append(run)
        expect = int(b) + 1


      #remove single base drops and starts & stops w/in 5 bases of start/stop
      if groups != [[]]:
        for g in groups:
          removed = False
          if len(g) == 1:
            groups.remove(g)
            removed = True
          if (int(g[-1]) - 6) < int(tstart) and removed == False:
            groups.remove(g)
            removed = True
          if (int(g[0]) + 6) > int(tend) and removed == False:
            groups.remove(g)
            removed = True
      else:
        groups = []


      collapsed_groups = []


      if len(groups) > 0:
        for g in groups:
            collapsed_groups.append([g[0], g[-1], len(g)])
      else:
        all_exons_status.append([ename, "no change"])

      if len(groups) > 0:

        close = False
        for g in collapsed_groups:
          if g[2] > (elen - 10):
            close = True
        if close:
          all_exons_status.append([ename, "cut completely"])
        else:
          if len(collapsed_groups) > 1:
            all_exons_status.append([ename, "cut multiple times", collapsed_groups])
          else:
            all_exons_status.append([ename, "cut once", collapsed_groups])

          #print(trans.defline)
          #print(elen)
          #print(collapsed_groups)


    tname = trans.defline.split('\t')[-1].split(';')[0].split('=')[1]
    exoncount = len(all_exons_status)
    #print(tname, "exon total:", len(all_exons_status))
    #print(all_exons_status)
    vals = [e[1] for e in all_exons_status]
    valset = list(set(vals))
    if len(valset) == 1 and "cut completely" in vals:
      #you do not create a new transcript. you delete the gene entirely.
      exon_cov_transcripts = []
      #elif len(valset) > 1 or "no change" not in vals: #there is a change
    elif valset != "no change":
      new_exons = {}
      exon_transcript_link = {}
      transnum = 1
      for e in all_exons_status:

        if e[1] == "cut completely":
          #no transcripts have this exon
          #continue?
          transnum += 1
        elif e[1] == "cut once":
          #cut it once

          e1start = e[0].split('-')[1]
          e1stop = e[2][0][0]

          e1name = chrom + '-' + str(e1start) + '-' + str(e1stop)

          new_exons[e1name] = [e1start, e1stop]
          exon_transcript_link[e1name] = transnum
          transnum += 1

          e2start = e[2][0][1]
          e2stop = e[0].split('-')[2]
          e2name = chrom + '-' + str(e2start) + '-' + str(e2stop)
          exon_transcript_link[e2name] = transnum
          new_exons[e2name] = [e2start, e2stop]

        elif e[1] == "cut multiple times": #this is the hard one
          newstarts = []
          newstops = []

          oristart = e[0].split('-')[1]
          oriend = e[0].split('-')[2]

          newstarts.append(oristart)

          for g in e[2]:
            stop = g[0]
            start = g[1]
            newstops.append(stop)
            newstarts.append(start)
          newstops.append(oriend)

          #for start in newstarts:
          for i in range(0, len(newstarts)):

            ename = chrom + '-' + str(newstarts[i]) + '-' + str(newstops[i])

            new_exons[ename] = [newstarts[i], newstops[i]]
            exon_transcript_link[ename] = transnum
            transnum += 1
          transnum = transnum -1

        elif e[1] == "no change":
          start = e[0].split('-')[1]
          stop = e[0].split('-')[2]
          new_exons[e[0]] = [start,stop]

          exon_transcript_link[e[0]] = transnum

      oldtname = trans.defline.split('\t')[-1].split('ID=')[1].split(';')[0]
      for i in range(0, transnum):
        #this many new transcripts
        if len(exons) > 0:
          newtname = oldtname.split('.')[0] + '-' + str(i + 1) + '.' + oldtname.split('.')[1]

          exons = [new_exons[e] for e in exon_transcript_link if exon_transcript_link[e] == i + 1]
          #sort exons
          #print(trans.defline)
          #print(all_exons_status)
          #print(exons)
          
          sorted_exons = sorted(exons, key = lambda x: int(x[0]))
          #print(sorted_exons)
          #print(trans.defline)
          #print(transnum)
          #print(sorted_exons)

          newstart = sorted_exons[0][0]
          newend = sorted_exons[0][1]
          newdefline = "\t".join(trans.defline.split('\t')[0:3]) + '\t' + str(newstart) + '\t' + str(newend) + '\t' + "\t".join(trans.defline.split('\t')[5:8]) + '\tID=' + newtname + ';Parent=' + newtname.split('.')[0] 
          newt = Transcript(newdefline, chrom, sorted_exons)
          exon_cov_transcripts.append(newt)

    else:
      #no changes!
      exon_cov_transcripts.append(t)



    for t in exon_cov_transcripts:
      print(t.defline) #commented for debug
      eID = t.defline.split('\t')[-1].replace("Parent=", "geneID=").replace("ID=", "Parent=")
      middle = "\t".join(t.defline.split('\t')[5:8])
      source = t.defline.split('\t')[1]
      for e in t.exons: #commented for debug
        print(t.chrom + '\t' + source + '\texon\t' + str(e[0]) + '\t' + str(e[1]) + '\t' + middle + '\t' + eID) #commented for debug

    #print()

    #print(vals)
    #valset = list(set(vals))
    #nochange = str(vals.count("no change"))
    #cutonce = str(vals.count("cut once"))
    #cutall = str(vals.count("cut completely"))
    #cutmulti = str(vals.count("cut multiple times"))
    #print(all_exons_status)
    #cutall_runs = False
    #cutall_run_ids = []
    #cutall_individuals = []
    # if int(cutall) > 1:
    #   #cutalls = vals.index("cut completely")
    #   cutalls = [e for e in range(0, len(vals) -1 ) if vals[e] == "cut completely"]
    #   run = []
    #   groups = [run]
    #   expect = None 
    #   for b in cutalls:
    #     if (int(b) == expect) or (expect is None):
    #       run.append(b)
    #     else:
    #       run = [b]
    #       groups.append(run)
    #     expect = int(b) + 1

    #   #cutall_groups = groups
    #   print(groups)

      #for g in collapsed_groups:

      #if cutall_runs:
        #print()
      #  cutall_print = ['-'.join(g) for g in cutall_run_ids]
        #print(all_exons_status)
      #  print(tname + '\t' + nochange + '\t' + cutonce + '\t' + cutall + '\t' + cutmulti + '\t' + str(exoncount) + '\t' + ','.join(cutall_print))

    #if int(nochange) != exoncount:

      #print(tname + '\t' + nochange + '\t' + cutonce + '\t' + cutall + '\t' + cutmulti + '\t' + str(exoncount))
        




 



  #drop any zero coverage exons



  #print exons bed for samtools depth

  #run samtools depth

  #read samtools depth into memory

  #split exons with zero coverage. if the first exon 








if __name__ == "__main__":
  main(sys.argv)
