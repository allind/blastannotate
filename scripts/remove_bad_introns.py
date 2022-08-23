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
        exons.append([start, end])
        #exons.append([start, end, 0])
      #elif feat == "CDS":
      #  start = line.split('\t')[3]
      #  end = line.split('\t')[4]
        #offset = int(line.split('\t')[7])
        #exons[[start, end]] = [start, end]
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

  def __init__(self, defline, chrom, exons, introns):
    self.defline = defline
    self.chrom = chrom
    self.exons = exons
    self.introns = introns

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
          trans = Transcript(defline, chrom, exons, introns)
          original_transcripts.append(trans)

          curr = []
        curr.append(line)
      else:
        curr.append(line)


  #clean up introns - if an intron is not supported by rnaseq merge the preceding and following exon together

  ok_introns = [line.strip('\n').replace('\t','-') for line in open(sys.argv[2])]
  clean_intron_transcripts = []
  for t in original_transcripts:
    chrom = t.chrom
    bad_introns = []

    for i in t.introns:
      ibad = False
      iname = chrom + '-' + i[0] + '-' + i[1]
      if iname not in ok_introns:
        bad_introns.append(iname)

    if len(bad_introns) > 0:
      tstart = t.defline.split('\t')[3]
      tend = t.defline.split('\t')[4]
      good_introns = []
      for i in t.introns:
        iname = chrom + '-' + i[0] + '-' + i[1]
        if iname not in bad_introns:
          good_introns.append(i)
      new_exons = []

      if len(good_introns) > 0:
      #create new exons
        start = tstart
        for i in range(0, len(good_introns) - 1):
          end = int(good_introns[i][0]) - 1
          new_exons.append([start, end])
          start = int(good_introns[i][1]) + 1
        new_exons.append([str(start), str(tend)])
      else:
        new_exons.append([tstart, tend])

      #print(t.defline)
      #print(t.exons)
      #t.exons = new_exons #update transcript
      
      #print(t.exons)
      #t.introns = good_introns
      new_t = Transcript(t.defline, chrom, new_exons, good_introns)
      clean_intron_transcripts.append(new_t)
    else:
      clean_intron_transcripts.append(t)


 # dest = open()
  seen_count = {}
  for t in clean_intron_transcripts:
    
    source = t.defline.split('\t')[1]
    middle = "\t".join(t.defline.split('\t')[5:8])
    #transid = t.defline.split('\t')[-1].split(';')[0].split('=')[1]
    oritrans = t.defline.split('\t')[-1].split(';')[0].split('=')[1]
    geneid = t.defline.split('\t')[-1].split(';')[1].split('=')[1]
    if geneid not in seen_count:
      transid = geneid + ".t1"
      seen_count[geneid] = 1
    else:
      transid = geneid + '.t' + str(seen_count[geneid] + 1)
      seen_count[geneid] += 1
    print(t.defline.replace(oritrans, transid).replace("geneID=", "Parent="))
    for e in t.exons:
      print(t.chrom + '\t' + source + '\texon\t' + str(e[0]) + '\t' + str(e[1]) + '\t' + middle + '\tParent=' + transid)# + ';gene_id=' + geneid)
    for i in t.introns:
      print(t.chrom + '\t' + source + '\tintron\t' + str(i[0]) + '\t' + str(i[1]) + '\t' + middle + '\tParent=' + transid)# + ';gene_id=' + geneid)





if __name__ == "__main__":
  main(sys.argv)
