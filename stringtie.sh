#!/bin/bash

bam=$1
prefix=$2
threads=$3

#~/applications/stringtie-2.2.1.Linux_x86_64/stringtie $bam --rf -m 100 -p $threads -v -o $prefix"_min25_stringtie.gtf"

#depth file
grep -P "\texon\t" $prefix"_min25_stringtie.gtf" | cut -f 1,3,4 > $prefix"_exons.bed"
samtools depth -aa -g SECONDARY,QCFAIL,DUP -b $prefix"_exons.bed" $bam | grep -P "\t0$" > $prefix"_exon_zeros.txt"

#split
python scripts/stringtie_depth_split.py $prefix"_exon_zeros.txt" $prefix"_min25_stringtie.gtf"  $prefix"_min25_stringtie_split.gtf" $prefix"_min25_stringtie_split.log"

#filter short
python scripts/filter_short_transcripts.py $prefix"_min25_stringtie_split.gtf" $prefix"_min25_stringtie_split_filtered.gtf" $prefix"_min25_stringtie_split_filtered.log"
