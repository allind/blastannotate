#!/bin/bash

bam=$1
prefix=$2
threads=$3

stringtie $bam --rf -m 100 -p $threads -v -o $prefix"_stringtie.gtf"

#depth file
grep -P "\texon\t" $prefix"_stringtie.gtf" | cut -f 1,4,5 > $prefix"_exons.bed"
samtools depth -aa -g SECONDARY,QCFAIL,DUP -b $prefix"_exons.bed" $bam | grep -P "\t0$" > $prefix"_exon_zeros.txt"

#split

python scripts/stringtie_depth_split_refactored.py $prefix"_depth_exon_below3.txt" $prefix"_stringtie.gtf"  $prefix"_stringtie_split_below3.gtf" $prefix"_stringtie_split_below3.log"

#filter short
python scripts/filter_short_transcripts.py $prefix"_stringtie_split.gtf" $prefix"_stringtie_split_filtered.gtf" $prefix"_stringtie_split_filtered.log"
