#!/bin/bash

prefix=$1
braker=$2
junc=$3
bam=$4


echo "adding introns"
gt gff3 -addintrons -retainids $braker > $prefix"_braker_with_introns.gff3"
grep -P "\tintron\t" $prefix"_braker_with_introns.gff3" | cut -f 1,4,5 > $prefix"_introns.bed"
cut -f 1,2,3 $junc > $prefix"_rnaseq_supported_juncs.bed"

echo "getting good introns"
grep -f $prefix"_introns.bed" $prefix"_rnaseq_supported_juncs.bed" > $prefix"_braker_introns_rnaseq_supported.bed"


echo "cleaning inrons"
/pollard/data/projects/alind/protist_sequencing/Blastocystis_BT1/bt1_lib1/assembly/blastannotate/scripts/remove_bad_introns.py $prefix"_braker_with_introns.gff3" $prefix"_braker_introns_rnaseq_supported.bed"  > $prefix"_braker_supported_introns.gff3"

echo "getting exon depth"
grep -P "\texon\t" $prefix"_braker_supported_introns.gff3" | cut -f 1,4,5 > $prefix"_braker_fixed_introns_exons.bed"
samtools depth -a -g SECONDARY,QCFAIL,DUP -b $prefix"_braker_fixed_introns_exons.bed" $bam | grep -P "\t0$" > $prefix"_braker_fixed_exons_zerodepth.txt"

echo "splitting exons"

/pollard/data/projects/alind/protist_sequencing/Blastocystis_BT1/bt1_lib1/assembly/blastannotate/split_exons.py $prefix"_braker_supported_introns.gff3" $prefix"_braker_fixed_exons_zerodepth.txt" > $prefix"_braker_supported_introns_split_exons.gff3"

