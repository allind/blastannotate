#!/bin/bash

gff=$1
bam=$2


prefix=$(echo $gff | cut -f 1,2 -d ".")

gffread -E -T --keep-genes --keep-comments -o $prefix".clean.gtf" $gff
htseq-count $bam $prefix".clean.gtf" > $prefix"_htseqcount.txt" 2> htseqcount.err
