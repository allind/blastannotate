#!/bin/bash

fwd=$1
rev=$2
genome=$3
threads=$4

prefix=$(echo $genome | rev | cut -f 1 -d "/" | rev | cut -f 1 -d ".")



#first alignment
STAR --readFilesCommand zcat --outSAMmapqUnique 10 --outFilterMismatchNmax 0 --alignIntronMax 40000 --alignMatesGapMax 40000 --alignIntronMin 15 --outSJfilterCountUniqueMin -1 10 10 10 --outSAMtype BAM SortedByCoordinate --genomeDir $prefix --outFileNamePrefix $prefix"_firstpass_" —readFilesIn $fwd $rev --runThreadN $threads

#filter junctions
python scripts/filter_sj_tabfile.py $prefix"_passone_SJ.out.tab" 25 > $prefix"_passone_min25_SJ.out.tab"

#re-align with strict junctions
STAR --readFilesCommand zcat  --outSAMmapqUnique 10 --outFilterMismatchNmax 0 --alignIntronMax 40000 --alignMatesGapMax 40000 --alignIntronMin 15 --alignSJoverhangMin 1000 --genomeDir $prefix --outFileNamePrefix $prefix"_secondpass_" --readFilesIn $fwd $rev --runThreadN $threads --sjdbFileChrStartEnd $prefix"_passone_min25_SJ.out.tab"

#process bam
#samtools sort -o $prefix"_secondpass_"


