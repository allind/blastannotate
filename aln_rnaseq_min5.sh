#!/bin/bash

fwd=$1
rev=$2
genome=$3
threads=$4

prefix=$(echo $genome | rev | cut -f 1 -d "/" | rev | cut -f 1 -d ".")

echo $fwd
echo $rev
echo $genome
echo $threads

#fix the stupid polypolish genome thing
#perl -pi -e 's/polypolish_polypolish_polypolish/polish/' $genome

#make db
#STAR --runMode genomeGenerate --runThreadN $threads --genomeDir $prefix --genomeFastaFiles $genome  --genomeSAindexNbases 11

#first alignment
#STAR --readFilesCommand zcat --outSAMmapqUnique 10 --outFilterMismatchNmax 0 --alignIntronMax 40000 --alignMatesGapMax 40000 --alignIntronMin 15 --outSJfilterCountUniqueMin -1 10 10 10 --genomeDir $prefix --outFileNamePrefix $prefix"_firstpass_" --readFilesIn $fwd $rev --runThreadN $threads

#filter junctions
python /pollard/data/projects/alind/protist_sequencing/Blastocystis_BT1/bt1_lib1/assembly/blastannotate/scripts/filter_sj_tabfile.py $prefix"_firstpass_SJ.out.tab" 5 > $prefix"_firstpass_min5_SJ.out.tab"

#re-align with strict junctions
STAR --readFilesCommand zcat  --outSAMmapqUnique 10 --outFilterMismatchNmax 0 --alignIntronMax 40000 --alignMatesGapMax 40000 --alignIntronMin 15 --alignSJoverhangMin 1000 --genomeDir $prefix --outFileNamePrefix $prefix"_secondpass_min5_" --readFilesIn $fwd $rev --runThreadN $threads --sjdbFileChrStartEnd $prefix"_firstpass_min5_SJ.out.tab"

#sort and compress sam
#samtools sort -@ $threads -o $prefix"_firstpass.sorted.bam" $prefix"_firstpass_Aligned.out.sam"
#rm $prefix"_firstpass_Aligned.out.sam"
samtools sort -@ $threads -o $prefix"_secondpass_min5.sorted.bam" $prefix"_secondpass_min5_Aligned.out.sam"
rm $prefix"_secondpass_min5_Aligned.out.sam"

#samtools index $prefix"_firstpass.sorted.bam"
samtools index $prefix"_secondpass_min5.sorted.bam"

#remove tmpdirs

