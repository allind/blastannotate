#!/bin/bash

braker=$1
junc=$2

echo $braker
echo $junc

#write gff3
gffread --keep-genes --keep-comments -o braker_fixgff.gff3 $braker

#add a header
echo "##gff-version 3" > x
cat x braker_fixgff.gff3 > y
mv y braker_fixgff.gff3
rm x

#add introns
gt gff3 -o "braker_with_introns.gff3" -addintrons -retainids -sort braker_fixgff.gff3


#add semicolons to the introns
perl -pi -e 's/$/;/' braker_with_introns.gff3

#intron bedfile
grep -P "\tintron\t" "braker_with_introns.gff3" |  cut -f 1,4,5,9 > "braker_introns.bed"
cut -f 1,2,3 $junc > rnaseq_introns.bed

bedtools subtract -a braker_introns.bed -b rnaseq_introns.bed > bad_braker_introns.bed

cut -f 4 bad_braker_introns.bed | cut -f 2 -d "=" > bad_braker_transcripts.txt

grep -P "\ttranscript\t" braker_with_introns.gff3 | cut -f 2 -d ";" | cut -f 2 -d "=" | sort | uniq -c | grep " 1 " > single_transcript_genes.txt
perl -pi -e 's/ 1 //' single_transcript_genes.txt
perl -pi -e 's/ *//' single_transcript_genes.txt
perl -pi -e 's/$/;/' single_transcript_genes.txt

cut -f 1 -d "." bad_braker_transcripts.txt > bad_braker_transcript_genenames.txt
perl -pi -e 's/$/;/' bad_braker_transcript_genenames.txt
grep -f bad_braker_transcript_genenames.txt single_transcript_genes.txt	| cut -f 1 -d ';' > bad_braker_genes.txt

perl -pi -e 's/;$//' braker_with_introns.gff3

grep -v -P "\tintron\t" braker_with_introns.gff3 | grep -v -f bad_braker_genes.txt | grep -v -f bad_braker_transcripts.txt > braker_intronfilter.gff3
gffread --keep-genes --keep-comments -o braker_intronfilter_fix.gff3 braker_intronfilter.gff3
