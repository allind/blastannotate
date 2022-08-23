#!/bin/bash

gff=$1
genome=$2
fwd=$3
rev=$4
prefix=$5
threads=$6

echo $prefix

grep "#PROT" $gff > $prefix"_pasa_prot.fasta"
perl -pi -e 's/#PROT />/' $prefix"_pasa_prot.fasta"
perl -pi -e 's/\t/\n/' $prefix"_pasa_prot.fasta"
#for some reason gffread doesn't work like this but the rest does
#gffread -w $prefix"_spliced_exons.fa" -g $genome $gff3

#run kallisto
kallisto index -i $prefix"_spliced_exons.idx" $prefix"_spliced_exons.fa"

kallisto quant -t $threads -i $prefix"_spliced_exons.idx" -o $prefix"_kallisto" -b 100 $fwd $rev

/pollard/data/projects/alind/protist_sequencing/Blastocystis_BT1/bt1_lib1/assembly/blastannotate/scripts/highest_abundance_isoform.py $prefix"_kallisto/abundance.tsv" > $prefix"_kallisto/highest_iso.txt"
/pollard/data/projects/alind/protist_sequencing/Blastocystis_BT1/bt1_lib1/assembly/blastannotate/scripts/get_highest_isoforms_from_protfile.py $prefix"_kallisto/highest_iso.txt" $prefix"_pasa_prot.fasta" > $prefix"_pasa_prot_highest_iso.fasta"
