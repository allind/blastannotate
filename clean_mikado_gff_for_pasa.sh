#!/bin/bash

gff=$1
prefix=$(echo $gff | cut -f 1 -d ".")

grep ncRNA_gene $gff | cut -f 9 | cut -f 1 -d ";" | cut -f 2 -d "=" > ncgenes.txt
/pollard/data/projects/alind/protist_sequencing/Blastocystis_BT1/bt1_lib1/assembly/blastannotate/scripts/filter_nc.py ncgenes.txt $gff > $prefix"_coding.gff3"

#validate
~/miniconda3/envs/pasa/opt/pasa-2.4.1/misc_utilities/pasa_gff3_validator.pl $prefix"_coding.gff3"
