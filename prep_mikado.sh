#!/bin/bash

genome=$1
trinitydn=$2
trinitygg=$3
stringtie=$4
threads=$5

#make gmap db

genome_name=$(echo $genome | rev | cut -f 1 -d "/" | rev | cut -f 1 -d ".")

~/applications/gmap-2021-12-17/bin/gmap_build -D gmap_db/ -d $genome_name $genome
~/applications/gmap-2021-12-17/bin/gmap -t $threads -D gmap_db -d $genome_name  -f gff3_match_cdna --gff3-fasta-annotation=2 $trinitygg > trinity_gg.gff3 2> gmap_tgg.err
#~/applications/gmap-2021-12-17/bin/gmap -t $threads -D gmap_db -d $genome_name  -f gff3_match_cdna --gff3-fasta-annotation=2 $trinitydn > trinity_denovo.gff3 2> gmap_tdn.err

#ln -s $stringtie stringtie.gtf

#cp configs
#cp /pollard/data/projects/alind/protist_sequencing/Blastocystis_BT1/bt1_lib1/assembly/blastannotate/mikado_files/* . 

