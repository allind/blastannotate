#!/bin/bash
#trinity is denovo assembled from untrimmed reads

trinity=$1
dbname=$2
dir=$3
genome=$4

#seqclean trinity

~/applications/seqclean-x86_64/seqclean $trinity -c 12

#create config
cp /pollard/data/projects/alind/protist_sequencing/Blastocystis_BT1/bt1_lib1/assembly/blastannotate/pasa_config/*.config .

dirregex=$(echo $dir | perl -pi -e 's/\//\\\//g')

perl -pi -e 's/db/'$dbname'/' *.config
perl -pi -e 's/path/'$dirregex'/' *.config
#change path manually

ln -s $genome .

#launch pasa directly bc conad env change. cmd is:

#~/miniconda3/pkgs/pasa-2.4.1-h1b792b2_1/opt/pasa-2.4.1/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g genome.fasta -t trinity.fasta.clean -T -u trinity.fasta --ALIGNERS blat --CPU 24
