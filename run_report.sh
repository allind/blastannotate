#!/bin/bash

#conda activate EukMS_report

fasta=$1
#outpath=$2
prefix=$2
threads=$3

#busco
#/pollard/home/alind/miniconda3/envs/EukMS_report/bin/busco -m prot --cpu $threads -l eukaryota -f -o $prefix --out_path $outpath -i $fasta

#mmseqs
mkdir mmseqdb_ncbi
/pollard/home/alind/miniconda3/envs/EukMS_report/bin/mmseqs createdb $fasta "mmseqdb_ncbi/"$prefix"_mmseqs_db"

/pollard/home/alind/miniconda3/envs/EukMS_report/bin/mmseqs linsearch "mmseqdb_ncbi/"$prefix"_mmseqs_db" /pollard/data/projects/alind/protist_sequencing/reference_euk_proteins/uniprot_refseqprot_blasto/mmseq_db/uniprot_refseq_blasto_cdhit100 "mmseqdb_ncbi/"$prefix"_ncbi_db0" tmp --split-memory-limit 62G -c 0.3 --cov-mode 1 --remove-tmp-files --threads $threads

/pollard/home/alind/miniconda3/envs/EukMS_report/bin/mmseqs convertalis "mmseqdb_ncbi/"$prefix"_mmseqs_db" /pollard/data/projects/alind/protist_sequencing/reference_euk_proteins/uniprot_refseqprot_blasto/mmseq_db/uniprot_refseq_blasto_cdhit100  "mmseqdb_ncbi/"$prefix"_ncbi_db0" $prefix"_ncbi.m8" --threads $threads


#eggnog
#mkdir eggnog
#/pollard/home/alind/miniconda3/envs/EukMS_report/bin/emapper.py -i $fasta --output $outpath"/eggnog/" --cpu $threads -m diamond --override --data_dir /pollard/data/protein_families/eggnog/
