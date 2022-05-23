#!/bin/bash
#run_mikado.sh [prefix for input/output files] [reference] [transcriptfile] [scoring file] [junctions file] [threads]
#have prepped:
prefix=$1
ref=$2
trans=$3
score=$4
junc=$4
threads=$5

#conda activate mikado

#transcripts file: transcripts_splitstringtie.txt

#make config file and list file
#mkdir $prefix
#mikado configure --list $trans --reference $ref --mode permissive --scoring $score --junctions $junc -bt /pollard/data/projects/alind/protist_sequencing/reference_euk_proteins/uniprot_refseqprot_blasto/uniprot_refseq_blasto_cdhit100.fasta $prefix"_configuration.yaml"
#replace stuff with prefix
#perl -pi -e 's/mikado/'$prefix'_mikado/g' $prefix"_configuration.yaml"

#mikado prepare

#mikado prepare --json-conf $prefix"_configuration.yaml"
#mv mikado_prepared.fasta $prefix"_mikado_prepared.fasta"
#mv mikado_prepared.gtf $prefix"_mikado_prepared.gtf"

#homology

#Transdecoder
#~/applications/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -t $prefix'_mikado_prepared.fasta'
#~/applications/TransDecoder-TransDecoder-v5.5.0/TransDecoder.Predict -t $prefix'_mikado_prepared.fasta'

#blast
#blastx -max_target_seqs 5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop" -db /pollard/data/projects/alind/protist_sequencing/reference_euk_proteins/uniprot_refseqprot_blasto/uniprot_refseq_blasto_cdhit100.fasta -out $prefix"_mikado_prepared.blast.tsv" -query $prefix"_mikado_prepared.fasta" -num_threads $threads


#serialize
#mikado serialise --json-conf $prefix"_configuration.yaml" --xml $prefix"_mikado_prepared.blast.tsv" --orfs $prefix"_mikado_prepared.fasta.transdecoder.gff3"  --blast_targets /pollard/data/projects/alind/protist_sequencing/reference_euk_proteins/uniprot_refseqprot_blasto/uniprot_refseq_blasto.fasta --junctions $junc

#pick

#mikado pick --configuration $prefix"_configuration.yaml" --subloci-out $prefix"_mikado.subloci.gff3" -od $prefix
