#!/bin/bash
#run_mikado.sh [prefix for input/output files] [reference] [junctions file] [threads]
#have prepped:
prefix=$1
ref=$2
junc=$3
threads=$4
scoringfile=$5

#make config file and list file
mikado configure --list transcripts.txt --reference $ref --mode permissive --scoring $scoringfile --junctions $junc -bt uniprot_refseq_cdhit100.fasta $prefix"_configuration.yaml"
#replace stuff with prefix
perl -pi -e 's/mikado/'$prefix'_mikado/g' $prefix"_configuration.yaml"

#mikado prepare

mikado prepare --json-conf $prefix"_configuration.yaml"
#mv mikado_prepared.fasta $prefix"_mikado_prepared.fasta"
#mv mikado_prepared.gtf $prefix"_mikado_prepared.gtf"

#homology - run these separately, then come back to this

#Run Transdecoder
TransDecoder.LongOrfs -t $prefix'_mikado_prepared.fasta'
TransDecoder-TransDecoder-v5.5.0/TransDecoder.Predict -t $prefix'_mikado_prepared.fasta'


#blast
diamond blastx --outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop" --db uniprot_refseq_cdhit100.fasta --out $prefix"_mikado_prepared.blast.tsv" --query $prefix"_mikado_prepared.fasta" --threads $threads
#blastx -max_target_seqs 5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop" -db /pollard/data/projects/alind/protist_sequencing/reference_euk_proteins/uniprot_refseqprot_blasto/uniprot_refseq_blasto_cdhit100.fasta -out $prefix"_mikado_prepared.blast.tsv" -query $prefix"_mikado_prepared.fasta" -num_threads $threads

#serialize
mikado serialise --json-conf $prefix"_configuration.yaml" --xml $prefix"_mikado_prepared.blast.tsv" --orfs $prefix"_mikado_prepared.fasta.transdecoder.gff3"  --blast_targets uniprot_refseq_cdhit100.fasta --junctions $junc

#pick

mikado pick --configuration $prefix"_configuration.yaml" --subloci-out $prefix"_mikado.subloci.gff3"

mv mikado.loci.gff3 $prefix"_mikado.loci.gff3"
mv mikado.loci.metrics.tsv $prefix"_mikado.loci.metrics.tsv"
mv mikado.loci.scores.tsv $prefix"_mikado.loci.scores.tsv"
gffread -g $ref -x $prefix"_mikado_cds.fna" $prefix"_mikado.loci.gff3"
