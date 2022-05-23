#!/bin/bash

nanopore=$1
assembly=$2
prefix=$3

minimap2 -t 50 -ax map-ont $assembly $nanopore 2> minimap1.out 1> $prefix_vs_metaflye.sam
~/applications/racon/build/bin/racon -m 8 -x -6 -g -8 -w 500 -t 50 $nanopore $prefix_vs_metaflye.sam $assembly 2> racon1.out 1> assembly_racon1.fasta

minimap2 -t 50 -ax map-ont assembly_racon1.fasta $nanopore 2> minimap2.out > $prefix_vs_racon1.sam
~/applications/racon/build/bin/racon -m 8 -x -6 -g -8 -w 500 -t 50 $nanopore $prefix_vs_racon1.sam assembly_racon1.fasta 2> racon2.out 1> assembly_racon2.fasta

minimap2 -t 50 -ax map-ont assembly_racon2.fasta $nanopore 2> minimap2.out > $prefix_vs_racon2.sam
~/applications/racon/build/bin/racon -m 8 -x -6 -g -8 -w 500 -t 50 $nanopore $prefix_vs_racon2.sam assembly_racon2.fasta 2> racon3.out 1> assembly_racon3.fasta

minimap2 -t 50 -ax map-ont assembly_racon3.fasta $nanopore 2> minimap3.out > $prefix_vs_racon3.sam
~/applications/racon/build/bin/racon -m 8 -x -6 -g -8 -w 500 -t 50 $nanopore $prefix_vs_racon3.sam assembly_racon3.fasta 2> racon4.out 1> assembly_racon4.fasta
