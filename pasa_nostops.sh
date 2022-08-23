#!/bin/bash
gff=$1

#grep PROT $1 > pasa_prots.fasta
#perl -pi -e 's/#PROT />/' pasa_prots.fasta
#perl -pi -e 's/\t/\n/' pasa_prots.fasta

/pollard/data/projects/alind/protist_sequencing/Blastocystis_BT1/bt1_lib1/assembly/blastannotate/scripts/pasa_nostops.py pasa_prots.fasta
