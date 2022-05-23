ori=$1
fwd=$2
rev=$3

oriname=$(echo $ori | cut -f 1 -d ".")

echo "polish 1"
bwa index $ori
bwa mem -t 25 -a $ori $fwd > aln_r1_1.sam
bwa mem -t 25 -a $ori $rev > aln_r2_1.sam
~/applications/polypolish/polypolish_insert_filter.py --in1 aln_r1_1.sam --in2 aln_r2_1.sam --out1 filt_r1_1.sam --out2 filt_r2_1.sam
~/applications/polypolish/polypolish $ori filt_r1_1.sam filt_r2_1.sam > $oriname"_poly1.fasta"


echo "polish2"
bwa index $oriname"_poly1.fasta"
bwa mem -t 25 -a $oriname"_poly1.fasta"  $fwd > aln_r1_2.sam
bwa mem -t 25 -a $oriname"_poly1.fasta" $rev > aln_r2_2.sam
~/applications/polypolish/polypolish_insert_filter.py --in1 aln_r1_2.sam --in2 aln_r2_2.sam --out1 filt_r1_2.sam --out2 filt_r2_2.sam
~/applications/polypolish/polypolish $oriname"_poly1.fasta" filt_r1_2.sam filt_r2_2.sam > $oriname"_poly2.fasta"

echo "polish3"
bwa index $oriname"_poly2.fasta"
bwa mem -t 25 -a $oriname"_poly2.fasta"   $fwd > aln_r1_3.sam
bwa mem -t 25 -a $oriname"_poly2.fasta"  $rev > aln_r2_3.sam
~/applications/polypolish/polypolish_insert_filter.py --in1 aln_r1_3.sam --in2 aln_r2_3.sam --out1 filt_r1_3.sam --out2 filt_r2_3.sam
~/applications/polypolish/polypolish $oriname"_poly2.fasta" filt_r1_3.sam filt_r2_3.sam > $oriname"_poly3.fasta"


