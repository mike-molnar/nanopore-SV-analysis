#! /bin/bash
CHROMOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

for CHROMOSOME in $CHROMOSOMES; do
    wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/$CHROMOSOME.fa.gz
    gunzip $CHROMOSOME.fa.gz
done
cat chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa chr18.fa chr19.fa chr20.fa chr21.fa chr22.fa chrX.fa chrY.fa > GRCh38.fa
rm chr*.fa
samtools faidx GRCh38.fa
