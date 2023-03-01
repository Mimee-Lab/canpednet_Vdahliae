#!/bin/bash 
 

source ~/miniconda3/etc/profile.d/conda.sh

conda activate wgs

bamPath=/gpfs/fs3/aafc/pilot/aafc_sjsr/verti_can/data/bam/bamlist.txt

while read -r file
do

name=$(echo $file | sed 's|.bam$||')

samtools flagstat ${name}.fixmate.markdup.sort.bam > ${name}.fixmate.markdup.sort.bam.flagstat


done < $bamPath