#!/bin/bash -l
#SBATCH --job-name=markdup
#SBATCH --output=markdup_bam.out
#SBATCH --partition=standard
#SBATCH --account=aafc_pilot
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=44
#SBATCH --mem-per-cpu=4000M
#SBATCH --comment="image=nrc/nrc_all_default_ubuntu-18.04-amd64_latest"
 

source ~/miniconda3/etc/profile.d/conda.sh

conda activate wgs

bamPath=/gpfs/fs3/aafc/pilot/aafc_sjsr/verti_can/data/bam/bamlist.txt

while read -r file
do

name=$(echo $file | sed 's|.bam$||')

samtools fixmate -m -@ 24 $file ${name}.fixmate.bam
samtools sort -@ 24 -m 3G -o ${name}.sort.bam ${name}.fixmate.bam
samtools markdup -s -@24 ${name}.sort.bam ${name}.fixmate.markdup.sort.bam
samtools index ${name}.fixmate.markdup.sort.bam
samtools flagstat ${name}.fixmate.markdup.sort.bam > ${name}.fixmate.markdup.sort.bam.flagstat

samtools quickcheck ${name}.fixmate.markdup.sort.bam && rm ${name}.bam ${name}.sort.bam ${name}.fixmate.bam \
    || echo "Final bam $name is not valid. Please check the temporary files"
    echo "Final bam is valid. Deleting temporary files"


done < $bamPath