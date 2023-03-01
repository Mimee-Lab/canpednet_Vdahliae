#!/bin/bash -l
#SBATCH --job-name=verti_fastp
#SBATCH --output=verti_fastp.out
#SBATCH --partition=standard
#SBATCH --account=aafc_pilot
#SBATCH --time=05:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2000M
#SBATCH --comment="image=nrc/nrc_all_default_ubuntu-18.04-amd64_latest"
 

source ~/miniconda3/etc/profile.d/conda.sh

conda activate wgs

fastqPath=/gpfs/fs3/aafc/pilot/aafc_sjsr/verti_can/data/raw/r1.path.fastq

while read -r file
do

name=$(echo $file | sed 's|/gpfs.*IDT_i5_[0-9]*.||' | sed 's|_R1.fastq.gz||')
R2=$(echo $file | sed 's|_R1.|_R2.|')

fastp \
-i $file \
-I $R2 \
-o /gpfs/fs3/aafc/pilot/aafc_sjsr/verti_can/data/trim/${name}_R1.fastq.gz \
-O /gpfs/fs3/aafc/pilot/aafc_sjsr/verti_can/data/trim/${name}_R2.fastq.gz \
-h /gpfs/fs3/aafc/pilot/aafc_sjsr/verti_can/data/trim/${name}_report.html \
-j /gpfs/fs3/aafc/pilot/aafc_sjsr/verti_can/data/trim/${name}_report.json \
-r \
-W 4 \
-M 20 \
--thread 12

done < $fastqPath