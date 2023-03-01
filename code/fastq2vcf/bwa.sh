#!/bin/bash -l
#SBATCH --job-name=verti_bwa
#SBATCH --output=verti_bwa3.out
#SBATCH --partition=standard
#SBATCH --account=aafc_pilot
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=44
#SBATCH --mem-per-cpu=4000M
#SBATCH --comment="image=nrc/nrc_all_default_ubuntu-18.04-amd64_latest"
 

source ~/miniconda3/etc/profile.d/conda.sh

conda activate wgs

fastqPath=/gpfs/fs3/aafc/pilot/aafc_sjsr/verti_can/data/trim/list3.txt
index=/gpfs/fs3/aafc/pilot/aafc_sjsr/verti_can/references/genome/Verticillium_dahliae.ASM15067v2.dna.toplevel.fa

while read -r file
do

name=$(echo $file | sed 's|_R1.fastq.gz||' | sed 's|/gpfs.*trim/||')
R2=$(echo $file | sed 's|_R1.|_R2.|')

bwa mem \
-M \
-t 24 \
-R $(echo "@RG\tID:$name\tSM:$name\tPL:ILLUMINA") \
$index \
$file \
$R2 | samtools sort - -@ 24 -n -m 3G -o /gpfs/fs3/aafc/pilot/aafc_sjsr/verti_can/data/bam/${name}.bam

done < $fastqPath