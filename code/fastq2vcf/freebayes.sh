#!/bin/bash 
#
# freebayes gpsc
#
#SBATCH --job-name=verti_vcf
#SBATCH --output=verti_vcf.out
#SBATCH --partition=standard
#SBATCH --account=aafc_pilot
#SBATCH --time=18:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=4000M
#SBATCH --comment="image=nrc/nrc_all_default_ubuntu-20.04-amd64_latest"

source ~/miniconda3/etc/profile.d/conda.sh

conda activate wgs

export TMPDIR=/gpfs/fs7/aafc/pilot/aafc_sjsr/tmp

bamList=/gpfs/fs7/aafc/pilot/aafc_sjsr/verti_can/data/bam/bamlist.txt

out_vcf=/gpfs/fs7/aafc/pilot/aafc_sjsr/verti_can/analysis/verti_wgs2_freebayes.vcf

freebayes-parallel <(fasta_generate_regions.py /gpfs/fs7/aafc/pilot/aafc_sjsr/verti_can/references/genome/old_vdls17/Verticillium_dahliae.ASM15067v2.dna.toplevel.fa.fai 100000) \
36 -f /gpfs/fs7/aafc/pilot/aafc_sjsr/verti_can/references/genome/old_vdls17/Verticillium_dahliae.ASM15067v2.dna.toplevel.fa \
-L $bamList \
--ploidy 1 \
--strict-vcf \
--min-alternate-fraction 0.2 \
--use-best-n-alleles 4 \
--min-alternate-total 30 \
--min-alternate-count 5 \
--min-coverage 500 > ${out_vcf}





