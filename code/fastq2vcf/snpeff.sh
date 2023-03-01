#!/bin/bash 
#
# snpeff gpsc
#
#SBATCH --job-name=verti_snpeff
#SBATCH --output=verti_test_snpeff.out
#SBATCH --partition=standard
#SBATCH --account=aafc_pilot
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=4000M
#SBATCH --comment="image=nrc/nrc_all_default_ubuntu-18.04-amd64_latest"


source ~/miniconda3/etc/profile.d/conda.sh

conda activate wgs

vcf=/gpfs/fs3/aafc/pilot/aafc_sjsr/verti_can/analysis/test_wgs_freebayes.vcf


out_name=/gpfs/fs3/aafc/pilot/aafc_sjsr/verti_can/analysis/verti_test_freebayes

snpEff eff \
-i vcf \
-t 20 \
Verticillium_dahliae \
$vcf > ${out_name}.snpeff.vcf








