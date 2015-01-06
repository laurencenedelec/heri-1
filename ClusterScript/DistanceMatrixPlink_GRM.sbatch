#!/bin/bash

#SBATCH --mail-user=duvanel@stanford.edu  
#SBATCH --mail-type=ALL
#SBATCH --qos=normal
#SBATCH --job-name=BP_Distance_GRM_plink
#SBATCH --time=48:00:00
#SBATCH --output=BP_Distance_GRM_plink.out
#SBATCH --error=BP_Distance_GRM_plink.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000

module load plink/1.9_b1

echo ""
echo "Compute genetic relationship matrix (GRM) between pairs of individuals with plink"
echo ""

plink --bfile /home/duvanel/BP/Data/GenoData1M/MendelAnalysis3 --make-rel square

