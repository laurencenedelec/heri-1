#!/bin/bash

#SBATCH --mail-user=duvanel@stanford.edu  
#SBATCH --mail-type=ALL
#SBATCH --qos=normal
#SBATCH --job-name=BP.GRM
#SBATCH --time=48:00:00
#SBATCH --output=BP.GRM.out
#SBATCH --error=BP.GRM.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000

module load plink/1.9_b1

echo ""
echo "Compute genetic relationship matrix (GRM) between pairs of individuals"
echo ""

plink --bfile /home/duvanel/BP/Data/GenoData1M/MendelAnalysis3 --make-rel square

