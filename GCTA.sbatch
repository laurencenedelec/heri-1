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

module load gcta

echo ""
echo "Compute genetic relationship matrix (GRM) between pairs of individuals"
echo ""

##gcta64  --bfile /home/duvanel/BP/Data/GenoData1M/MendelAnalysis3  --make-grm  --autosome --out BP.GRM --thread-num 16
gcta64  --bfile /home/duvanel/BP/Data/GenoData1M/MendelAnalysis3 --autosome --maf 0.01 --grm-adj 0 --grm-cutoff 0.025 --make-grm  --out BP.GRM --thread-num 16
