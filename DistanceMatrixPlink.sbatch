#!/bin/bash

#SBATCH --mail-user=duvanel@stanford.edu  
#SBATCH --mail-type=ALL
#SBATCH --qos=normal
#SBATCH --job-name=BP_IBS_Distance
#SBATCH --time=48:00:00
#SBATCH --output=BP_IBS_Distance.out
#SBATCH --error=BP_IBS_Distance.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000

module load plink

echo ""
echo "Compute IBS distance matrix with plink"
echo ""

## plink --bfile /home/duvanel/BP/Data/GenoData1M/MendelAnalysis3 --distance square 1-ibs flat-missing
## plink --bfile /home/duvanel/BP/Data/GenoData1M/MendelAnalysis3 --matrix
plink --bfile /home/duvanel/BP/Data/GenoData1M/MendelAnalysis3 --genome
