#!/bin/bash

#SBATCH --mail-user=duvanel@stanford.edu  
#SBATCH --mail-type=ALL
#SBATCH --qos=normal
#SBATCH --job-name=BP.Distance.GRM
#SBATCH --time=48:00:00
#SBATCH --output=BP.Distance.GRM.out
#SBATCH --error=BP.Distance.GRM.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000

module load plink

echo ""
echo "Compute IBS distance matrix"
echo ""

######plink --bfile /home/duvanel/BP/Data/GenoData1M/MendelAnalysis3 --distance square 1-ibs flat-missing
########plink --bfile /home/duvanel/BP/Data/GenoData1M/MendelAnalysis3 --matrix
plink --bfile /home/duvanel/BP/Data/GenoData1M/MendelAnalysis3 --genome
