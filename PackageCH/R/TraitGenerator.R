
####**********************************************************************
####  Written and Developed by: 
####**********************************************************************
####
####    Julien Duvanel, Copyright 2014
####    email: duvanel@stanford.edu
####
####**********************************************************************
####**********************************************************************

####**********************************************************************
####
####  Generate traits
####
####**********************************************************************

# bim file
files <- "/Users/julien/Dropbox/Ecole/EPFL/5eme annee/ma4 - pdm/PackageCH/data/MendelAnalysis3"

# Load data
bim.data <- read.table(file = paste0(files, ".bim"))
nbr.snps <- nrow(bim.data)

# Get randomly 5 SNPs
id.snps <- sample(x=seq(1, nbr.snps), size=5)

name.snps <- bim.data$V2[id.snps]

character.snps <- paste0(unlist(name.snps),collapse=",")

# Cluster Script path
sbatch.file <- "/home/duvanel/ClusterScript/ExtractSNPS.sbatch"

# Write into a sbatch file
sink(sbatch.file)
cat("#!/bin/bash
    
    #set a job name
    #SBATCH --job-name=plinkExtractSNP
    
    #a file for job output, you can check job progress
    #SBATCH --output=plinkExtractSNP.out
    
    # a file for errors from the job
    #SBATCH --error=plinkExtractSNP.err
    
    #time you think you need; default is one hour
    #in minutes in this case
    #SBATCH --time=600:00
    
    #quality of service; think of it as job priority
    #SBATCH --qos=normal
    
    #number of nodes you are requesting
    #SBATCH --nodes=1
    
    #memory per node; default is 4000 MB per CPU
    #SBATCH --mem=4000
    
    #get emailed about job BEGIN, END, and FAIL 
    #SBATCH --mail-type=ALL
    
    #who to send email to; please change to your email
    #SBATCH  --mail-user=duvanel@stanford.edu
    
    #task to run per node; each node has 16 cores
    #SBATCH --ntasks-per-node=16
    
    module load R
    plink --bfile", files," --snps ", character.snps)
sink()

# Give chmod +x
system(paste0("chmod +x ", sbatch.file))

# Execute using cluster
#system(paste0("sbatch ", sbatch.file))
