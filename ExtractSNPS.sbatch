#!/bin/bash
        
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
        
        module load plink
        plink --bfile /home/duvanel/BP/Data/GenoData1M/MendelAnalysis3  --snps  kgp22815272,rs7816543,rs2838128,kgp11032426,rs701089  --recodeA --out  /home/duvanel/git/PackageCH/data/SNPs_14102014_104636 
