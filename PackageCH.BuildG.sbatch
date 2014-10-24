#!/bin/bash 

#set a job name  
#SBATCH --job-name=PackageCH.BuildG

#a file for job output, you can check job progress
#SBATCH --output=PackageCH.BuildG.out

# a file for errors from the job
#SBATCH --error=PackageCH.BuildG.err

#time you think you need; default is one hour
#in minutes in this case
#SBATCH --time=120:00

#quality of service; think of it as job priority
#SBATCH --qos=normal

#number of nodes you are requesting
#SBATCH --nodes=1

#memory per node; default is 8000 MB per CPU
#SBATCH --mem=8000

#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=ALL

#who to send email to; please change to your email
#SBATCH  --mail-user=duvanel@stanford.edu

#task to run per node; each node has 16 cores
#SBATCH --ntasks-per-node=4

module load R
srun -n 1 Rscript ~/git/PackageCH/Main_2.R
