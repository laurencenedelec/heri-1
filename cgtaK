#!/bin/bash                                                                                       
#SBATCH --mail-user=nedelec@math.stanford.edu                                                   
#SBATCH --mail-type=ALL                                                                           
#SBATCH --qos=normal                                                                              
#SBATCH --job-name=cgtaKrun                                                                      
#SBATCH --time=2:00:00                                                                           
#SBATCH --output=cgta.out
#SBATCH --error=cgta.err                                                                         
#SBATCH --cpus-per-task=16                                                                        
#SBATCH --mem=64000  



module load gcta

gcta64 --bfile ~/NFG/raw/NFBC/Mat/NFBC_dbGaP_20091127 --autosome --maf 0.01 --make-grm --out ~/NFG/raw/Kcgta --thread-num 10


#cgtaKrun Ended, Run time 00:01:56