#!/bin/bash                                                                    
#SBATCH --mail-user=nedelec@math.stanford.edu                                                                                             
#SBATCH --mail-type=ALL                                                                                                                   
#SBATCH --qos=normal                                                                                                                      
#SBATCH --time=1:00:00                                                                                                                    
#SBATCH --job-name=desi                                                                                                           
#SBATCH --cpus-per-task=16                                                                                                             
#SBATCH --mem=64000                                                                                                                
#SBATCH --output=desi.out                                                                                                         
#SBATCH --error=desi.err                                                                                                          
module load R                                                                                                                             
srun -n 1 Rscript ~/NFG/ler/desi.R  
