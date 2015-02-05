#!/bin/bash                                                                                      
#SBATCH --mail-user=nedelec@math.stanford.edu                                                    
#SBATCH --mail-type=ALL                                                                          
#SBATCH --qos=normal                                                                             
#SBATCH --job-name=heri                                                                       
#SBATCH --time=6:30:00                                                                           
#SBATCH --output=heri3.out                                                                        
#SBATCH --error=heri3.err                                                                          
#SBATCH --cpus-per-task=16                                                                       
#SBATCH --mem=64000                                                                               
module load R                                                                                    
srun -n 1 Rscript ~/NFG/ler/heritabilite.R                                                                
