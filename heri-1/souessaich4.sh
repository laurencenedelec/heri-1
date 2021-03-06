#!/bin/ksh                                                                      
    SHFILE=essai15.sh
     cat > $SHFILE <<EOF                                                                                                    
#!/bin/bash
#SBATCH --mail-user=nedelec@math.stanford.edu   
#SBATCH --mail-type=ALL
#SBATCH --qos=normal
#SBATCH --time=48:00:00   
#SBATCH --job-name=essai15
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
#SBATCH --output=essai15.out
#SBATCH --error=essai15.err
module load R
srun -n 1 Rscript ~/NFG/ler/essaich4.R 
EOF
chmod +x $SHFILE                                                                                                               
sbatch $SHFILE
                                                                                                         
 
                                                                                                          
