#!/bin/ksh
for posi in  4 5 6  ; do
    SHFILE=essaipdcorm-$posi.sh
     cat > $SHFILE <<EOF
#!/bin/bash
#SBATCH --mail-user=nedelec@math.stanford.edu   
#SBATCH --mail-type=ALL
#SBATCH --qos=normal
#SBATCH --job-name=essaipdcorm-$posi
#SBATCH --time=48:00:00
#SBATCH --output=essaimibs-$posi.out
#SBATCH --error=essaipdcorm-$posi.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000
#################
module load R
srun -n 1 Rscript ~/NFG/ler/essaipdcor.R $posi
EOF
chmod +x $SHFILE
sbatch $SHFILE
done



##srun -n 1 Rscript $@

