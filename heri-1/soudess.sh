#!/bin/ksh
for posi in 1 4 7 ; do
    SHFILE=dessin-$posi.sh
    cpus=1
    mem=4800
    cat > $SHFILE <<EOF
#!/bin/bash
#SBATCH --mail-user=nedelec@math.stanford.edu   
#SBATCH --mail-type=ALL
#SBATCH --qos=normal
#SBATCH --time=1:00:00
#SBATCH --job-name=dessin-$posi
#SBATCH --cpus-per-task=$cpus
#SBATCH --mem=$mem
#SBATCH --output=dessin-$posi.out
#SBATCH --error=dessin-$posi.err
module load R
srun -n 1 Rscript ~/NFG/ler/dessinpv.R $posi
EOF

chmod +x $SHFILE
sbatch $SHFILE
done

##srun -n 1 Rscript $@



