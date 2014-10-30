#!/bin/ksh
for posi in  4 5 7 8 9    ; do
    SHFILE=select-$posi.sh
    cpus=16
    mem=64000
     cat > $SHFILE <<EOF
#!/bin/bash
#SBATCH --mail-user=nedelec@math.stanford.edu   
#SBATCH --mail-type=ALL
#SBATCH --qos=normal
#SBATCH --job-name=select-$posi
#SBATCH --time=48:00:00
#SBATCH --output=select-$posi.out
#SBATCH --error=select-$posi.err
#SBATCH --mem=$mem
#SBATCH --cpus-per-task=$cpus   
module load R
srun -n 1 Rscript ~/NFG/ler/select.R $posi
EOF

chmod +x $SHFILE
sbatch $SHFILE
done



##srun -n 1 Rscript $@

