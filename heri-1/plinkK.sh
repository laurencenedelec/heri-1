#!/bin/bash                                                                                       
#SBATCH --mail-user=nedelec@math.stanford.edu                                                   
#SBATCH --mail-type=ALL                                                                           
#SBATCH --qos=normal                                                                              
#SBATCH --job-name=plinkKrun                                                                      
#SBATCH --time=2:00:00                                                                           
#SBATCH --output=plinkK.out
#SBATCH --error=plinkK.err                                                                         
#SBATCH --cpus-per-task=16                                                                        
#SBATCH --mem=64000  

module load plink/1.9_b1
#module load plink

##compute the simmilarity
#plink --bfile ~/NFG/NFBC/Mat/NFBC_dbGaP_20091127 --noweb  --matrix --out ~/NFG/result/NFGsimi

##compute the distance
## plink --bfile ~/NFG/NFBC/Mat/NFBC_dbGaP_20091127 --noweb --cluster --distance-matrix --out ~/NFG/result/plinkdis

## Start Time       = 05/14/2014 12:33:25
## End Time         = 05/17/2014 15:33:28
## CPU              = 3:02:57:08
## Max vmem         = 1.077G

## compute the ibd from ibs and ibs 
 #plink --noweb --bfile ~/NFG/NFBC/Mat/NFBC_dbGaP_20091127 --genome --out ~/NFG/result/NFGibsibd

 #gawk ' { print $2,$4,$10,$12 } ' ~/NFG/result/NFGibsibd.genome > ~/NFG/result/NFGredibsibd.genome

##Start Time       = 05/14/2014 12:34:40
##End Time         = 05/17/2014 15:54:00
## CPU              = 3:03:15:22
## Max vmem         = 1019.473M

## compute the frequency of the SNP
#plink --bfile ~/NFG/NFBC/Mat/NFBC_dbGaP_20091127 --noweb --freq  --maf 0.01 --hwe 0.0001 --geno 0.05 --out  ~/NFG/result/freqNFG

##compute the distance with the frequency as weight for plink 1.9
##problem avec le retour chariot
 cd ~/NFG
#plink --bfile NFBC/Mat/NFBC_dbGaP_20091127 --maf 0.01 --hwe 0.0001 --geno 0.05 --distance-exp 0.5 --distance ibs --read-freq result/freqNFG.frq --out result/plinkdisex

##plink --bfile NFBC/Mat/NFBC_dbGaP_20091127 --maf 0.01  --hwe 0.0001 --geno 0.05 --distance-exp 0,5 --distance ibs --out result/plinkdisex

plink --bfile NFBC/Mat/NFBC_dbGaP_20091127 --maf 0.01  --hwe 0.0001 --geno 0.05 --make-rel square --out result/plinkgcta                 
