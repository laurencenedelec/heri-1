#!/bin/bash                                                                                                                                          
#$ -cwd                                                                                                                                              
#$ -S /bin/bash                                                                                                                                      
#$ -j y                                                                                                                                              

# email address to send notices to                                                                                                                   
#$ -M nedelec@math.stanford.edu                                                                                                                      
#$-m bea                                                                                                                                             

## set a name to make it easy to pick out of qstat output                                                                                            
#$ -N finrun                                                                                                                                         

        module load plink

##  plink --bfile ~/NFG/raw/NFBC/Mat/NFBC_dbGaP_20091127 --noweb --cluster --matrix                                                                    
## Clean genotype data using criteria that p-val for HWE < 0.0001,                                                                                   
## MAF > 0.01 and                                                                                                                                    
## call rate > 95%                                                                                                 
##extract part of chr 1,4,19 for testing dcor                                  

plink --noweb --from-kb 65000 --to-kb 83000 --chr 4 --bfile ~/NFG/raw/NFBC/Mat/NFBC_dbGaP_20091127  --maf 0.01 --hwe 0.0001 --geno 0.05 --make-bed --out ~/NFG/raw/Hd4cl
plink --noweb --from-kb 40000 --to-kb 60000 --chr 19 --bfile ~/NFG/raw/NFBC/Mat/NFBC_dbGaP_20091127   --maf 0.01 --hwe 0.0001 --geno 0.05 --make-bed  --out ~/NFG/raw/Hd19cl 
plink --noweb --from-kb 47000 --to-kb 63000 --chr 1 --bfile ~/NFG/raw/NFBC/Mat/NFBC_dbGaP_20091127 --maf 0.01 --hwe 0.0001 --geno 0.05 --make-bed  --out ~/NFG/raw/Hd1cl  


## Additional step to filter out copy number variants                                                                                                
plink --noweb --bfile ~/NFG/raw/Hd19cl --write-snplist --out ~/NFG/raw/Hd19cl
grep "cnv" ~/NFG/raw/Hd19cl.snplist > cnvs.txt
plink --noweb --bfile ~/NFG/raw/Hd19cl --exclude cnvs.txt --make-bed  --out ~/NFG/raw/Hd19cl
plink --noweb --bfile ~/NFG/raw/Hd19cl --recodeA --out ~/NFG/raw/Hd19clA

plink --noweb --bfile ~/NFG/raw/Hd4cl --write-snplist --out ~/NFG/raw/Hd4cl
grep "cnv" ~/NFG/raw/Hd4cl.snplist > cnvs.txt
plink --noweb --bfile ~/NFG/raw/Hd4cl --exclude cnvs.txt --make-bed  --out ~/NFG/raw/Hd4cl
plink --noweb --bfile ~/NFG/raw/Hd4cl --recodeA	 --out ~/NFG/raw/Hd4clA

plink --noweb --bfile ~/NFG/raw/Hd1cl --write-snplist --out ~/NFG/raw/Hd1cl
grep "cnv" ~/NFG/raw/Hd1cl.snplist > cnvs.txt
plink --noweb --bfile ~/NFG/raw/Hd1cl --exclude cnvs.txt --make-bed  --out ~/NFG/raw/Hd1cl
plink --noweb --bfile ~/NFG/raw/Hd1cl --recodeA	 --out ~/NFG/raw/Hd1clA

##For heritability  full size SNP
plink --noweb  --bfile ~/NFG/raw/NFBC/Mat/NFBC_dbGaP_20091127 --maf 0.01 --hwe 0.0001 --geno 0.05 --make-bed  --out ~/NFG/raw/allgene
plink --noweb --bfile ~/NFG/raw/allgene --write-snplist --out ~/NFG/raw/allgene
plink --noweb --bfile ~/NFG/raw/allgene --recodeA  --out ~/NFG/raw/allgeneA

##For comparing with standart p-value full size chr4
plink --noweb  --chr 4 --bfile ~/NFG/raw/NFBC/Mat/NFBC_dbGaP_20091127  --maf 0.01 --hwe 0.0001 --geno 0.05 --make-bed --out ~/NFG/raw/Hd4fullcl
plink --noweb --bfile ~/NFG/raw/Hd4fullcl --write-snplist --out ~/NFG/raw/Hd4fullcl
grep "cnv" ~/NFG/raw/Hd4fullcl.snplist > cnvs.txt
plink --noweb --bfile ~/NFG/raw/Hd4fullcl --exclude cnvs.txt --make-bed  --out ~/NFG/raw/Hd4fullcl
plink --noweb --bfile ~/NFG/raw/Hd4fullcl --recodeA  --out ~/NFG/raw/Hd4fullclA

##the full 19 
plink --noweb  --chr 19 --bfile ~/NFG/raw/NFBC/Mat/NFBC_dbGaP_20091127  --maf 0.01 --hwe 0.0001 --geno 0.05 --make-bed --out ~/NFG/raw/Hd19fullcl
plink --noweb --bfile ~/NFG/raw/Hd19fullcl --write-snplist --out ~/NFG/raw/Hd19fullcl
grep "cnv" ~/NFG/raw/Hd19fullcl.snplist > cnvs.txt
plink --noweb --bfile ~/NFG/raw/Hd19fullcl --exclude cnvs.txt --make-bed  --out ~/NFG/raw/Hd19fullcl
plink --noweb --bfile ~/NFG/raw/Hd19fullcl --recodeA  --out ~/NFG/raw/Hd19fullclA


##the full 1
plink --noweb --chr 1 --bfile ~/NFG/raw/NFBC/Mat/NFBC_dbGaP_20091127 --maf 0.01 --hwe 0.0001 --geno 0.05 --make-b\
ed --out ~/NFG/raw/Hd1fullcl
plink --noweb --bfile ~/NFG/raw/Hd1fullcl --write-snplist --out ~/NFG/raw/Hd1fullcl
grep "cnv" ~/NFG/raw/Hd1fullcl.snplist > cnvs.txt
plink --noweb --bfile ~/NFG/raw/Hd1fullcl --exclude cnvs.txt --make-bed  --out ~/NFG/raw/Hd1fullcl
plink --noweb --bfile ~/NFG/raw/Hd1fullcl --recodeA  --out ~/NFG/raw/Hd1fullclA

