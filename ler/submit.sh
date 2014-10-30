## qsub -l mem_free=5G ~/bin/rjob.sh ~/NFG/ler/ppermdcor.R 1
## qsub -l mem_free=5G ~/bin/rjob.sh ~/NFG/ler/ppermdcor.R 2
## qsub -l mem_free=5G ~/bin/r2job.sh  ~/NFG/ler/ppermdcor.R 3
## qsub -l mem_free=5G ~/bin/r2job.sh ~/NFG/ler/ppermdcor.R 4
## qsub -l mem_free=5G ~/bin/r2job.sh ~/NFG/ler/ppermdcor.R 5
## qsub -l mem_free=5G ~/bin/rjob.sh ~/NFG/ler/ppermdcor.R 6
## qsub -l mem_free=5G ~/bin/rjob.sh  ~/NFG/ler/ppermdcor.R 7
## qsub -l mem_free=5G ~/bin/rjob.sh ~/NFG/ler/ppermdcor.R  8
## qsub -l mem_free=5G ~/bin/rjob.sh ~/NFG/ler/ppermdcor.R 9
##
## qsub -l mem_free=5G ~/bin/rjob.sh ~/NFG/ler/stand.R 1
## qsub -l mem_free=5G ~/bin/rjob.sh ~/NFG/ler/stand.R 2
## qsub -l mem_free=5G ~/bin/rjob.sh  ~/NFG/ler/stand.R 3
## qsub -l mem_free=5G ~/bin/rjob.sh ~/NFG/ler/stand.R 10                                                                               
## qsub -l mem_free=5G ~/bin/rjob.sh ~/NFG/ler/stand.R 11                                                                              
## qsub -l mem_free=5G ~/bin/rjob.sh ~/NFG/ler/stand.R 12   
##
##qsub -l mem_free=5G ~/bin/rjob.sh ~/NFG/ler/dcovpppheno4.R  11
##qsub -l mem_free=5G ~/bin/r2job.sh ~/NFG/ler/dcorpartchr.R 1
##qsub -l mem_free=5G ~/bin/r2job.sh ~/NFG/ler/dcorpartchr.R 2                                                     
##qsub -l mem_free=5G ~/bin/r2job.sh  ~/NFG/ler/dcorpartchr.R 3
##qsub -l mem_free=5G ~/bin/r2job.sh ~/NFG/ler/essai.R 1
## qsub -l mem_free=30G ~/bin/r2job.sh ~/NFG/ler/essai.R 2
##qsub -l mem_free=5G ~/bin/r2job.sh  ~/NFG/ler/essai.R 3
##qsub -l mem_free=5G ~/bin/r2job.sh ~/NFG/ler/essai.R 4
##qsub -l mem_free=5G ~/bin/r2job.sh ~/NFG/ler/essai.R 5
##qsub -l mem_free=5G ~/bin/r2job.sh  ~/NFG/ler/essai.R 6
##qsub -l mem_free=5G ~/bin/r2job.sh ~/NFG/ler/essai.R 7
##qsub -l mem_free=5G ~/bin/r2job.sh ~/NFG/ler/essai.R 8
##qsub -l mem_free=5G ~/bin/r2job.sh  ~/NFG/ler/essai.R 9
## qsub -l mem_free=5G ~/bin/rjob.sh ~/NFG/ler/dcorpartchr.R 4                                                    
## qsub -l mem_free=5G ~/bin/rjob.sh ~/NFG/ler/dcorpartchr.R 5
## qsub -l mem_free=5G ~/bin/rjob.sh ~/NFG/ler/ dcorpartchr.R 6                                                   
## qsub -l mem_free=5G ~/bin/rjob.sh  ~/NFG/ler/dcorpartchr.R 7                                                   
## qsub -l mem_free=5G ~/bin/rjob.sh ~/NFG/ler/dcorpartchr.R  8                                                   
##  qsub -l mem_free=5G ~/bin/rjob.sh ~/NFG/ler/dcorpartchr.R 9
qsub -l mem_free=3G ~/bin/r2job.sh ~/NFG/ler/select.R 4
##qsub -l mem_free=3G ~/bin/r2job.sh ~/NFG/ler/select.R 5
##qsub -l mem_free=3G ~/bin/r2job.sh  ~/NFG/ler/select.R 6
srun -n 1 Rscript ~/bin/rjob.sh select.R 4



## time it take to run any of these dcorpartchr

## Start Time       = 05/12/2014 11:19:29
## End Time         = 05/15/2014 20:25:22
## CPU              = 3:08:54:44
## Max vmem         = 3.720G


## Wallclock Time   = 3:09:25:59
## CPU              = 3:09:14:57
## Max vmem         = 3.758G

## time it take if only on pheno is select                                                                          

## Start Time       = 05/09/2014 11:36:36                                                                           
# End Time         = 05/11/2014 00:10:45                                                                            
## CPU              = 1:12:33:36                                                                                    
## Max vmem         = 2.420G     