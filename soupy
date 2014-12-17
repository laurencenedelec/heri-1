#!/usr/bin/python

import os
import datetime 

howmany = 240

startjob = datetime.datetime.now() 

for i in xrange(1, howmany+1): 
  num = i/float( howmany)
  job = 'qsub -N poviogen-%03d ~/poviogen.submit&nbsp;%.3f&nbsp;%03d'&nbsp;% (i, num, i)
  os.system(job)
  print job

endjob = datetime.datetime.now()

print 'time to submit all jobs: ', endjob-startjob