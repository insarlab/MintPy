#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import sys
import os
try:
  files=sys.argv[1].split(',')
  print files
except:
  print '''

  usage: Matching_all.py 'file1.h5,file2.h5,...,filen.h5'

'''  
  sys.exit(1)

f1=files[0]
f2=files[1]

cmd='Matching.py '+f1+' '+f2
print cmd
os.system(cmd)
f1=f1.split('.h5')[0]+'_'+f2.split('.h5')[0]+'.h5'
cmd='mv Matched.h5 '+f1
os.system(cmd)

if len(files)>2:
  for i in range(2,len(files)):
    print '''
    ---------------------------------

    '''
    f2=files[i]
    cmd='Matching.py '+f1+' '+f2
    os.system(cmd)
    print cmd
    f1=f1.split('.h5')[0]+'_'+f2.split('.h5')[0]+'.h5'
    cmd='mv Matched.h5 '+f1
    os.system(cmd)


