#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import sys
import time
import datetime

from  numpy import pi as pi
from numpy import shape
import h5py
import getopt
import random
import matplotlib.pyplot as plt


def Usage():
  print '''
    ****************************************************************************************
    ****************************************************************************************
    Simulating a set of interferograms based on real interferograms and an existing displacement velocity field.

    simulation.py -i interferograms -v velocity -u unwrapping_error_percentage -m mask for unwrapping error addition -y susbet in y -x subste in x -o output name
    ****************************************************************************************
    ****************************************************************************************
'''

def main(argv):
 
  unwrapMaskFile='Mask.h5' 
  unw_error_perc= 0.0
  outName= 'simulatedIgrams.h5'
  try:
      opts, args = getopt.getopt(argv,"h:i:v:u:m:o:x:y:")

  except getopt.GetoptError:
      Usage() ; sys.exit(1)

  for opt,arg in opts:
      if opt in ("-h","--help"):
        Usage()
        sys.exit()
      elif opt == '-i':
        igramFile = arg
      elif opt == '-v':
        velocityFile = arg
      elif opt == '-u':
        unw_error_perc = float(arg)  #percentage of interferograms with unwrapping error
      elif opt == '-o':
        outName = arg
      elif opt == '-m':
        unwrapMaskFile=arg
      elif opt == '-x':
        xsub = [int(i) for i in arg.split(':')]
        xsub.sort()
      elif opt == '-y':
        ysub = [int(i) for i in arg.split(':')]
        ysub.sort()

  try:
     h5file=h5py.File(igramFile,'r')
     h5vel=h5py.File(velocityFile,'r')
     h5mask=h5py.File(unwrapMaskFile,'r')
  except:
     Usage() ; sys.exit(1)


  kv=h5vel.keys()
  vset=h5vel[kv[0]].get(kv[0])
  rate=vset[0:vset.shape[0],0:vset.shape[1]]
#####################################################
  igramList=h5file['interferograms'].keys()
  km=h5mask.keys()
  mset=h5mask[km[0]].get(km[0])
  mask=mset[0:mset.shape[0],0:mset.shape[1]]
  indx=range(len(igramList))
  if unw_error_perc >=1.0:
      unw_error_perc = unw_error_perc/100.0
  num_unw_err_igrams=int(unw_error_perc*len(igramList))
  indx_unw_error=random.sample(indx,num_unw_err_igrams)
  unw_err_list=[igramList[i] for i in indx_unw_error]
  unwrapError=mask*2.0*pi

  print unw_err_list
#####################################################
  try:
    ysub
  except:
    ysub=[0,mset.shape[0]]
  
  try:
    xsub
  except:
    xsub=[0,mset.shape[1]]  
#####################################################
  h5simulate=h5py.File(outName,'w') 
  gg = h5simulate.create_group('interferograms')
  wvl=float(h5file['interferograms'][igramList[0]].attrs['WAVELENGTH'])
  print wvl
  for igram in igramList:
      print igram
      date = h5file['interferograms'][igram].attrs['DATE12'].split('-')
      if date[0][0] == '9':
        date[0] = '19'+date[0]
      else:
        date[0] = '20'+date[0]
      if date[1][0] == '9':
        date[1] = '19'+date[1]
      else:
        date[1] = '20'+date[1]
      t1=datetime.datetime(*time.strptime(date[0],"%Y%m%d")[0:5])
      t2=datetime.datetime(*time.strptime(date[1],"%Y%m%d")[0:5])
      dt=(t2-t1)
      dt=float(dt.days)/365.0 
      group = gg.create_group(igram)
      unw=-rate*dt*4.0*pi/wvl
      if igram in unw_err_list:
         print 'adding unwrapping error to: ' + igram
         unw=unw+unwrapError

      
      dset = group.create_dataset(igram, data=unw[ysub[0]:ysub[1],xsub[0]:xsub[1]], compression='gzip')
      for key, value in h5file['interferograms'][igram].attrs.iteritems():
          group.attrs[key] = value
      if igram in unw_err_list:
          group.attrs['unwrap_error'] = 'yes'
      else:
          group.attrs['unwrap_error'] = 'no'

      group.attrs['FILE_LENGTH']=shape(unw[ysub[0]:ysub[1],xsub[0]:xsub[1]])[0]
      group.attrs['WIDTH']=shape(unw[ysub[0]:ysub[1],xsub[0]:xsub[1]])[1]

  h5MASK=h5py.File('Mask.h5','r')
  MASK=h5MASK[h5MASK.keys()[0]].get(h5MASK.keys()[0])
      
  h5MASKSim=h5py.File('simulatedMask.h5','w')
  ggg=h5MASKSim.create_group('mask')
  ggg.create_dataset('mask',data=MASK[ysub[0]:ysub[1],xsub[0]:xsub[1]])
 
  gm=h5simulate.create_group('mask')
  gm.create_dataset('mask',data=MASK[ysub[0]:ysub[1],xsub[0]:xsub[1]])

  h5file.close()
  h5MASKSim.close()
  h5MASK.close()

  

if __name__ == '__main__':
                  
  main(sys.argv[1:]) 


