#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import sys
import os
import numpy as np
import h5py
import getopt
import _pysar_utilities as ut
from scipy.linalg import pinv as pinv


def Usage():
  print '''
  Finding and correcting the unwrapping errors based on the triangular consistency of interferograms (phase closure)

  Usage: 

    unwrap_error.py Seeded_LoadedData_BajaT499EnvD2.h5 mask.h5

    unwrap_error.py Seeded_LoadedData_BajaT499EnvD2.h5

  mask.h5 : is a mask file to specify tgose pixels which user wants to correct for unwrapping errors.
  

'''

def main(argv):
   try:
      file=argv[0]
   except:
      Usage();sys.exit(1)

   h5file=h5py.File(file)
   ifgramList = h5file['interferograms'].keys()
   sx = int(h5file['interferograms'][ifgramList[0]].attrs['WIDTH'])
   sy = int(h5file['interferograms'][ifgramList[0]].attrs['FILE_LENGTH'])
   curls,Triangles,C=ut.get_triangles(h5file)
   A,B = ut.design_matrix(h5file)   
   ligram,lv=np.shape(B)
   lcurls=np.shape(curls)[0]
   print 'Number of all triangles: '+  str(lcurls)
   print 'Number of interferograms: '+ str(ligram)
#   print curls

   curlfile='curls.h5'
   if not os.path.isfile(curlfile):
       ut.generate_curls(curlfile,h5file,Triangles,curls)
  
    
   thr=0.50
   curls=np.array(curls)

   n1=curls[:,0]
   n2=curls[:,1]
   n3=curls[:,2]

   numPixels=sy*sx
   print 'reading interferograms...'   
   data = np.zeros((ligram,numPixels),np.float32)
   for ni in range(ligram):
      dset=h5file['interferograms'][ifgramList[ni]].get(ifgramList[ni])
      d = dset[0:dset.shape[0],0:dset.shape[1]]
      data[ni] = d.flatten(1)   

   print np.shape(data)

   print 'reading curls ...' 
   h5curl=h5py.File(curlfile)
   curlList=h5curl['interferograms'].keys()
   curlData = np.zeros((lcurls,numPixels),np.float32)
   for ni in range(lcurls):
      dset=h5curl['interferograms'][curlList[ni]].get(curlList[ni])
      d = dset[0:dset.shape[0],0:dset.shape[1]]
      curlData[ni] = d.flatten(1)
   pi=np.pi
   EstUnwrap=np.zeros((ligram,numPixels),np.float32)


   try:
      maskFile=argv[1]
      h5Mask=h5py.File(maskFile)
      dset = h5Mask['mask'].get('mask')
      Mask=dset[0:dset.shape[0],0:dset.shape[1]]
   except:
      dset = h5file['mask'].get('mask')
      Mask=dset[0:dset.shape[0],0:dset.shape[1]]
   
   Mask=Mask.flatten(1)

   for ni in range(numPixels):
      #dU = np.zeros([ligram,1])
      #print np.shape(dU)
      #print np.shape(data[:,ni])
    if Mask[ni]==1:

      dU = data[:,ni]
    #  nan_ndx = dataPoint == 0.
      unwCurl = np.array(curlData[:,ni])
     # print unwCurl

      ind=np.abs(unwCurl)>=thr
      N1=n1[ind]
      N2=n2[ind]
      N3=n3[ind]
      
      indC=np.abs(unwCurl)<thr
      Nc1=n1[indC]
      Nc2=n2[indC]
      Nc3=n3[indC]

      N=np.hstack([N1,N2,N3])
      UniN=np.unique(N)

      Nc=np.hstack([Nc1,Nc2,Nc3])
      
      UniNc=np.unique(Nc)

      inter=list(set(UniNc) & set(UniN)) # intersetion
      UniNc= list(UniNc)
      for x in inter:
         UniNc.remove(x)
      
      D=np.zeros([len(UniNc),ligram])
      for i in range(len(UniNc)):
         D[i,UniNc[i]]=1
   
      AAA=np.vstack([-2*pi*C,D])
     # AAA1=np.hstack([AAA,np.zeros([AAA.shape[0],lv])])
     # AAA2=np.hstack([-2*pi*np.eye(ligram),B]) 
     # AAAA=np.vstack([AAA1,AAA2])
      AAAA=np.vstack([AAA,0.25*np.eye(ligram)])

     # print '************************'
     # print np.linalg.matrix_rank(C)
     # print np.linalg.matrix_rank(AAA) 
     # print np.linalg.matrix_rank(AAAA)
     # print '************************'

     # LLL=list(np.dot(C,dU)) + list(np.zeros(np.shape(UniNc)[0]))# + list(dU)
     # ind=np.isnan(AAA)
     # M1=pinv(AAA)      
     # M=np.dot(M1,LLL)
     # EstUnwrap[:,ni]=np.round(M[0:ligram])*2.0*np.pi

      ##########
      # with Tikhonov regularization:
      AAAA=np.vstack([AAA,0.25*np.eye(ligram)])
      LLL=list(np.dot(C,dU)) + list(np.zeros(np.shape(UniNc)[0])) + list(np.zeros(ligram))
      ind=np.isnan(AAAA)
      M1=pinv(AAAA)
      M=np.dot(M1,LLL)
      EstUnwrap[:,ni]=np.round(M[0:ligram])*2.0*np.pi
    #  print M[0:ligram]
    #  print np.round(M[0:ligram])

    else:
      EstUnwrap[:,ni]=np.zeros([ligram])
      if not np.remainder(ni,10000): print 'Processing point: %7d of %7d ' % (ni,numPixels)

   dataCor=data+EstUnwrap
   unwCorFile=file.replace('.h5','')+'_unwCor.h5'
   h5unwCor=h5py.File(unwCorFile,'w') 
   gg = h5unwCor.create_group('interferograms') 
   for i in range(ligram):
       group = gg.create_group(ifgramList[i])
      # dset = group.create_dataset(ifgramList[i], data=np.reshape(dataCor[i,:],[sy,sx]), compression='gzip')
       dset = group.create_dataset(ifgramList[i], data=np.reshape(dataCor[i,:],[sx,sy]).T, compression='gzip')
      
       for key, value in h5file['interferograms'][ifgramList[i]].attrs.iteritems():
          group.attrs[key] = value

   MASK=h5file['mask'].get('mask')
   gm = h5unwCor.create_group('mask')
   dset = gm.create_dataset('mask', data=MASK, compression='gzip')

   try:
      Coh=h5file['meanCoherence'].get('meanCoherence')
      gc = h5unwCor.create_group('meanCoherence')
      dset = gc.create_dataset('meanCoherence', data=Coh, compression='gzip')
   except:
      print ''

   h5unwCor.close()
   h5file.close()
   h5curl.close() 

if __name__ == '__main__':

  main(sys.argv[1:])

