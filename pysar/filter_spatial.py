#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import sys
import os
import numpy as np
import getopt
import h5py
import _readfile as readfile
import _writefile as writefile

try:
    from skimage.filter import roberts, sobel,canny,gaussian_filter
    
except:
    print '++++++++++++++++++++++++++++++++++++++++++++'
    print ''
    print 'Could not import skimage'
    print 'To use filter.py you must install skimage'
    print 'See: http://scikit-image.org/'
    print ''
    print '++++++++++++++++++++++++++++++++++++++++++++'
    sys.exit(1)

def Usage():

   print '''

   ***************************************************************
   Spatial filtering of the time-series or velocity ...

   Usage:

    filter.py -f  file -t filter_type  -p parameters

    file: PySAR h5 files [interferogram, time-series, velocity], 
          roi_pac files [dem, unw, hgt, ]
          or image files [jpeg, jpg, png]

    filter_type:
                  lowpass_gaussian [-p defines the sigma] 
                  highpass_gaussian [-p defines the sigma]
                  lowpass_avg  [-p defines the size of the kernel]
                  highpass_avg  [-p defines the size of the kernel]
                  sobel
                  roberts
                  canny              

   Example:
             filter.py -f timeseries.h5 -t lowpass_avg -p 5  
             filter.py -f velocity.h5 -t lowpass_avg -p 5             
             filter.py -f velocity.h5 -t sobel
             filter.py -f velocity.h5 -t highpass_gaussian -p 3

   ***************************************************************

'''

def filter(data,filtType,par):

    if filtType == "sobel":
       filt_data = sobel(data)
    elif filtType == "roberts":
       filt_data = roberts(data)
    elif filtType ==  "canny":
       filt_data = canny(data)
    elif filtType ==  "lowpass_avg":
       from scipy import ndimage
       p=int(par)
       kernel = np.ones((p,p),np.float32)/(p*p)
       filt_data = ndimage.convolve(data, kernel)
    elif filtType ==  "lowpass_gaussian":
        
       s=float(par)
       filt_data = gaussian_filter(data, sigma=s)

    elif filtType ==  "highpass_gaussian":

       s=float(par)
       lp_data = gaussian_filter(data, sigma=s)
       filt_data = data - lp_data

    elif filtType ==  "highpass_avg":
       from scipy import ndimage
       p=int(par)
       kernel = np.ones((p,p),np.float32)/(p*p)
       lp_data = ndimage.convolve(data, kernel)
       filt_data = data - lp_data

    #elif filtType ==  "gradient":
       
    return filt_data

def multilook(ifg,lksy,lksx):

    
    rows,cols=ifg.shape
    rows_lowres=np.floor(rows/lksy)
    cols_lowres=np.floor(cols/lksx)


   # %ifg_lowres=NaN(rows_lowres,cols_lowres)
    
    thr = np.floor(lksx*lksy/2)
    
    
    ifg_Clowres=np.zeros((rows,cols_lowres))
    ifg_lowres=np.zeros((rows_lowres,cols_lowres))
    
    
    for c in range(int(cols_lowres)):
        ifg_Clowres[:,c]=np.sum(ifg[:,(c)*lksx:(c+1)*lksx],1)
        
    
    
    for r in range(int(rows_lowres)):
        
        ifg_lowres[r,:]=np.sum(ifg_Clowres[(r)*lksy:(r+1)*lksy,:],0)
    
    
    
    ifg_lowres=ifg_lowres/(lksy*lksx) 
    return ifg_lowres

def main(argv):

  try:
    opts, args = getopt.getopt(argv,"h:f:t:p:")
    
  except getopt.GetoptError:
    Usage() ; sys.exit(1)

  if opts==[]:
    Usage() ; sys.exit(1)
  for opt,arg in opts:
    if opt in ("-h","--help"):
      Usage()
      sys.exit()
    elif opt == '-f':
      file = arg
    elif opt == '-t':
      filtType=arg
    elif opt == '-p':
      par=arg


#  try:  
#    file=argv[0]
#    alks=float(argv[1])
#    rlks=float(argv[2])
#  except:
#    Usage();sys.exit(1)



  ext = os.path.splitext(file)[1]
  outName=file.split('.')[0]+'_'+filtType+ext
  try:
    par
  except:
    par=[]

  print '+++++++++++++++++++++++++++'
  print 'Filter type : '+filtType
  print 'parameters : ' + str(par)
  print '+++++++++++++++++++++++++++'
###############################################
  if ext == '.int' or ext == '.slc':
    a,p,r = readfile.read_complex64(file)
    plks=multilook(p,alks,rlks)
    alks=multilook(a,alks,rlks)


    r['FILE_LENGTH']=str(dlks.shape[0])
    r['WIDTH']=str(dlks.shape[1])
    r['XMAX']=str(int(r['WIDTH']) - 1)
    r['YMAX']=str(int(r['FILE_LENGTH']) - 1)
    try:
       r['Y_STEP']=str(float(r['Y_STEP'])*alks)
       r['X_STEP']=str(float(r['X_STEP'])*rlks)
    except:
       Geo=0

    f = open(outName+'.rsc','w')
    for k in r.keys():
       f.write(k+'    '+r[k]+'\n')
    f.close()   

  elif ext == '.unw' or ext == '.cor' or ext == '.hgt':
    a,p,r = readfile.read_float32(file)
    plks=multilook(p,alks,rlks)
    alks=multilook(a,alks,rlks)
    
    writefile.write_float32(plks,outName)

    r['FILE_LENGTH']=str(dlks.shape[0])
    r['WIDTH']=str(dlks.shape[1])
    r['XMAX']=str(int(r['WIDTH']) - 1)
    r['YMAX']=str(int(r['FILE_LENGTH']) - 1)
    
    try:
       r['Y_STEP']=str(float(r['Y_STEP'])*alks)
       r['X_STEP']=str(float(r['X_STEP'])*rlks)
    except:
       Geo=0

    f = open(outName+'.rsc','w')
    for k in r.keys():
       f.write(k+'    '+r[k]+'\n')
    f.close()

  elif ext == ('.dem'):
    d,r = readfile.read_dem(file)
    dlks=multilook(d,alks,rlks)

    print 'writing '+outName
    writefile.write_dem(dlks,outName)
    
    r['FILE_LENGTH']=str(dlks.shape[0])
    r['WIDTH']=str(dlks.shape[1])
    r['XMAX']=str(int(r['WIDTH']) - 1)
    r['YMAX']=str(int(r['FILE_LENGTH']) - 1)

    try:
      r['Y_STEP']=str(float(r['Y_STEP'])*alks)
      r['X_STEP']=str(float(r['X_STEP'])*rlks)
    except:
      Geo=0

    f = open(outName+'.rsc','w')
    for k in r.keys():
       f.write(k+'    '+r[k]+'\n')
    f.close()

  elif ext in ['.jpeg','jpg','png']:

    import Image
    im = Image.open(file)

    width = im.size[0] / int(rlks)
    height = im.size[1] / int(alks)

    imlks = im.resize((width, height), Image.NEAREST)
    print 'writing ' + outName
    imlks.save(outName)

    try:
      r=readfile.read_rsc_file(file+'.rsc')
    except:
      sys.exit(1)

    r['FILE_LENGTH']=str(height)
    r['WIDTH']=str(width)
    r['XMAX']=str(int(r['WIDTH']) - 1)
    r['YMAX']=str(int(r['FILE_LENGTH']) - 1)
    try:
      r['Y_STEP']=str(float(r['Y_STEP'])*alks)
      r['X_STEP']=str(float(r['X_STEP'])*rlks)
    except:
      Geo=0
    
    f = open(outName+'.rsc','w')
    for k in r.keys():
       f.write(k+'    '+r[k]+'\n')
    f.close()


  elif ext == ('.h5'):

    h5file=h5py.File(file,'r')
   # outName=file.split('.')[0]+'_a'+str(int(alks))+'lks_r'+str(int(rlks))+'lks.h5'
    h5file_lks=h5py.File(outName,'w')
  
    if 'interferograms' in h5file.keys():
      print 'Filtering the interferograms in space'
      gg = h5file_lks.create_group('interferograms')
      igramList=h5file['interferograms'].keys()
      for igram in igramList:
        print igram
        unwSet = h5file['interferograms'][igram].get(igram)
        unw=unwSet[0:unwSet.shape[0],0:unwSet.shape[1]]
        unw=filter(unw,filtType,par)
        group = gg.create_group(igram)
        dset = group.create_dataset(igram, data=unw, compression='gzip')
        for key, value in h5file['interferograms'][igram].attrs.iteritems():
            group.attrs[key] = value

      dset1=h5file['mask'].get('mask')
      mask=dset1[0:dset1.shape[0],0:dset1.shape[1]]
      group=h5file_lks.create_group('mask')
      dset = group.create_dataset('mask', data=mask, compression='gzip')

    elif 'timeseries' in h5file.keys():
      print 'Filtering the time-series'
      group = h5file_lks.create_group('timeseries')
      dateList=h5file['timeseries'].keys()
      for d in dateList:
        print d
        dset1 = h5file['timeseries'].get(d)
        data=dset1[0:dset1.shape[0],0:dset1.shape[1]]
        data=filter(data,filtType,par)
        
        dset = group.create_dataset(d, data=data, compression='gzip')      

      for key,value in h5file['timeseries'].attrs.iteritems():
        group.attrs[key] = value


      try:
        dset1 = h5file['mask'].get('mask')
        Mask = dset1[0:dset1.shape[0],0:dset1.shape[1]]
       # Masklks=multilook(Mask,alks,rlks)
        group=h5file_lks.create_group('mask')
        dset = group.create_dataset('mask', data=Mask, compression='gzip')
      except:
        print 'Filterd file does not include the maske'


    elif 'temporal_coherence' in h5file.keys() or 'velocity' in h5file.keys() or 'mask' in h5file.keys():
      k=h5file.keys()
      print 'filtering the '+ k[0]
   
      group=h5file_lks.create_group(k[0])    
      dset1 = h5file[k[0]].get(k[0])
      data = dset1[0:dset1.shape[0],0:dset1.shape[1]]
      data=filter(data,filtType,par)
      dset = group.create_dataset(k[0], data=data, compression='gzip')
      for key , value in h5file[k[0]].attrs.iteritems():
         group.attrs[key]=value

    h5file.close()
    h5file_lks.close()

if __name__ == '__main__':

  main(sys.argv[1:])



