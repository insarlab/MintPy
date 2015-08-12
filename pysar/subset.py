#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# Add LatLon option by Yunjun, Mar 2015
# Add 'coherence' option, Yunjun, Jul 2015
#

import sys
from  numpy import shape
import getopt
import h5py
import os
import _readfile as readfile
import _writefile as writefile


def Usage():
  print '''
  ***************************************************************************************************
  ***************************************************************************************************

  generating a subset of the dataset

  Usage:
       subset.py -f file -y subset in rows      -x subset in cloumns    -o output name
       subset.py -f file -l subset in latitude  -L subset in longitude  -o output name 

       file: PySAR h5 files [interferograms, coherence, timeseries, velocity, temporal_coherence, rmse, mask] 
             roi_pac files [.unw, .cor, .dem, .hgt]
             image files [jpeg, jpg, png]

  Example:

       subset.py -f Seeded_LoadedData_SanAndreasT356EnvD.h5 -y '400:1500' -x '200:600' -o subste.h5

       subset.py -f Seeded_LoadedData_SanAndreasT356EnvD.h5 -y '400:1500' -x '200:600'

       subset.py -f velocity.h5   -y 1000:1500 -x 350:450
       
       subset.py -f rmse_velocity.h5   -y 1000:1500 -x 350:450

       subset.py -f geo_velocity.h5   -y 1000:1500 -x 350:450
       
       subset.py -f 030405_090801.unw   -y 1000:1500 -x 350:450

       subset.py -f geo_velocity.h5  -l 31.8800:31.9530 -L 130.8510:130.9250

       subset.py -f geo_LoadedData.h5 -l 31.8800:31.9530 -L 130.8510:130.9250

       subset.py -f geo_Coherence.h5 -l 31.8800:31.9530 -L 130.8510:130.9250

  ***************************************************************************************************
  ***************************************************************************************************
  '''


def main(argv):

  #outName='subsetIgrams.h5'
  try:
      opts, args = getopt.getopt(argv,"h:f:x:y:o:l:L:")

  except getopt.GetoptError:
      print 'Error while getting args'
      Usage() ; sys.exit(1)

  for opt,arg in opts:
      if opt in ("-h","--help"):
        Usage()
        sys.exit()
      elif opt == '-f':
        File = arg
      elif opt=='-y':
        ysub=[int(i) for i in arg.split(':')]
        ysub.sort()
      elif opt=='-x':
        xsub = [int(i) for i in arg.split(':')]
        xsub.sort()
      elif opt=='-o':
        outName=arg
      elif opt=='-l':
        Latsub=[float(i) for i in arg.split(':')]
        Latsub.sort()
      elif opt=='-L':
        Lonsub = [float(i) for i in arg.split(':')]
        Lonsub.sort()


#####################################################

  try:
    File
    xsub
    ysub
  except:
    try:
      File
      Latsub
      Lonsub
    except:
      Usage();sys.exit(1)

  try:
    outName
  except:
    outName='subset_'+File

  ext = os.path.splitext(File)[1]

  if ext == '.h5':
    try:
      h5file=h5py.File(File,'r')
    except:
      Usage() ; sys.exit(1)
    k=h5file.keys()

    # convert LatLon to xy for geocoded file
    try:
      Latsub
      Lonsub
      if 'X_FIRST' in h5file[k[0]].attrs.keys():
        xsub=[0]*2
        ysub=[0]*2
        xsub[0]=int((Lonsub[0]-float(h5file[k[0]].attrs['X_FIRST']))/float(h5file[k[0]].attrs['X_STEP']))
        xsub[1]=int((Lonsub[1]-float(h5file[k[0]].attrs['X_FIRST']))/float(h5file[k[0]].attrs['X_STEP']))
        ysub[0]=int((Latsub[1]-float(h5file[k[0]].attrs['Y_FIRST']))/float(h5file[k[0]].attrs['Y_STEP']))
        ysub[1]=int((Latsub[0]-float(h5file[k[0]].attrs['Y_FIRST']))/float(h5file[k[0]].attrs['Y_STEP']))
        print 'Subseting geocoded',ext,' file with Latitude and Longitude...'
      elif 'X_FIRST' in h5file[k[0]][h5file[k[0]].keys()[0]].attrs.keys():	# for geocoded interferograms/coherence
        igramList=h5file[k[0]].keys()
        xsub=[0]*2
        ysub=[0]*2
        xsub[0]=int((Lonsub[0]-float(h5file[k[0]][igramList[0]].attrs['X_FIRST']))/float(h5file[k[0]][igramList[0]].attrs['X_STEP']))
        xsub[1]=int((Lonsub[1]-float(h5file[k[0]][igramList[0]].attrs['X_FIRST']))/float(h5file[k[0]][igramList[0]].attrs['X_STEP']))
        ysub[0]=int((Latsub[1]-float(h5file[k[0]][igramList[0]].attrs['Y_FIRST']))/float(h5file[k[0]][igramList[0]].attrs['Y_STEP']))
        ysub[1]=int((Latsub[0]-float(h5file[k[0]][igramList[0]].attrs['Y_FIRST']))/float(h5file[k[0]][igramList[0]].attrs['Y_STEP']))
        print 'Subseting geocoded',ext,' file with Latitude and Longitude...'
      else:
        print 'Not geocoded file, cannot be subseted with LatLon.'
        Usage() ; sys.exit(1)
    except:
      Geo=0
#    k=h5file.keys()

    if 'interferograms' in k:
  
      igramList=h5file['interferograms'].keys()
      h5out=h5py.File(outName,'w')
      gg=h5out.create_group('interferograms')
      for igram in igramList:
        print igram
        dset1=h5file['interferograms'][igram].get(igram)
        group=gg.create_group(igram)
        dset=group.create_dataset(igram, data=dset1[ysub[0]:ysub[1],xsub[0]:xsub[1]], compression='gzip')
        for key, value in h5file['interferograms'][igram].attrs.iteritems():
          group.attrs[key] = value
        group.attrs['FILE_LENGTH']=shape(dset1[ysub[0]:ysub[1],xsub[0]:xsub[1]])[0]
        group.attrs['WIDTH']=shape(dset1[ysub[0]:ysub[1],xsub[0]:xsub[1]])[1]
        group.attrs['subset_x0']=xsub[0]
        group.attrs['subset_x1']=xsub[1]
        group.attrs['subset_y0']=ysub[0]
        group.attrs['subset_y1']=ysub[1]

        if 'X_FIRST' in h5file['interferograms'][igram].attrs.keys():
          group.attrs['X_FIRST']=float(h5file['interferograms'][igram].attrs['X_FIRST']) + xsub[0]*float(h5file['interferograms'][igram].attrs['X_STEP'])
          group.attrs['Y_FIRST']=float(h5file['interferograms'][igram].attrs['Y_FIRST']) + ysub[0]*float(h5file['interferograms'][igram].attrs['Y_STEP'])  
 
      gm=h5out.create_group('mask')
      try:
        Mset=h5file['mask'].get('mask')
        dset=gm.create_dataset('mask', data=Mset[ysub[0]:ysub[1],xsub[0]:xsub[1]], compression='gzip')
      except:
        print 'No group for mask found! It may cause problem in other processing steps.'

      try:    
        Cset=h5file['meanCoherence'].get('meanCoherence')                  
        gm=h5out.create_group('meanCoherence')
        dset=gm.create_dataset('meanCoherence', data=Cset[ysub[0]:ysub[1],xsub[0]:xsub[1]], compression='gzip')
      except:
        print 'No average coherence found in the File'

    elif k[0] in ('coherence','wrapped'):
      corList=h5file[k[0]].keys()
      h5out=h5py.File(outName,'w')
      gg=h5out.create_group(k[0])
      for cor in corList:
        print cor
        dset1=h5file[k[0]][cor].get(cor)
        group=gg.create_group(cor)
        dset=group.create_dataset(cor, data=dset1[ysub[0]:ysub[1],xsub[0]:xsub[1]], compression='gzip')
        for key, value in h5file[k[0]][cor].attrs.iteritems():
          group.attrs[key] = value
        group.attrs['FILE_LENGTH']=shape(dset1[ysub[0]:ysub[1],xsub[0]:xsub[1]])[0]
        group.attrs['WIDTH']=shape(dset1[ysub[0]:ysub[1],xsub[0]:xsub[1]])[1]
        group.attrs['subset_x0']=xsub[0]
        group.attrs['subset_x1']=xsub[1]
        group.attrs['subset_y0']=ysub[0]
        group.attrs['subset_y1']=ysub[1]

        if 'X_FIRST' in h5file[k[0]][cor].attrs.keys():
          group.attrs['X_FIRST']=float(h5file[k[0]][cor].attrs['X_FIRST']) + xsub[0]*float(h5file[k[0]][cor].attrs['X_STEP'])
          group.attrs['Y_FIRST']=float(h5file[k[0]][cor].attrs['Y_FIRST']) + ysub[0]*float(h5file[k[0]][cor].attrs['Y_STEP'])


    elif 'timeseries' in h5file.keys():
      
      dateList=h5file['timeseries'].keys()
      h5out=h5py.File(outName,'w')
      group=h5out.create_group('timeseries')
      for d in dateList:
        print d
        dset1=h5file['timeseries'].get(d)
        dset=group.create_dataset(d, data=dset1[ysub[0]:ysub[1],xsub[0]:xsub[1]], compression='gzip')
      for key, value in h5file['timeseries'].attrs.iteritems():
        group.attrs[key] = value
      group.attrs['FILE_LENGTH']=shape(dset1[ysub[0]:ysub[1],xsub[0]:xsub[1]])[0]
      group.attrs['WIDTH']=shape(dset1[ysub[0]:ysub[1],xsub[0]:xsub[1]])[1]
      group.attrs['subset_x0']=xsub[0]
      group.attrs['subset_x1']=xsub[1]
      group.attrs['subset_y0']=ysub[0]
      group.attrs['subset_y1']=ysub[1] 

      if 'X_FIRST' in h5file['timeseries'].attrs.keys():
          group.attrs['X_FIRST']=float(h5file['timeseries'].attrs['X_FIRST']) + xsub[0]*float(h5file['timeseries'].attrs['X_STEP'])
          group.attrs['Y_FIRST']=float(h5file['timeseries'].attrs['Y_FIRST']) + ysub[0]*float(h5file['timeseries'].attrs['Y_STEP'])
      h5file.close()
      h5out.close()

    elif 'temporal_coherence' in h5file.keys() or 'velocity' in h5file.keys() or 'mask' in h5file.keys() or 'rmse' in h5file.keys():
      print 'writing  >>>  ' +outName
      dset=h5file[k[0]].get(k[0]) 
      data=dset[ysub[0]:ysub[1],xsub[0]:xsub[1]]
      hfout=h5py.File(outName,'w')
      group= hfout.create_group(k[0])
      group.create_dataset(k[0],data=data,compression='gzip')
    
      for key,value in h5file[k[0]].attrs.iteritems():
         group.attrs[key]=value
      group.attrs['FILE_LENGTH']=data.shape[0]
      group.attrs['WIDTH']=data.shape[1]
      group.attrs['XMIN']=0
      group.attrs['XMAX']=data.shape[1]-1
      group.attrs['YMIN']=0
      group.attrs['YMAX']=data.shape[0]-1
      group.attrs['subset_x0']=xsub[0]
      group.attrs['subset_x1']=xsub[1]
      group.attrs['subset_y0']=ysub[0]
      group.attrs['subset_y1']=ysub[1]
      if 'X_FIRST' in h5file[k[0]].attrs.keys():
         group.attrs['X_FIRST']=float(h5file[k[0]].attrs['X_FIRST']) + xsub[0]*float(h5file[k[0]].attrs['X_STEP'])
         group.attrs['Y_FIRST']=float(h5file[k[0]].attrs['Y_FIRST']) + ysub[0]*float(h5file[k[0]].attrs['Y_STEP'])
      h5file.close()
      hfout.close()
 
  elif ext in ['.unw','.cor','.hgt']:
    a,p,r = readfile.read_float32(File)
  
    try:
      Latsub
      Lonsub
      try:
        r['X_FIRST']
        xsub=[0]*2
        ysub=[0]*2
        xsub[0]=int((Lonsub[0]-float(r['X_FIRST']))/float(r['X_STEP']))
        xsub[1]=int((Lonsub[1]-float(r['X_FIRST']))/float(r['X_STEP']))
        ysub[0]=int((Latsub[1]-float(r['Y_FIRST']))/float(r['Y_STEP']))
        ysub[1]=int((Latsub[0]-float(r['Y_FIRST']))/float(r['Y_STEP']))
        print 'Subseting geocoded',ext,' file with Latitude and Longitude...'
      except:
        print 'Not geocoded file, cannot be subseted with LatLon.'
        Usage() ; sys.exit(1)
    except:
      Geo=0
  
    a=a[ysub[0]:ysub[1],xsub[0]:xsub[1]] 
    p=p[ysub[0]:ysub[1],xsub[0]:xsub[1]]

    print 'writing >>> '+outName
    writefile.write_float32(p,outName)
    
    r['FILE_LENGTH']=str(p.shape[0])
    r['WIDTH']=str(p.shape[1])
    r['XMAX']=str(int(r['WIDTH']) - 1)
    r['YMAX']=str(int(r['FILE_LENGTH']) - 1)
    r['subset_x0']=str(xsub[0]) 
    r['subset_x1']=str(xsub[1])
    r['subset_y0']=str(ysub[0]) 
    r['subset_y1']=str(ysub[1])
    try:
       r['Y_FIRST']=str(float(r['Y_FIRST'])+ysub[0]*float(r['Y_STEP']))
       r['X_FIRST']=str(float(r['X_FIRST'])+xsub[0]*float(r['X_STEP']))
    except:
       Geo=0

    f = open(outName+'.rsc','w')
    for k in r.keys():
       f.write(k+'    '+r[k]+'\n')
    f.close()

  elif ext== '.dem':
    d,r = readfile.read_dem(File)

    try:
      Latsub
      Lonsub
#      print Latsub
      try:
        r['X_FIRST']
        xsub=[0]*2
        ysub=[0]*2
        xsub[0]=int((Lonsub[0]-float(r['X_FIRST']))/float(r['X_STEP']))
        xsub[1]=int((Lonsub[1]-float(r['X_FIRST']))/float(r['X_STEP']))
        ysub[0]=int((Latsub[1]-float(r['Y_FIRST']))/float(r['Y_STEP']))
        ysub[1]=int((Latsub[0]-float(r['Y_FIRST']))/float(r['Y_STEP']))
        print 'Subseting',ext,' file with Latitude and Longitude...'
      except:
        print 'Not geocoded file, cannot be subseted with LatLon.'
        Usage() ; sys.exit(1)
    except:
      Geo=0

    d=d[ysub[0]:ysub[1],xsub[0]:xsub[1]]  

    print 'writing >>> '+outName
    writefile.write_dem(d,outName)
    
    r['FILE_LENGTH']=str(d.shape[0])
    r['WIDTH']=str(d.shape[1])
    r['XMAX']=str(int(r['WIDTH']) - 1)
    r['YMAX']=str(int(r['FILE_LENGTH']) - 1)
    r['subset_x0']=str(xsub[0]) 
    r['subset_x1']=str(xsub[1])
    r['subset_y0']=str(ysub[0]) 
    r['subset_y1']=str(ysub[1])

    try:
       r['Y_FIRST']=str(float(r['Y_FIRST'])+ysub[0]*float(r['Y_STEP']))
       r['X_FIRST']=str(float(r['X_FIRST'])+xsub[0]*float(r['X_STEP']))
    except:
       Geo=0

    f = open(outName+'.rsc','w')
    for k in r.keys():
       f.write(k+'    '+r[k]+'\n')
    f.close()

  elif ext in ['.jpeg','jpg','png']:
    import Image
    im = Image.open(File)

    try:
      r=readfile.read_rsc_file(File+'.rsc')
    except:
      sys.exit(1)

    try:
      Latsub
      Lonsub
      try:
        r['X_FIRST']
        xsub=[0]*2
        ysub=[0]*2
        xsub[0]=int((Lonsub[0]-float(r['X_FIRST']))/float(r['X_STEP']))
        xsub[1]=int((Lonsub[1]-float(r['X_FIRST']))/float(r['X_STEP']))
        ysub[0]=int((Latsub[1]-float(r['Y_FIRST']))/float(r['Y_STEP']))
        ysub[1]=int((Latsub[0]-float(r['Y_FIRST']))/float(r['Y_STEP']))
        print 'Subseting geocoded',ext,' file with Latitude and Longitude...'
      except:
        print 'Not geocoded file, cannot be subseted with LatLon.'
        Usage() ; sys.exit(1)
    except:
      Geo=0

    box = (xsub[0],ysub[0],xsub[1],ysub[1])
    output_img = im.crop(box)
    print 'writing >>> '+outName
    output_img.save(outName)
#    try:
#      r=readfile.read_rsc_file(File+'.rsc')
#    except:
#      sys.exit(1)
     
    r['FILE_LENGTH']=str(ysub[1]-ysub[0])
    r['WIDTH']=str(xsub[1]-xsub[0])
    r['XMAX']=str(int(r['WIDTH']) - 1)
    r['YMAX']=str(int(r['FILE_LENGTH']) - 1)
    r['subset_x0']=str(xsub[0]) 
    r['subset_x1']=str(xsub[1])
    r['subset_y0']=str(ysub[0]) 
    r['subset_y1']=str(ysub[1])
    try:
       r['Y_FIRST']=str(float(r['Y_FIRST'])+ysub[0]*float(r['Y_STEP']))
       r['X_FIRST']=str(float(r['X_FIRST'])+xsub[0]*float(r['X_STEP']))
    except:
       Geo=0

    f = open(outName+'.rsc','w')
    for k in r.keys():
       f.write(k+'    '+r[k]+'\n')
    f.close()


if __name__ == '__main__':

  main(sys.argv[1:])


