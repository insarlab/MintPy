#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# Yunjun, Mar 2015: Add lat/lon option
# Yunjun, Jul 2015: Add 'coherence' option
# Yunjun, Sep 2015: Read .par file for jpg
#                   Add '.mli' and '.slc' option for Gamma product
#                   Make x/y/l/L option independent
#                   Add min/max check for input x/y/lat/lon
#                   Merge all files into PySAR, ROI_PAC, Image and GAMMA
#                   Merge 'interferograms','coherence','wrapped' into one
#                   Merge ROI_APC, Image and GAMMA into one
# Yunjun, Oct 2015: Support '.trans' file

import os
import sys
import getopt

from  numpy import shape
import h5py

import pysar._readfile as readfile
import pysar._writefile as writefile


def Usage():
  print '''
****************************************************************
****************************************************************

  generating a subset of the dataset

  Usage:
       subset.py -f file -y subset_row      -x subset_column    -o output_name
       subset.py -f file -l subset_latitude -L subset_longitude -o output_name 

       file: PySAR h5 files [interferograms, coherence, timeseries, velocity,
                             temporal_coherence, rmse, mask] 
             roi_pac files  [.unw, .cor, .dem, .hgt]
             image files    [jpeg, jpg, png, bmp]
             GAMMA files    [.mli, .slc]

  Example:

       subset.py -f velocity.h5       -y 1000:1500
       subset.py -f velocity.h5       -y 1000:1500       -x 350:450
       subset.py -f LoadedData.h5     -y 400:1500        -x 200:600 -o subste.h5
       subset.py -f 030405_090801.unw -y 1000:1500       -x 350:450
       subset.py -f 101016.mli.jpg    -x 760:1060        -y 620:960
       subset.py -f 101016.slc        -x 760:1060        -y 620:960
       subset.py -f geomap.trans      -l 31.4:31.5       -L 130.5:130.6
       subset.py -f geo_velocity.h5   -l 31.8800:31.9530 -L 130.8510:130.9250


****************************************************************
****************************************************************
  '''


def main(argv):

  if len(sys.argv)>2:
    try:
      opts, args = getopt.getopt(argv,"h:f:x:y:o:l:L:")
    except getopt.GetoptError:
      print 'Error while getting args'
      Usage() ; sys.exit(1)

    for opt,arg in opts:
      if   opt in ("-h","--help"):    Usage() ; sys.exit()
      elif opt == '-f':   File = arg
      elif opt == '-y':   ysub = [int(i) for i in arg.split(':')];      ysub.sort()
      elif opt == '-x':   xsub = [int(i) for i in arg.split(':')];      xsub.sort()
      elif opt == '-o':   outName=arg
      elif opt == '-l':   Latsub = [float(i) for i in arg.split(':')];  Latsub.sort()
      elif opt == '-L':   Lonsub = [float(i) for i in arg.split(':')];  Lonsub.sort()

  else:   Usage(); sys.exit(1)

  try:     outName
  except:  outName='subset_'+File

  ext = os.path.splitext(File)[1]

############################################################################
#################################  PySAR  ##################################

  if ext == '.h5':
    try:      h5file=h5py.File(File,'r')
    except:   Usage() ; sys.exit(1)
    k=h5file.keys()
    if 'interferograms' in k: k[0] = 'interferograms'
    elif 'coherence'    in k: k[0] = 'coherence'
    elif 'timeseries'   in k: k[0] = 'timeseries'
    if k[0] in ('interferograms','coherence','wrapped'):
       atr  = h5file[k[0]][h5file[k[0]].keys()[0]].attrs
    elif k[0] in ('dem','velocity','mask','temporal_coherence','rmse','timeseries'):
       atr  = h5file[k[0]].attrs

    ############# Subset Option #############
    width=int(atr['WIDTH'])
    length=int(atr['FILE_LENGTH'])

    try:
      Latsub
      try:
        lat_step = float(atr['Y_STEP'])
        lat1 = float(atr['Y_FIRST'])
        lat0 = lat1 + length*lat_step
        if Latsub[0]<lat0:  Latsub[0]=lat0; print 'WARNING: input latitude < min ('+str(lat0)+')! Set it to min.'
        if Latsub[1]>lat1:  Latsub[1]=lat1; print 'WARNING: input latitude > max ('+str(lat1)+')! Set it to max.'
        print 'subset in latitude - '+str(Latsub[0])+':'+str(Latsub[1])
        ysub=[0]*2
        ysub[0] = int((Latsub[1]-lat1)/lat_step)
        ysub[1] = int((Latsub[0]-lat1)/lat_step)
      except:
        print 'Not geocoded file, cannot be subseted with LatLon.'
        Usage() ; return
    except:
      try:
        ysub
        if ysub[0]<0:       ysub[0]=0;      print 'WARNING: input y < min (0)! Set it to min.'
        if ysub[1]>length:  ysub[1]=length; print 'WARNING: input y > max ('+str(length)+')! Set it to max.'
        print 'subset in y direction - '+str(ysub[0])+':'+str(ysub[1])
      except: ysub = [0,length]

    try:
      Lonsub
      try:
        lon_step = float(atr['X_STEP'])
        lon0 = float(atr['X_FIRST'])
        lon1 = lon0 + width*lon_step
        if Lonsub[0]<lon0:  Lonsub[0]=lon0; print 'WARNING: input longitude < min ('+str(lon0)+')! Set it to min.'
        if Lonsub[1]>lon1:  Lonsub[1]=lon1; print 'WARNING: input longitude > max ('+str(lon1)+')! Set it to max.'
        print 'subset in longitude - '+str(Lonsub[0])+':'+str(Lonsub[1])
        xsub=[0]*2
        xsub[0]=int((Lonsub[0]-lon0)/lon_step)
        xsub[1]=int((Lonsub[1]-lon0)/lon_step)
      except:
        print 'Not geocoded file, cannot be subseted with LatLon.'
        Usage() ; return
    except:
      try:
        xsub
        if xsub[0]<0:      xsub[0]=0;     print 'WARNING: input x < min (0)! Set it to min.'
        if xsub[1]>width:  xsub[1]=width; print 'WARNING: input x > max ('+str(width)+')! Set it to max x.'
        print 'subset in x direction - '+str(xsub[0])+':'+str(xsub[1])
      except: xsub = [0,width]

    if ysub[0]>length or ysub[1]<0 or xsub[0]>length or xsub[1]<0:
      print 'ERROR: input index is out of data range!'
      print 'range in rdr: x - 0:'+str(width)+'    y - 0:'+str(length)
      try: print 'range in geo: lat - '+str(lat0)+':'+str(lat1)+'    lon - '+str(lon0)+':'+str(lon1)
      except: Geo=0
      sys.exit(1)

    ######## Data Read, Crop and Write #######
    ##### N dset, N attributes
    if k[0] in ('interferograms','coherence','wrapped'):
      print 'writing  >>>  ' +outName
      h5out=h5py.File(outName,'w')
      gg=h5out.create_group(k[0])

      igramList=h5file[k[0]].keys()
      for igram in igramList:
        print igram
        dset1=h5file[k[0]][igram].get(igram)
        group=gg.create_group(igram)
        dset=group.create_dataset(igram, data=dset1[ysub[0]:ysub[1],xsub[0]:xsub[1]], compression='gzip')

        for key, value in h5file[k[0]][igram].attrs.iteritems():    group.attrs[key] = value
        group.attrs['FILE_LENGTH']=ysub[1]-ysub[0]
        group.attrs['WIDTH']      =xsub[1]-xsub[0]
        try:
          sub_x0_ori = int(group.attrs['subset_x0'])
          group.attrs['subset_x0'] = xsub[0]+sub_x0_ori
          group.attrs['subset_x1'] = xsub[1]+sub_x0_ori
        except:
          group.attrs['subset_x0']=xsub[0]
          group.attrs['subset_x1']=xsub[1]
        try:
          sub_y0_ori = int(group.attrs['subset_y0'])
          group.attrs['subset_y0'] = ysub[0]+sub_y0_ori
          group.attrs['subset_y1'] = ysub[1]+sub_y0_ori
        except:
          group.attrs['subset_y0']=ysub[0]
          group.attrs['subset_y1']=ysub[1]
        if 'X_FIRST' in atr.keys():
          group.attrs['X_FIRST']=float(atr['X_FIRST']) + xsub[0]*float(atr['X_STEP'])
          group.attrs['Y_FIRST']=float(atr['Y_FIRST']) + ysub[0]*float(atr['Y_STEP'])  

      ## support of old format
      try:
        Mset=h5file['mask'].get('mask')
        gm=h5out.create_group('mask')
        dset=gm.create_dataset('mask', data=Mset[ysub[0]:ysub[1],xsub[0]:xsub[1]], compression='gzip')
      except:  print 'No group for mask found in the file.'
      try:    
        Cset=h5file['meanCoherence'].get('meanCoherence')                  
        gm=h5out.create_group('meanCoherence')
        dset=gm.create_dataset('meanCoherence', data=Cset[ysub[0]:ysub[1],xsub[0]:xsub[1]], compression='gzip')
      except:  print 'No group for meanCoherence found in the file'

      h5file.close()
      h5out.close()

    ##### N/1 dset, 1 attributes
    elif k[0] in ['timeseries', 'temporal_coherence', 'velocity', 'mask', 'rmse']:
      print 'writing  >>>  ' +outName
      h5out=h5py.File(outName,'w')
      group=h5out.create_group(k[0])

      if k[0] == 'timeseries':      
        dateList=h5file[k[0]].keys()
        for d in dateList:
          print d
          dset1=h5file[k[0]].get(d)
          dset=group.create_dataset(d, data=dset1[ysub[0]:ysub[1],xsub[0]:xsub[1]], compression='gzip')
      elif k[0] in ['temporal_coherence', 'velocity', 'mask', 'rmse']:
        dset1=h5file[k[0]].get(k[0]) 
        dset=group.create_dataset(k[0],data=dset1[ysub[0]:ysub[1],xsub[0]:xsub[1]],compression='gzip')

      ## Update attributes
      for key, value in h5file[k[0]].attrs.iteritems():        group.attrs[key] = value
      group.attrs['FILE_LENGTH']=ysub[1]-ysub[0]
      group.attrs['WIDTH']      =xsub[1]-xsub[0]
      try:
        sub_x0_ori = int(group.attrs['subset_x0'])
        group.attrs['subset_x0'] = xsub[0]+sub_x0_ori
        group.attrs['subset_x1'] = xsub[1]+sub_x0_ori
      except:
        group.attrs['subset_x0']=xsub[0]
        group.attrs['subset_x1']=xsub[1]
      try:
        sub_y0_ori = int(group.attrs['subset_y0'])
        group.attrs['subset_y0'] = ysub[0]+sub_y0_ori
        group.attrs['subset_y1'] = ysub[1]+sub_y0_ori
      except:
        group.attrs['subset_y0']=ysub[0]
        group.attrs['subset_y1']=ysub[1]
      if 'X_FIRST' in atr.keys():
        group.attrs['X_FIRST']=float(atr['X_FIRST']) + xsub[0]*float(atr['X_STEP'])
        group.attrs['Y_FIRST']=float(atr['Y_FIRST']) + ysub[0]*float(atr['Y_STEP'])

      h5file.close()
      h5out.close()


    else:
      print 'Error: group of HDF5 file not recogized!'
      h5file.close()
      Usage() ; sys.exit(1)


############################################################################
#########################  ROI_PAC / Image / GAMMA  ########################

  elif ext in ['.unw','.cor','.hgt','.dem','.trans'] or\
       ext in ['.jpeg','.jpg','.png','.ras','.bmp'] or\
       ext in ['.mli','.slc']:

    try:     atr = readfile.read_rsc_file(File + '.rsc')
    except:
      try:     atr = readfile.read_par_file(File + '.par')
      except:  atr = readfile.read_par_file(os.path.splitext(File)[0] + '.par')

    ############# Subset Option #############
    try:     width  = int(atr['WIDTH']);          length = int(atr['FILE_LENGTH'])
    except:  width  = int(atr['range_samples:']); length = int(atr['azimuth_lines:'])

    try:
      Latsub
      try:
        lat_step = float(atr['Y_STEP'])
        lat1 = float(atr['Y_FIRST'])
        lat0 = lat1 + length*lat_step
        if Latsub[0]<lat0:  Latsub[0]=lat0; print 'WARNING: input latitude < min ('+str(lat0)+')! Set it to min.'
        if Latsub[1]>lat1:  Latsub[1]=lat1; print 'WARNING: input latitude > max ('+str(lat1)+')! Set it to max.'
        print 'subset in latitude - '+str(Latsub[0])+':'+str(Latsub[1])
        ysub=[0]*2
        ysub[0] = int((Latsub[1]-lat1)/lat_step)
        ysub[1] = int((Latsub[0]-lat1)/lat_step)
      except:
        print 'Not geocoded file, cannot be subseted with LatLon.'
        Usage() ; return
    except:
      try:
        ysub
        if ysub[0]<0:       ysub[0]=0;      print 'WARNING: input y < min (0)! Set it to min.'
        if ysub[1]>length:  ysub[1]=length; print 'WARNING: input y > max ('+str(length)+')! Set it to max.'
        print 'subset in y direction - '+str(ysub[0])+':'+str(ysub[1])
      except: ysub = [0,length]

    try:
      Lonsub
      try:
        lon_step = float(atr['X_STEP'])
        lon0 = float(atr['X_FIRST'])
        lon1 = lon0 + width*lon_step
        if Lonsub[0]<lon0:  Lonsub[0]=lon0; print 'WARNING: input longitude < min ('+str(lon0)+')! Set it to min.'
        if Lonsub[1]>lon1:  Lonsub[1]=lon1; print 'WARNING: input longitude > max ('+str(lon1)+')! Set it to max.'
        print 'subset in longitude - '+str(Lonsub[0])+':'+str(Lonsub[1])
        xsub=[0]*2
        xsub[0]=int((Lonsub[0]-lon0)/lon_step)
        xsub[1]=int((Lonsub[1]-lon0)/lon_step)
      except:
        print 'Not geocoded file, cannot be subseted with LatLon.'
        Usage() ; return
    except:
      try:
        xsub
        if xsub[0]<0:      xsub[0]=0;     print 'WARNING: input x < min (0)! Set it to min.'
        if xsub[1]>width:  xsub[1]=width; print 'WARNING: input x > max ('+str(width)+')! Set it to max x.'
        print 'subset in x direction - '+str(xsub[0])+':'+str(xsub[1])
      except: xsub = [0,width]

    if ysub[0]>length or ysub[1]<0 or xsub[0]>length or xsub[1]<0:
      print 'ERROR: input index is out of data range!'
      print 'range in rdr: x - 0:'+str(width)+'    y - 0:'+str(length)
      try: print 'range in geo: lat - '+str(lat0)+':'+str(lat1)+'    lon - '+str(lon0)+':'+str(lon1)
      except: Geo=0
      sys.exit(1)

    ######## Data Read, Crop and Write #######
    print 'writing >>> '+outName
    box = (xsub[0],ysub[0],xsub[1],ysub[1])
    if ext in ['.unw','.cor','.hgt']:
       a,p,r = readfile.read_float32(File,box)
       #p = p[ysub[0]:ysub[1],xsub[0]:xsub[1]]
       writefile.write_float32(p,outName)
    elif ext == '.dem':
       p,r = readfile.read_dem(File)
       p = p[ysub[0]:ysub[1],xsub[0]:xsub[1]]
       writefile.write_dem(p,outName)
    elif ext == '.trans':
       a,p,r = readfile.read_float32(File,box)
       #a = a[ysub[0]:ysub[1],xsub[0]:xsub[1]]
       #p = p[ysub[0]:ysub[1],xsub[0]:xsub[1]]
       writefile.write_float32(a,p,outName)
    elif ext in ['.jpeg','.jpg','.png','.ras','.bmp']: 
       import Image
       im  = Image.open(File)
       box = (xsub[0],ysub[0],xsub[1],ysub[1])
       output_img = im.crop(box)
       output_img.save(outName)
    elif ext == '.mli':
       d,r = readfile.read_gamma_float(File)
       d = d[ysub[0]:ysub[1],xsub[0]:xsub[1]]
       writefile.write_gamma_float(d,outName)
    elif ext == '.slc':
       d,r = readfile.read_gamma_scomplex(File,box)
       writefile.write_gamma_scomplex(d,outName)

    ########### Update .rsc file #############
    atr['FILE_LENGTH'] = str(ysub[1]-ysub[0])
    atr['WIDTH']       = str(xsub[1]-xsub[0])
    atr['XMAX'] = str(width - 1)
    atr['YMAX'] = str(length - 1)
    try:
       sub_x0_ori = int(atr['subset_x0'])
       atr['subset_x0'] = str(xsub[0] + sub_x0_ori)
       atr['subset_x1'] = str(xsub[1] + sub_x0_ori)
    except:
       atr['subset_x0'] = str(xsub[0])
       atr['subset_x1'] = str(xsub[1])
    try:
       sub_y0_ori = int(atr['subset_y0'])
       atr['subset_y0'] = str(ysub[0] + sub_y0_ori)
       atr['subset_y1'] = str(ysub[1] + sub_y0_ori)
    except:
       atr['subset_y0'] = str(ysub[0])
       atr['subset_y1'] = str(ysub[1])
    if 'X_FIRST' in atr.keys():
       atr['Y_FIRST']=str(float(atr['Y_FIRST'])+ysub[0]*float(atr['Y_STEP']))
       atr['X_FIRST']=str(float(atr['X_FIRST'])+xsub[0]*float(atr['X_STEP']))

    f = open(outName+'.rsc','w')
    for k in atr.keys():
       f.write(k+'    '+atr[k]+'\n')
    f.close()

###########################################################################

  else:
    print 'File extension not recogized.'
    Usage() ; sys.exit(1)

###########################################################################

if __name__ == '__main__':
  main(sys.argv[1:])



