#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Dec 2015: Add find_lat_lon(), add 'ref_lat','ref_lon' for geocoded file
# Yunjun, Jan 2016: Add input option for template file
# Yunjun, Apr 2016: Add maskFile input option
# Yunjun, Jun 2016: Add seed_attributes(), support to all file types
#                   Add reference file option


import os
import sys
import getopt

import numpy as np
import h5py

import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._pysar_utilities as ut


########################################## Sub Functions #############################################
###############################################################
def random_selection(M):
  import random
  nrow,ncol=np.shape(M)

  y = random.choice(range(nrow))
  x = random.choice(range(ncol))

  while M[y,x]==0:
     y = random.choice(range(nrow))
     x = random.choice(range(ncol))

  return y,x

###############################################################
def nearest(x, tbase,xstep):
  ## """ find nearest neighbour """
  dist = np.sqrt((tbase -x)**2)
  if min(dist) <= np.abs(xstep):
     indx=dist==min(dist)
  else:
     indx=[]
  return indx

###############################################################
#def get_lat_lon(h5file):
def get_lat_lon(atr):

   Width    = float(atr['WIDTH'])
   Length   = float(atr['FILE_LENGTH'])
   ullon    = float(atr['X_FIRST'])
   ullat    = float(atr['Y_FIRST'])
   lon_step = float(atr['X_STEP'])
   lat_step = float(atr['Y_STEP'])
   lon_unit = atr['Y_UNIT']
   lat_unit = atr['X_UNIT']

   lllat = ullat+Length*lat_step
   urlon = ullon+Width*lon_step
   lat = np.arange(ullat,lllat,lat_step)
   lon = np.arange(ullon,urlon,lon_step)
   return lat,lon,lat_step,lon_step

###############################################################
def find_row_column(Lon,Lat,atr):
  # finding row and column numbers of the GPS point
  lat,lon,lat_step,lon_step = get_lat_lon(atr)
  idx= nearest(Lon, lon, lon_step)
  idy= nearest(Lat, lat, lat_step)
  if idx !=[] and idy != []:
     IDX=np.where(idx==True)[0][0]
     IDY=np.where(idy==True)[0][0]
  else:
     IDX=np.nan
     IDY=np.nan
  return IDY, IDX

################################################
def find_lat_lon(y,x,atr):
  # finding lat and lon of point
  mlat,mlon,lat_step,lon_step = get_lat_lon(atr)
  lat = mlat[y]
  lon = mlon[x]
  return lat, lon

###############################################################
def Seeding(h5file,h5file_Seeded,y,x,lat=[],lon=[]):

  igramList = h5file['interferograms'].keys()
  igramList = sorted(igramList)
  gg = h5file_Seeded.create_group('interferograms')
  for igram in igramList:
      print igram
      unw = h5file['interferograms'][igram].get(igram)
      unw=unw-unw[y][x]
      group = gg.create_group(igram)
      dset = group.create_dataset(igram, data=unw, compression='gzip')
      for key, value in h5file['interferograms'][igram].attrs.iteritems():
          group.attrs[key] = value
      group.attrs['ref_x'] = x
      group.attrs['ref_y'] = y
      if bool(lat) and bool(lon):
         group.attrs['ref_lat'] = lat
         group.attrs['ref_lon'] = lon

  # support old format
  try:
      mask = h5file['mask'].get('mask')
      gm = h5file_Seeded.create_group('mask')
      dset = gm.create_dataset('mask', data=mask, compression='gzip')
      for key, value in mask.attrs.iteritems():
          gm.attrs[key] = value
  except:  print 'No group of mask found in loaded interferograms'

  try:
      meanCoherence = h5file['meanCoherence'].get('meanCoherence')
      gc = h5file_Seeded.create_group('meanCoherence')
      dset = gc.create_dataset('meanCoherence', data=meanCoherence, compression='gzip')
      for key, value in meanCoherence.attrs.iteritems():
          gc.attrs[key] = value
  except:  print 'No group of meanCoherence found in loaded interferograms'

def seed_attributes(atr_in,x,y):
  atr = dict()
  for key, value in atr_in.iteritems():  atr[key] = str(value)
  
  atr['ref_y']=y
  atr['ref_x']=x
  try:
      atr['X_FIRST']
      lat,lon = find_lat_lon(y,x,atr)
      atr['ref_lat']=lat
      atr['ref_lon']=lon
      geocoord='yes'
  except: geocoord='no'

  return atr


#########################################  Usage  ##############################################
def Usage():
  print '''
    ****************************************************************************************
    **************************************************************************************** 
    Referencing all interferograms to the same pixel.
  
    usage:
      seed_data.py -f filename -t templateFile               [ -M MaskFile -m method ]
      seed_data.py -f filename -y lineNumber -x pixelNumber  [ -M MaskFile -m method ]
      seed_data.py -f filename -l latitude   -L longitude    [ -M MaskFile -m method ]
      seed_data.py filename

 
    -f : the interferograms, time-series or velocity file saved in hdf5 format.
    -m : the method used for seeding. if x and y are in the input options method is not required.
         method options are: 

         'manual': if method is manual the stack of interferograms
                   is displayed and asks the user to choose a pixel as the reference pixel.

         'auto' or no method specified: A pixel with maximum average coherence is considered as the reference pixel. 
                                        if avergae coherence is not available a random pixel which has  phase value in 
                                        all interferograms is selected.
         ##### Priority:
             Input mask file > pysar.mask.file > existed Modified_Mask.h5 > existed Mask.h5 > 'mask' group in file

    -y : line number  of the reference pixel (if already is known)
    -x : pixel number of the reference pixel (if already is known)
    -l : latitude     of the reference pixel (for geocoded file)
    -L : longitude    of the reference pixel (for geocoded file)
    -t : template file with setting of reference point information
         Example: pysar.seed.yx   = 1160,300
                  pysar.seed.lalo = 33.1,130.0
         If both pysar.seed.lalo and pysar.seed.yx are setted, pysar.seed.lalo will be used.

        ##### Priority: 
        For radar coord file, only y/x is possible: -y/x > pysar.seed.yx
        For geo   coord file, all option work: -l/L > -y/x > pysar.seed.lalo > pysar.seed.yx
    -r : reference file.


    -M : Mask file to check if all interferograms have valid phase value in the reference pixel.
    -o : output file name [Seeded_file by default]

    example:
       seed_data.py LoadedData.h5
       seed_data.py -f LoadedData.h5 -t ShikokuT417F650_690AlosA.template
       seed_data.py -f LoadedData.h5 -t ShikokuT417F650_690AlosA.template -o Seeded_LoadedData.h5
       seed_data.py -f LoadedData.h5 -y 257   -x 151
       seed_data.py -f LoadedData.h5 -l 34.45 -L -116.23
       seed_data.py -f LoadedData.h5 -M Mask.h5
       seed_data.py -f LoadedData.h5 -m manual
       seed_data.py -f velocity.h5 -y 450 -x 230

    ****************************************************************************************
    ****************************************************************************************
    '''


#######################################  Main Function  ########################################
def main(argv):

  method = 'auto'
  #maskFile = 'Mask.h5'

  ##### Get Input Args
  if len(sys.argv) > 2:
     try:
        opts, args = getopt.getopt(argv,"h:f:m:y:x:l:L:M:t:o:r:")
     except getopt.GetoptError:
        Usage() ; sys.exit(1)

     for opt,arg in opts:
        if   opt in ("-h","--help"):   Usage();  sys.exit()
        elif opt == '-f':        file     = arg
        elif opt == '-m':        method   = arg
        elif opt == '-y':        y        = int(arg)
        elif opt == '-x':        x        = int(arg)
        elif opt == '-l':        latr     = float(arg)
        elif opt == '-L':        lonr     = float(arg)
        elif opt == '-t':        templateFile = arg
        elif opt == '-M':        maskFile = arg
        elif opt == '-o':        outFile  = arg
        elif opt == '-r':        refFile  = arg

  elif len(sys.argv)==2:
     if   argv[0]=='-h':            Usage(); sys.exit(1)
     elif os.path.isfile(argv[0]):  file = argv[0]
     else:  print 'Input file does not existed: '+argv[0];  sys.exit(1)
  elif len(sys.argv)<2:             Usage(); sys.exit(1)

  try:     file
  except:  Usage() ; sys.exit(1)

  try:    outFile
  except: outFile = 'Seeded_'+file

  ##### Read Seed Inputs
  try:
      refFile
      atr_ref = readfile.read_attributes(refFile)
  except: pass

  try:
      templateFile
      templateContents = readfile.read_template(templateFile)
  except: pass

  try:
      latr
      lonr
  except:
      try:
          latr = float(atr_ref['ref_lat'])
          lonr = float(atr_ref['ref_lon'])
      except:
          try: latr,lonr = [float(i) for i in templateContents['pysar.seed.lalo'].split(',')]
          except: pass
  try:
      y
      x
  except:
      try:
          y = int(atr_ref['ref_y'])
          x = int(atr_ref['ref_x'])
      except:
          try: y,x       = [int(i)   for i in templateContents['pysar.seed.yx'].split(',')]
          except: pass

  ##### Judge: use lalo / yx  
  atr = readfile.read_attributes(file)
  h5file = h5py.File(file)
  k=h5file.keys()
  print '************** Reference Point ******************'
  try:
     atr['X_FIRST']      ## Geocoded
     latr
     lonr
     y,x=find_row_column(lonr,latr,atr)
     print 'Reference point: lat = '+str(latr)+', lon = '+str(lonr)
     print 'Its corresponding row and column number: y = '+str(y)+', x = '+str(x)
  except:
     print 'Skipping lat/lon reference point.'
     print 'Continue with the y/x reference point.' 

  ##### Read Mask File 
  ## Priority: Input mask file > pysar.mask.file > existed Modified_Mask.h5 >
  ##                             existed Mask.h5 > 'mask' group in file
  try:       maskFile
  except:
     try:    maskFile = templateContents['pysar.mask.file']
     except:
         if   os.path.isfile('Modified_Mask.h5'):  maskFile = 'Modified_Mask.h5'
         elif os.path.isfile('Mask.h5'):           maskFile = 'Mask.h5'
         else: print 'No mask found!';
  try:  M,Matr = readfile.read(maskFile);  print 'mask: '+maskFile
  except: print 'Can not open mask file: '+maskFile; sys.exit(1)

################################
  import matplotlib.pyplot as plt

  if 'interferograms' in k:
    try:
       x
       y
       numIfgrams = len(h5file['interferograms'].keys())
       if numIfgrams == 0.:     print "There is no data in the file";    sys.exit(1)

       ## Checking the reference pixel
       if M[y][x]==0:
          print '*************************************************************************'    
          print 'ERROR:'
          print 'The slecetd refernce pixel has NaN value in one or more interferograms!'
          print 'Chhose another pixel as the reference pixel.'
          print '*************************************************************************'
          sys.exit(1)
       else:
          print 'Referencing all interferograms to the same pixel at:' +\
                ' y= '+str(y)+' , x= '+str(x)+':'
          h5file_Seeded = h5py.File(outFile,'w');
          if ut.radar_or_geo(file) == 'geo':
              lat,lon = find_lat_lon(y,x,atr)
              Seeding(h5file,h5file_Seeded,y,x,lat,lon)
          else:
              Seeding(h5file,h5file_Seeded,y,x)
          print 'Done!'
          h5file_Seeded.close()          
          h5file.close()
          
    ################################
   
    except:
      if method=='manual':
         print 'manual selection of the reference point'

         igramList = h5file['interferograms'].keys()
         stack = ut.stacking(h5file)
         stack[M==0]=np.nan

         fig = plt.figure();    ax=fig.add_subplot(111);    ax.imshow(stack)
         print 'Click on a pixel that you want to choose as the refernce pixel \
                in the time-series analysis and then close the displayed velocity.'

         SeedingDone='no' 
         def onclick(event):
            if event.button==1:
               print 'click'
               x = int(event.xdata)
               y = int(event.ydata)
               if not np.isnan(stack[y][x]):
                  
                  print 'Referencing all interferograms to the same pixel at:' +\
                        ' y= '+str(y)+' , x= '+str(x)+':'
                  #outFile= 'Seeded_'+file; print 'writing >>> '+outFile
                  h5file_Seeded = h5py.File(outFile,'w')
                  if ut.radar_or_geo(file) == 'geo':
                      lat,lon = find_lat_lon(y,x,atr)
                      Seeding(h5file,h5file_Seeded,y,x,lat,lon)
                  else:
                      Seeding(h5file,h5file_Seeded,y,x)
                  print 'Done!'   
                  h5file_Seeded.close()
                  SeedingDone='yes'
                  plt.close() # this gic=ves an error message "segementation fault". Although it can be ignored, there should be a better way to close without error message!
               else:
                  print ''
                  print 'warning:'
                  print 'The slecetd refernce pixel has NaN value for some interferograms'
                  print 'Choose another pixel as the reference pixel'

         cid = fig.canvas.mpl_connect('button_press_event', onclick)
         plt.show()
         h5file.close()
         h5mask.close()

         if SeedingDone=='no':
            print '''
          **********************************     
          WARNING: interferograms are not referenced to the same pixel yet!
          **********************************
         '''
      else:
         print 'Automatic selection of the reference pixel!'

         #Mset=h5file['mask'].get('mask')
         #M=Mset[0:Mset.shape[0],0:Mset.shape[1]]
         ind0=M==0
         if ind0.sum()==M.shape[0]*M.shape[1]:
            print 'Error:'
            print 'There is no pixel that has valid phase value in all interferograms.' 
            print 'Check the interferograms!'
            print 'Seeding failed'
            sys.exit(1)            

         try:
           Cset=h5file['meanCoherence'].get('meanCoherence')
           C=Cset[0:Cset.shape[0],0:Cset.shape[1]]
           C=C*M;    print 'finding a pixel with maximum avergae coherence'
           y,x=np.unravel_index(np.argmax(C), C.shape)
         except:
           y,x=random_selection(M)

         print 'Referencing all interferograms to the same pixel at:' + ' y= '+str(y)+' , x= '+str(x)+':'
         #outFile= 'Seeded_'+file; print 'writing >>> '+outFile
         h5file_Seeded = h5py.File(outFile,'w')
         if ut.radar_or_geo(file) == 'geo':
             lat,lon = find_lat_lon(y,x,atr)
             Seeding(h5file,h5file_Seeded,y,x,lat,lon)
         else:
             Seeding(h5file,h5file_Seeded,y,x)
         print 'Done!'
         h5file_Seeded.close()
         h5file.close()
         #plt.imshow(C)
         #plt.plot(x,y,'^',ms=10)
         #plt.show()         

  #####################################
  elif 'timeseries' in k:
     print 'Seeding time-series'
     try:      print 'Seeding time-series epochs to : y=' + str(y) + ' x='+str(x)
     except:   print 'y and x coordinates of the Seed point are required!';           sys.exit(1);
     #outFile= 'Seeded_'+file; print 'writing >>> '+outFile 
     h5file_Seeded = h5py.File(outFile,'w')
     group=h5file_Seeded.create_group('timeseries')
     dateList=h5file['timeseries'].keys()
     for d in dateList:
        print d
        dset1=h5file['timeseries'].get(d)
        data=dset1[0:dset1.shape[0],0:dset1.shape[1]]
        dset = group.create_dataset(d, data=data-data[y,x], compression='gzip')
     
     atr = seed_attributes(h5file['timeseries'].attrs,x,y)
     for key,value in atr.iteritems():
        group.attrs[key] = value


  #####################################
  else:
     V, atr = readfile.read(file)
     try:
        V=V-V[y,x]
        atr = seed_attributes(atr,x,y)
        writefile.write(V,atr,outFile)
     except:
        print"Choose the reference point on the screen"
        fig = plt.figure()
        ax=fig.add_subplot(111)
        ax.imshow(V)

        print 'Click on a pixel that you want to choose as the refernce pixel:'
        SeedingDone='no'
        def onclick(event):
            if event.button==1:
               print 'click'
               x = int(event.xdata)
               y = int(event.ydata)
               print V[1][1]
               if not np.isnan(V[y][x]):

                  print 'Referencing all interferograms to the same pixel at:' + ' y= '+str(y)+' , x= '+str(x)+':'
                  V=V-V[y][x]
                  atr = seed_attributes(atr,x,y)
                  writefile.write(V,atr,outFile)
                  SeedingDone='yes'
                  plt.close() # this gic=ves an error message "segementation fault". Although it can be ignored, there should be a better way to close without error message!
               else:
                  print ''
                  print 'warning:'
                  print 'The slecetd refernce pixel has NaN value'
                  print 'Choose another pixel as the reference pixel'


        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show()

  try: h5file.close()
  except: pass
  try: h5file_Seeded.close()
  except: pass
  try: h5mask.close()
  except: pass

  #else: print 'Unrecognized file type: '+k[0]

  print 'writing >>> '+outFile


################################################################################################

if __name__ == '__main__':
  main(sys.argv[1:])


