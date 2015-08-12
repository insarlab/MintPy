#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import sys
import os
import re
import time
import datetime

import numpy as np
import h5py

import matplotlib.pyplot as plt
import getopt
import _pysar_utilities as ut

def random_selection(M):
  import random
  nrow,ncol=np.shape(M)
  
  y = random.choice(range(nrow))
  x = random.choice(range(ncol))
 
  while M[y,x]==0:
     
     y = random.choice(range(nrow))      
     x = random.choice(range(ncol))

  return y,x

def nearest(x, tbase,xstep):
 # """ find nearest neighbour """
  dist = np.sqrt((tbase -x)**2)
  if min(dist) <= np.abs(xstep):
     indx=dist==min(dist)
  else:
     indx=[]
  return indx

def get_lat_lon(h5file):

   k=h5file.keys()

   if 'interferograms' in k:

     ifgramList = h5file['interferograms'].keys()
     Width=float(h5file['interferograms'][ifgramList[0]].attrs['WIDTH'])
     Length= float(h5file['interferograms'][ifgramList[0]].attrs['FILE_LENGTH'])
     ullon=float(h5file['interferograms'][ifgramList[0]].attrs['X_FIRST'])
     ullat=float(h5file['interferograms'][ifgramList[0]].attrs['Y_FIRST'])
     lon_step=float(h5file['interferograms'][ifgramList[0]].attrs['X_STEP'])
     lat_step=float(h5file['interferograms'][ifgramList[0]].attrs['Y_STEP'])
     lon_unit=h5file['interferograms'][ifgramList[0]].attrs['Y_UNIT']
     lat_unit=h5file['interferograms'][ifgramList[0]].attrs['X_UNIT']

   elif 'timeseries' in k:
 #  Length,Width = np.shape(insarData)
     Width=float(h5file['timeseries'].attrs['WIDTH'])
     Length= float(h5file['timeseries'].attrs['FILE_LENGTH'])
     ullon=float(h5file['timeseries'].attrs['X_FIRST'])
     ullat=float(h5file['timeseries'].attrs['Y_FIRST'])
     lon_step=float(h5file['timeseries'].attrs['X_STEP'])
     lat_step=float(h5file['timeseries'].attrs['Y_STEP'])
     lon_unit=h5file['timeseries'].attrs['Y_UNIT']
     lat_unit=h5file['timeseries'].attrs['X_UNIT']

   elif 'velocity' in k:
 #  Length,Width = np.shape(insarData)
     Width=float(h5file['velocity'].attrs['WIDTH'])
     Length= float(h5file['velocity'].attrs['FILE_LENGTH'])
     ullon=float(h5file['velocity'].attrs['X_FIRST'])
     ullat=float(h5file['velocity'].attrs['Y_FIRST'])
     lon_step=float(h5file['velocity'].attrs['X_STEP'])
     lat_step=float(h5file['velocity'].attrs['Y_STEP'])
     lon_unit=h5file['velocity'].attrs['Y_UNIT']
     lat_unit=h5file['velocity'].attrs['X_UNIT']
  
   lllat=ullat+Length*lat_step
   urlon=ullon+Width*lon_step
   lat=np.arange(ullat,lllat,lat_step)
   lon=np.arange(ullon,urlon,lon_step)
   return lat,lon,lat_step,lon_step


def find_row_column(Lon,Lat,h5file):
   ################################################
  # finding row and column numbers of the GPS point
 
  lat,lon,lat_step,lon_step = get_lat_lon(h5file)
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


def Seeding(h5file,h5file_Seeded,y,x):

  igramList=h5file['interferograms'].keys()
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

  gm = h5file_Seeded.create_group('mask')
  mask = h5file['mask'].get('mask')
  dset = gm.create_dataset('mask', data=mask, compression='gzip')
  
  try:
     
     meanCoherence = h5file['meanCoherence'].get('meanCoherence')
     gc = h5file_Seeded.create_group('meanCoherence')
     dset = gc.create_dataset('meanCoherence', data=meanCoherence, compression='gzip')

  except:
     print 'The Loaded interferograms does not contain the average coherence'

  
def Usage():
  print '''
    ****************************************************************************************
    **************************************************************************************** 
    Referencing all interferograms to the same pixel.
  
    usage:

      seed_data.py -f filename -m method -y lineNumber -x pixelNumber -l latitude -L longitude -s MaskFile
 
    -f : the interferograms, time-series or velocity file saved in hdf5 format.
    -m : the method used for seeding. if x and y are in the input options method is not required.
         method options are: 

         'manual': if method is manual the stack of interferograms
                   is displayed and asks the user to choose a pixel as the reference pixel.

         'auto' or no method specified: A pixel with maximum average coherence is considered as the reference pixel. 
                                        if avergae coherence is not available a random pixel which has  phase value in 
                                        all interferograms is selected.

    -y : line number of the reference pixel (if already is known)
    -x : pixel number of the reference pixel (if already is known)
    -s : Mask file to check if all interferograms have valid phase value in the reference pixel.

    example:
       seed_data.py -f LoadedData_SanAndreasT356EnvD.h5
       seed_data.py -f LoadedData_SanAndreasT356EnvD.h5 -y 257 -x 151
       seed_data.py -f LoadedData_SanAndreasT356EnvD.h5 -m manual
       seed_data.py -f LoadedData_SanAndreasT356EnvD.h5 -l 34.45 -L -116.23
       seed_data.py -f velocity.h5 -y 450 -x 230
       seed_data.py -f velocity.h5 

    ****************************************************************************************
    ****************************************************************************************
    '''

######################################
def main(argv):
  method = 'auto'
  maskFile = 'Mask.h5'
  try:
      opts, args = getopt.getopt(argv,"h:f:m:y:x:l:L:s:")

  except getopt.GetoptError:
      Usage() ; sys.exit(1)

  for opt,arg in opts:
      if opt in ("-h","--help"):
        Usage()
        sys.exit()
      elif opt == '-f':
        file = arg
      elif opt == '-m':
        method = arg
      elif opt == '-y':
        y = int(arg)
      elif opt == '-x':
        x = int(arg)
      elif opt == '-l':
        latr = float(arg)
      elif opt == '-L':
        lonr = float(arg)
      elif opt == '-s':
        maskFile = arg

  try:
     file
  except:
     Usage() ; sys.exit(1)

#  if os.path.isfile('Seeded_'+file):
#      print ''
#      print 'Seeded_'+file+ '  already exists.'
#      print ''
#      sys.exit(1)

################################ 
  h5file = h5py.File(file)
  k=h5file.keys()

  try:
     print 'Finding the row and column number for the lat/lon'
     y,x=find_row_column(lonr,latr,h5file)
     print 'The y and x found for lat lon : ' +str(y) + ' , ' + str(x)
     
  except:
    print 'Skipping lat/lon reference point.'
    print 'Continue with the y/x reference point.' 


  if 'interferograms' in k:
   Mset = h5file['mask'].get('mask')
   M = Mset[0:Mset.shape[0],0:Mset.shape[1]]
   try:
    x
    y
  
  #  h5file = h5py.File(file)
    numIfgrams = len(h5file['interferograms'].keys())
    if numIfgrams == 0.:
       print "There is no data in the file"
       sys.exit(1)
  
   # h5mask=h5py.File(maskFile,'r')
   # Mset=h5mask[h5mask.keys()[0]].get(h5mask.keys()[0])
   # M=Mset[0:Mset.shape[0],0:Mset.shape[1]]
    print 'Checking the reference pixel'
    if M[y][x]==0:
          print '*************************************************************************'    
          print 'ERROR:'
          print 'The slecetd refernce pixel has NaN value in one or more interferograms!'
          print 'Chhose another pixel as the reference pixel.'
          print '*************************************************************************'
          sys.exit(1)

    else:
          print 'Referencing all interferograms to the same pixel at:' + ' y= '+str(y)+' , x= '+str(x)+':'
          h5file_Seeded = h5py.File('Seeded_'+file,'w')
          Seeding(h5file,h5file_Seeded,y,x)
          print 'Done!'
          h5file_Seeded.close()          
          h5file.close()
         # h5mask.close()
          
################################
   
   except:

     # h5mask=h5py.File(maskFile,'r')
     # Mset=h5mask[h5mask.keys()[0]].get(h5mask.keys()[0])
     # M=Mset[0:Mset.shape[0],0:Mset.shape[1]]
      if method=='manual':
         print 'manual selection of the reference point'
         
         h5file = h5py.File(file)
         igramList = h5file['interferograms'].keys()
         stack = ut.stacking(h5file)
         stack[M==0]=np.nan

         fig = plt.figure()
         ax=fig.add_subplot(111)
         ax.imshow(stack)
         
         print 'Click on a pixel that you want to choose as the refernce pixel in the time-series analysis and then close the displayed velocity.'

         SeedingDone='no' 
         def onclick(event):
            if event.button==1:
               print 'click'
               x = int(event.xdata)
               y = int(event.ydata)
               if not np.isnan(stack[y][x]):
                  
                  print 'Referencing all interferograms to the same pixel at:' + ' y= '+str(y)+' , x= '+str(x)+':'
                  h5file_Seeded = h5py.File('Seeded_'+file)   
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

        # Mset=h5file['mask'].get('mask')
        # M=Mset[0:Mset.shape[0],0:Mset.shape[1]]
         
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
         
           C=C*M
           print 'finding a pixel with maximum avergae coherence'
           y,x=np.unravel_index(np.argmax(C), C.shape)

         except:
           y,x=random_selection(M)
             

         print 'Referencing all interferograms to the same pixel at:' + ' y= '+str(y)+' , x= '+str(x)+':'
         h5file_Seeded = h5py.File('Seeded_'+file,'w')
         Seeding(h5file,h5file_Seeded,y,x)
         print 'Done!'
         h5file_Seeded.close()
         h5file.close()
        # plt.imshow(C)
        # plt.plot(x,y,'^',ms=10)
        # plt.show()         
         
  elif 'timeseries' in k:
     print 'Seeding time-series'
     try:
        print 'Seeding time-series epochs to : y=' + str(y) + ' x='+str(x)
     except:
        print 'y and x coordinates of the Seed point are required!'
        sys.exit(1);
 
     h5file_Seeded = h5py.File('Seeded_'+file,'w')
     group=h5file_Seeded.create_group('timeseries')
     dateList=h5file['timeseries'].keys()
     for d in dateList:
        print d
        dset1=h5file['timeseries'].get(d)
        data=dset1[0:dset1.shape[0],0:dset1.shape[1]]
        dset = group.create_dataset(d, data=data-data[y,x], compression='gzip')
     
     for key,value in h5file['timeseries'].attrs.iteritems():
        group.attrs[key] = value
     group.attrs['ref_y']=y
     group.attrs['ref_x']=x

  elif 'velocity' in k:
     Vset=h5file['velocity'].get('velocity')
     V=Vset[0:Vset.shape[0],0:Vset.shape[1]]
     try:
     #   Vset=h5file[h5file.keys()[0]].get(h5file.keys()[0])
     #   V=Vset[0:Vset.shape[0],0:Vset.shape[1]]
        V=V-V[y,x]
        print y
        print x
        outFile= 'seeded_'+file
        h5file2 = h5py.File(outFile,'w')
        group=h5file2.create_group('velocity')
        dset = group.create_dataset('velocity', data=V, compression='gzip')
        for key, value in h5file[k[0]].attrs.iteritems():
            group.attrs[key] = value
        group.attrs['ref_y']=y
        group.attrs['ref_x']=x  
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
               Vset=h5file[h5file.keys()[0]].get(h5file.keys()[0])
               V=Vset[0:Vset.shape[0],0:Vset.shape[1]]
               print V[1][1]
               if not np.isnan(V[y][x]):

                  print 'Referencing all interferograms to the same pixel at:' + ' y= '+str(y)+' , x= '+str(x)+':'

                  h5file2 = h5py.File('seeded_'+file,'w')
                  V=V-V[y][x]
                  group=h5file2.create_group('velocity')
                  dset = group.create_dataset('velocity', data=V, compression='gzip')
                  for key, value in h5file[k[0]].attrs.iteritems():
                      group.attrs[key] = value
                  group.attrs['ref_y']=y
                  group.attrs['ref_x']=x
                  print 'Done!'
                  h5file2.close()
                  SeedingDone='yes'
                  plt.close() # this gic=ves an error message "segementation fault". Although it can be ignored, there should be a better way to close without error message!
               else:
                  print ''
                  print 'warning:'
                  print 'The slecetd refernce pixel has NaN value'
                  print 'Choose another pixel as the reference pixel'


        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show()
       # h5file.close()
     #   h5mask.close()         
######################################
if __name__ == '__main__':

  main(sys.argv[1:])


