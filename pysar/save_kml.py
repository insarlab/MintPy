#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# add 'timeseries'/'wrapped' option, Yunjun, Jul 2015
#

import os
import sys
import h5py
import numpy as np
import matplotlib as mpl              # FA 7/2015: allows plot generation without running an X server
mpl.use('Agg')  
import matplotlib.pyplot as plt
#import matplotlib.pyplot as plt
#import matplotlib.mpl as mpl
#from numpy import round,pi
try:
    from pykml.factory import KML_ElementMaker as KML
except:
    print 'pykml should be installed!'
    sys.exit(1)

from lxml import etree
import getopt
#from numpy import nanmin,nanmax


def rewrap(unw):
   rewrapped = unw - np.round(unw/(2*np.pi)) * 2*np.pi
   return rewrapped

def Usage():
    print '''
***************************************************************
***************************************************************    
  generating  kml kmz files. (needs geocoded files )

  Usage: save_kml.py -f file -m min -M max -d epoch_date -c color_map -i no(yes)
  
  file: a geocoded PySAR product
  -m minmum value
  -M maximum value
  -d date of interferogram or time-series epoch to be converted to kml
     For interferogram, like 971220-990703; for timeseries, like 060924 or 20060924
  -c colormap, jet as default
  -i inverse the colormap
  -w re-wrapping the interferogram [default : yes]
  -r dpi (dots per inch) [default = 500]
  Example:
         
         save_kml.py -f geo_velocity.h5 -m -0.015 -M 0.015 -i yes -c jet -r 250
         save_kml.py -f LoadedData_ChamanT256EnvA6.h5 -d 971220-990703 
         save_kml.py -f Wrapped_ChamanT256EnvA6.h5 -d 971220-990703 
         save_kml.py -f timeseries.h5 -d 20060924         

***************************************************************
***************************************************************
'''

def main(argv):

  color_map='jet'
  disp_opposite = 'no'
  disp_colorbar='yes'
  rewrapping='yes'
  dpi=500

  try:
      opts, args = getopt.getopt(argv,"f:m:M:d:c:w:i:r:")

  except getopt.GetoptError:
      Usage() ; sys.exit(1)
 
  for opt,arg in opts:

      if opt == '-f':
        File = arg
      elif opt == '-m':
        Vmin = float(arg)
      elif opt == '-M':
        Vmax = float(arg)
      elif opt == '-d':
        epoch_date=arg
      elif opt == '-c':
        color_map=arg
      elif opt == '-i':
        disp_opposite=arg
      elif opt == '-w':
        rewrapping=arg
      elif opt == '-r':
        dpi=int(arg)

  try:
    h5file=h5py.File(File,'r')
    k=h5file.keys()
    outName=File.split('.')[0]
  except:
    Usage() ; sys.exit(1)

  print 'Input file is '+k[0]
  ccmap=plt.get_cmap(color_map)


  if k[0] in ('interferograms','wrapped'):

    ifgramList=h5file[k[0]].keys()
    for i in range(len(ifgramList)):
       if epoch_date in ifgramList[i]:
          epoch_number = i
    print ifgramList[epoch_number]
    outName=epoch_date
           
    dset = h5file[k[0]][ifgramList[epoch_number]].get(ifgramList[epoch_number])
    data = dset[0:dset.shape[0],0:dset.shape[1]]

    if k[0] == 'wrapped':
       rewrapping = 'no'
       Vmin = -np.pi    
       Vmax = np.pi

    if rewrapping=='yes':
       data=rewrap(data)
       Vmin = -np.pi	#[-pi,pi] for wrapped interferograms
       Vmax = np.pi
    else:
       try:
          Vmin
       except:
          Vmin = np.nanmin(data)
       try:
          Vmax
       except:
          Vmax = np.nanmax(data)

    try:
       West	=float(h5file[k[0]][ifgramList[epoch_number]].attrs['X_FIRST'])
       North	=float(h5file[k[0]][ifgramList[epoch_number]].attrs['Y_FIRST'])
       lon_step	=float(h5file[k[0]][ifgramList[epoch_number]].attrs['X_STEP'])
       lat_step	=float(h5file[k[0]][ifgramList[epoch_number]].attrs['Y_STEP'])
       lon_unit	=h5file[k[0]][ifgramList[epoch_number]].attrs['Y_UNIT']
       lat_unit	=h5file[k[0]][ifgramList[epoch_number]].attrs['X_UNIT']
       South	=North+lat_step*(data.shape[0]-1)
       East	=West+lon_step*(data.shape[1]-1)
       geocoord='yes'
       print 'Input file is Geocoded.'
    except:
       print '%%%%%%%%%%'
       print 'Error:'
       print 'The input file is not geocoded'
       print ''
       print '%%%%%%%%%%'
       Usage();sys.exit(1)

  elif 'timeseries' in k:
    epochList=h5file['timeseries'].keys()
    for i in range(len(epochList)):
       if epoch_date in epochList[i]:
          epoch_number = i

    ref_date=h5file['timeseries'].attrs['ref_date']
    if len(epoch_date)==8:
       outName=ref_date[2:]+'-'+epoch_date[2:]
    else:
       outName=ref_date[2:]+'-'+epoch_date

    dset = h5file['timeseries'].get(epochList[epoch_number])
    data = dset[0:dset.shape[0],0:dset.shape[1]]
    if rewrapping=='yes':
       data=rewrap(data)

    try:
       Vmin
    except:
       Vmin = np.nanmin(data)
    try:
       Vmax
    except:
       Vmax = np.nanmax(data)

    try:
       West     = float(h5file['timeseries'].attrs['X_FIRST'])
       North    = float(h5file['timeseries'].attrs['Y_FIRST'])
       lon_step = float(h5file['timeseries'].attrs['X_STEP'])
       lat_step = float(h5file['timeseries'].attrs['Y_STEP'])
       lon_unit = h5file['timeseries'].attrs['Y_UNIT']
       lat_unit = h5file['timeseries'].attrs['X_UNIT']
       South = North + lat_step*(data.shape[0]-1)
       East  = West  + lon_step*(data.shape[1]-1)
       South = North + lat_step*(data.shape[0]-1)
       East  = West  + lon_step*(data.shape[1]-1)
       geocoord='yes'
       print 'Input file is Geocoded.'
    except:
       print '%%%%%%%%%%'
       print 'Error:'
       print 'The input file is not geocoded'
       print ''
       print '%%%%%%%%%%'
       Usage();sys.exit(1)


  else:		# one dataset format: velocity, mask, temporal_coherence, rmse, std, etc.

    dset = h5file[k[0]].get(k[0])
    data=dset[0:dset.shape[0],0:dset.shape[1]]
    if disp_opposite in('yes','Yes','Y','y','YES'):
      data=-1*data

    xref=h5file[k[0]].attrs['ref_x']
    yref=h5file[k[0]].attrs['ref_y']

    try:
       Vmin
    except:
       Vmin = np.nanmin(data)
    try:
       Vmax
    except:
       Vmax = np.nanmax(data)

    try:
       West=float(h5file[k[0]].attrs['X_FIRST'])
       North=float(h5file[k[0]].attrs['Y_FIRST'])
       lon_step=float(h5file[k[0]].attrs['X_STEP'])
       lat_step=float(h5file[k[0]].attrs['Y_STEP'])
       lon_unit=h5file[k[0]].attrs['Y_UNIT']
       lat_unit=h5file[k[0]].attrs['X_UNIT']
       South=North+lat_step*(data.shape[0]-1)
       East=West+lon_step*(data.shape[1]-1)     
       geocoord='yes'
       print 'Input file is Geocoded.'
    except:
       print '%%%%%%%%%%'
       print 'Error:'
       print 'The input file is not geocoded'
       print ''
       print '%%%%%%%%%%'
       Usage();sys.exit(1)



#######################################################
  print 'Making png file ...'   
  length = data.shape[0]
  width = data.shape[1]
  fig = plt.figure()
  fig = plt.figure(frameon=False)
 # fig.set_size_inches(width/1000,length/1000)
  ax = plt.Axes(fig, [0., 0., 1., 1.], )
  ax.set_axis_off()
  fig.add_axes(ax)
  
  aspect = width/(length*1.0)
#  ax.imshow(data,aspect='normal')
  
  
  try:
     ax.imshow(data,aspect='normal',vmax=Vmax,vmin=Vmin)
  except:
     ax.imshow(data,aspect='normal')

  ax.set_xlim([0,width])
  ax.set_ylim([length,0])

 # figName = k[0]+'.png'
  figName = outName + '.png'  
  plt.savefig(figName,pad_inches=0.0,dpi=dpi)
 # plt.show()
#############################################################
#Making colorbar
  ######Making the colorbar
  pc = plt.figure(figsize=(1,4))
  axc = pc.add_subplot(111)
  cmap = mpl.cm.jet
  norm = mpl.colors.Normalize(vmin=Vmin*1000, vmax=Vmax*1000)
  clb = mpl.colorbar.ColorbarBase(axc,cmap=cmap,norm=norm, orientation='vertical')
  clb.set_label('mm/yr')
  pc.subplots_adjust(left=0.25,bottom=0.1,right=0.4,top=0.9)
  pc.savefig('colorbar.png',transparent=True,dpi=300)

#############################################################
  print 'generating kml file'
  doc = KML.kml(KML.Folder(KML.name('PySAR product')))
  slc = KML.GroundOverlay(KML.name(figName),KML.Icon(KML.href(figName)),KML.TimeSpan(KML.begin('2003'),KML.end('2010')),KML.LatLonBox(KML.north(str(North)),KML.south(str(South)),KML.east(str(East)),KML.west(str(West))))
  doc.Folder.append(slc)

#############################
  print 'adding colorscale'  
  latdel = North-South
  londel = East-West
  slc1 = KML.GroundOverlay(KML.name('colorbar'),KML.Icon(KML.href('colorbar.png')),KML.LatLonBox(KML.north(str(North-latdel/2.+0.5)), KML.south(str(South+latdel/2.0-0.5)), KML.east(str(West-0.2*londel)), KML.west(str(West-0.4*londel))),KML.altitude('9000'),KML.altitudeMode('absolute'))
  doc.Folder.append(slc1)

#############################

  kmlstr = etree.tostring(doc, pretty_print=True) 
 # kmlname=k[0]+'.kml'
  kmlname = outName + '.kml'
  print 'writing '+kmlname
  kmlfile = open(kmlname,'w')
  kmlfile.write(kmlstr)
  kmlfile.close()

 # kmzName = k[0]+'.kmz'
  kmzName = outName + '.kmz'
  print 'writing '+kmzName
 # cmdKMZ = 'zip ' + kmzName +' '+ kmlname +' ' + figName 
  cmdKMZ = 'zip ' + kmzName +' '+ kmlname +' ' + figName + ' colorbar.png'
  os.system(cmdKMZ)

  cmdClean = 'rm '+kmlname
  os.system(cmdClean)
  cmdClean = 'rm '+figName
  os.system(cmdClean)
  cmdClean = 'rm colorbar.png'
  os.system(cmdClean)

if __name__ == '__main__':

  main(sys.argv[1:])


