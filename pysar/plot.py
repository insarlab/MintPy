#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

#from matplotlib import colors
import getopt
import numpy as np
import h5py
import _readfile as readfile
import _pysar_utilities as ut
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import sys
import os
from matplotlib.colors import LinearSegmentedColormap

def Usage():
   print '''
   plotting the geocoded PySAR product

   plot.py -f velocity.h5  -d dem -m min -M max -x subset -y subst -o outName  -i inverse colormap display (yes or no) -c colomap

   -x : 'xmin:xmax'
   -y : 'ymin:ymax'
   -c : all colormaps in matplotlib is supported (see http://matplotlib.org/examples/pylab_examples/show_colormaps.html)
   
   Examples:

     plot.py -f geo_velocity.h5 -d Sonoran.dem -m -0.01 -M 0.01 -i yes -o plotVelocity.png -c pysar_hsv 

   '''

def main(argv):
     color_map='jet'
     disp_opposite='no'
     try:
        opts, args = getopt.getopt(argv,"h:f:d:o:x:y:m:M:i:c:")

     except getopt.GetoptError:
        Usage() ; sys.exit(1)

     if opts==[]:
        Usage() ; sys.exit(1)
     for opt,arg in opts:
      if opt in ("-h","--help"):
        Usage()
        sys.exit()
      elif opt == '-f':
        File = arg
      elif opt == '-d':
        demFile=arg
      elif opt=='-m':
        Vmin=float(arg)
      elif opt=='-M':
        Vmax=float(arg)
      elif opt == '-x':
        winx=arg.split(':')
      elif opt == '-y':
        winy = arg.split(':')
      elif opt == '-o':
        outName = arg
      elif opt == '-i':
        disp_opposite = arg
      elif opt == '-c':
        color_map=arg 


     h5file=h5py.File(File,'r')
     k=h5file.keys()
     print k[0]
    # ccmap=plt.get_cmap(color_map)

     ################################################
     cdict1 = {'red':   ((0.0, 0.0, 0.0),
                   (0.5, 0.0, 0.0),
                   (0.6, 1.0, 1.0),
                   (0.8, 1.0, 1.0),
                   (1.0, 0.5, 0.5)),

         'green': ((0.0, 0.0, 0.0),
                   (0.2, 0.0, 0.0),
                   (0.4, 1.0, 1.0),
                   (0.6, 1.0, 1.0),
                   (0.8, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.5, .5),
                   (0.2, 1.0, 1.0),
                   (0.4, 1.0, 1.0),
                   (0.5, 0.0, 0.0),
                   (1.0, 0.0, 0.0),)
        }
     if color_map =='pysar_hsv':
         ccmap = LinearSegmentedColormap('BlueRed1', cdict1)
     else:
         ccmap=plt.get_cmap(color_map)
     
     print 'colormap is : '+ color_map
 
     
     ################################################
     dset = h5file[k[0]].get(k[0])
     data=dset[0:dset.shape[0],0:dset.shape[1]]
     if disp_opposite in('yes','Yes','Y','y','YES'):
       data=-1*data
       
     try:
       xref=h5file[k[0]].attrs['ref_x']
       yref=h5file[k[0]].attrs['ref_y']
     except:
       print 'No reference point'

     try:
       ullon=float(h5file[k[0]].attrs['X_FIRST'])
       ullat=float(h5file[k[0]].attrs['Y_FIRST'])
       lon_step=float(h5file[k[0]].attrs['X_STEP'])
       lat_step=float(h5file[k[0]].attrs['Y_STEP'])
       lon_unit=h5file[k[0]].attrs['Y_UNIT']
       lat_unit=h5file[k[0]].attrs['X_UNIT']
       llcrnrlon=ullon
       llcrnrlat=ullat+lat_step*data.shape[0]
       urcrnrlon=ullon+lon_step*data.shape[1]
       urcrnrlat=ullat
       geocoord='yes'
       print 'Input file is Geocoded'
     except:
       geocoord='no'

     fig = plt.figure()
     ax = fig.add_axes([0.1,0.1,0.8,0.8])
     m = Basemap(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                 resolution='l', area_thresh=1., projection='cyl',suppress_ticks=False,ax=ax)

     print demFile
     demFile
     if os.path.basename(demFile).split('.')[1]=='hgt':
           amp,dem,demRsc = readfile.read_float32(demFile)
     elif os.path.basename(demFile).split('.')[1]=='dem':
           dem,demRsc = readfile.read_dem(demFile)

#################################################################

     try:
         winx
         wx=[int(i) for i in win_x.split()]
         dem=dem[:,wx[0]:wx[1]]
         data=data[:,wx[0]:wx[1]]

         ullon=float(h5file[k[0]].attrs['X_FIRST'])+wx[0]
         llcrnrlon=ullon
         urcrnrlon=ullon+lon_step*data.shape[1]

     except:
         print ''

     try:
         winy
         wy=[int(i) for i in winy.split()]
         dem=dem[wy[0]:wy[1],:]
         data=data[wy[0]:wy[1],:]
     except:
         print ''
     
################################################################
     fig = plt.figure()
     ax = fig.add_axes([0.1,0.1,0.8,0.8])
     m = Basemap(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                 resolution='l', area_thresh=1., projection='cyl',suppress_ticks=False,ax=ax)
     cmap_dem=plt.get_cmap('gray')
     m.imshow(ut.hillshade(np.flipud(dem),50.0),cmap=cmap_dem)

      
     try:
       im=m.imshow(np.flipud(data),vmin=Vmin,vmax=Vmax,cmap=ccmap)
      # cb = m.colorbar(im,"right", size="5%", pad='2%')
     except:
       im=m.imshow(np.flipud(data))
      # cb = m.colorbar(im,"right", size="5%", pad='2%')
    # m.bluemarble()
#     cb = m.colorbar(im,"right", size="5%", pad='2%')
    # parallels = np.arange(31.,34,0.5)
    # m.drawparallels(parallels,labels=[1,0,0,1],linewidth=0.0)
    # meridians = np.arange(-115.,-112.,0.5)
    # m.drawmeridians(meridians,labels=[1,0,0,1],linewidth=0.0) 
    # m.drawmapscale()
    # m = Basemap(llcrnrlon=-110.,llcrnrlat=0.,urcrnrlon=-20.,urcrnrlat=57.,
     #       projection='lcc',lat_1=20.,lat_2=40.,lon_0=-60.,
      #      resolution ='l',area_thresh=1000.)

    # m.drawcoastlines()
    # m.drawcountries()
    # m.drawmapboundary(fill_color='#99ffff')
   #  m.fillcontinents(color='#cc9966',lake_color='#99ffff')
    # m.drawparallels(np.arange(10,70,20),labels=[1,1,0,0])
    # m.drawmeridians(np.arange(-100,0,20),labels=[0,0,0,1])
   #  plt.title('Atlantic Hurricane Tracks (Storms Reaching Category 4, 1851-2004)')


     try:
        figName = outName 
     except:
        outName=os.path.basename(File).replace('.h5','')
        figName = outName + '.png'
     plt.savefig(figName,pad_inches=0.0) 
   #  plt.show()

     h5file.close()
if __name__ == '__main__':

  main(sys.argv[1:])



