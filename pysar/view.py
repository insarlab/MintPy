#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import getopt
import h5py
import _readfile as readfile
import _pysar_utilities as ut
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap

def add_inner_title(ax, title, loc, size=None, **kwargs):
    from matplotlib.offsetbox import AnchoredText
    from matplotlib.patheffects import withStroke
    if size is None:
        size = dict(size=plt.rcParams['legend.fontsize'])
    at = AnchoredText(title, loc=loc, prop=size,
                      pad=0., borderpad=0.5,
                      frameon=False, **kwargs)
    ax.add_artist(at)
    at.txt._text.set_path_effects([withStroke(foreground="w", linewidth=3)])
    return at

def rewrap(unw):
     
   rewrapped = unw - np.round(unw/(2*np.pi)) * 2*np.pi
   return rewrapped


def Usage():
  print '''
****************************************************************
****************************************************************
  This function displays the PySAR products:
  velocity.h5, temporal_coherence.h5, rmse.h5, timeseries.h5, Loaded_interferograms.h5 ...

  -f: the hdf5 file to display
  -m: minimum bound of the colorscale (default is the minimum value of the data set)
  -M: Maximum bound of the colorscale (default is the maximum value of the data set)
  -l: flip left-right: yes or no (deafult is no)
  -u: flip up-down: yes or no (default is no)
  -s: font size (default is 8)
  -c: colormaps of matplotlib found in (http://matplotlib.org/examples/pylab_examples/show_colormaps.html) (default is jet).
      some options are: seismic, bwr, spectral, jet, ... 
  -e: display epoch e of the time-series (if not specified then all epochs are diplayed). It's the Same for interferograms.
  -d: display a specific date of the time-series epochs (if not specified then all epochs are diplayed)
  -r: row number of figures in each window (used only when display the time-series epochs or interferograms. Default: 5)
  -p: column number of figures in each window (used only when display the time-series epochs or interferograms. Default: 8)
  -i: width space between subplots (default 0.1)
  -j: height hspace between subplots (default 0.1)
  -w: rewrap data to display the time-series epochs and interferograms: yes or no (default is yes)
  -t: title of time-series epochs or interferograms. options are 'in' and 'out'. (default is out)
  -R: display reference point (yes or no) default is yes
  -a: color of marker for the reference point (k,g,r,b,...)
  -b: symbol of marker for the reference point (s,o,^,v,p,x,'*'...)
  -k: size of marker for the reference point ()
  -x: subste in x direction 
  -y: subset in y direction
  
  -S: save the figure (yes or no)
  -D: dem file
  -G: consider the geo coordinates for geocoded files: yes or no (default is yes)
  -O: Opposite sign. if yes then data is multiplied by -1 and then displayed (default is no)
  Usage:  
           view.py -f file.h5 
           view.py -f file.h5 -m minValue -M maxValue
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  Example:
           view.py -f velocity.h5
           view.py -f velocity.h5 -m -0.02 -M 0.02 -l yes -c bwr
           
           view.py -f timeseries.h5
           view.py -f timeseries.h5 -d 20030502
           view.py -f timeseries.h5 -e 5 
           view.py -f Loaded_SanAndreasT356EnvD.h5
           view.py -f Loaded_SanAndreasT356EnvD.h5 -e 5 -w no
           view.py -f timeseries.h5 -r 5 -p 8 -s 8 -i 0.1 -j 0.1 -w no         
           view.py -f LoadedData_SanAndreasT356EnvD.h5 -t in

           view.py -f velocity_demCor_masked.h5 -D radar_8rlks.hgt -l yes -m -0.01 -M 0.01 -O yes -y '1 900' -S yes -c hsv -a r -b ^ -k 5

   Display a subset:
           view.py -f velocity.h5 -x '100 600' -y '200 800'

   Not display reference point:
           view.py -f velocity.h5 -R no -O yes
           

****************************************************************
****************************************************************           
  '''

def main(argv):
  try:
    opts, args = getopt.getopt(argv,"h:D:O:G:S:f:m:M:l:u:s:c:e:d:r:p:w:i:j:t:R:a:b:k:x:y:")
    
  except getopt.GetoptError:
    Usage() ; sys.exit(1)

  flip_lr='no'
  flip_ud='no'
  disp_geo = 'no'
  font_size=8
  color_map='jet'
  figs_rows=5
  figs_cols=8
  rewrapping='yes'
  allData2display='yes'
  Wspace = 0.1
  Hspace = 0.1
  title = 'out'
#  title = 'None'
  showRef = 'yes' 
  ref_color='k'
  ref_symbol='s'
  ref_size =10
  dip_opposite = 'no'
  saveFig='no'

  if opts==[]:
    Usage() ; sys.exit(1)
  for opt,arg in opts:
    if opt in ("-h","--help"):
      Usage()
      sys.exit()
    elif opt == '-f':
      File = arg
    elif opt == '-D':
      demFile=arg
    elif opt == '-m':
      min = float(arg)
    elif opt == '-M':
      max = float(arg)
    elif opt == '-l':
      flip_lr = arg
    elif opt == '-u':
      flip_ud = arg
    elif opt == '-s':
      font_size = int(arg)
    elif opt == '-c':
      color_map = arg
    elif opt == '-e':
      epoch_number = int(arg)
      allData2display='no'
    elif opt == '-d':
      epoch_date = arg
      allData2display='no'
    elif opt == '-r':
      figs_rows = int(arg)
    elif opt == '-p':
      figs_cols = int(arg)
    elif opt == '-w':
      rewrapping = arg
    elif opt == '-i':
      Wspace = float(arg)
    elif opt == '-j':
      Hspace = float(arg)
    elif opt == '-t':
      title = arg
    elif opt == '-R':
      showRef = arg
    elif opt == '-a':
      ref_color = arg
    elif opt == '-b':
      ref_symbol = arg
    elif opt == 'k':
      ref_size=int(arg)
    elif opt == '-x':
      win_x = arg
    elif opt == '-y':
      win_y = arg
    elif opt == '-G':
      disp_geo = arg
    elif opt == '-O':
      dip_opposite=arg
    elif opt=='-S':
      saveFig=arg 
  

  h5file=h5py.File(File,'r')
  k=h5file.keys()
  print k
  if color_map == 'hsv':
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
        
     ccmap = LinearSegmentedColormap('BlueRed1', cdict1)
        
     ################################################
  else:
     ccmap=plt.get_cmap(color_map)
  
####################################################################
####################################################################
 # if k[0]=='velocity' or k[0]=='temporal_coherence' or k[0]=='rmse':
  if len(k)==1 and k[0] in ('dem','velocity','mask','temporal_coherence','rmse'):  
         
     dset = h5file[k[0]].get(k[0])
     data=dset[0:dset.shape[0],0:dset.shape[1]]
     if dip_opposite in('yes','Yes','Y','y','YES'):
       data=-1*data

     try:
       xref=h5file[k[0]].attrs['ref_x']
       yref=h5file[k[0]].attrs['ref_y']
     except:
       print 'No reference point'

# Yunjun, Mar 2015
     try:
       xref=xref-h5file[k[0]].attrs['subset_x0']
       yref=yref-h5file[k[0]].attrs['subset_y0']
     except:
       print 'No subset'

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
     
     
     

     try:
       win_x
       wx=[int(i) for i in win_x.split()]
       data=data[:,wx[0]:wx[1]]
       xref = xref-wx[0]
     except:
       print 'No subste in x direction' 
     try:
       win_y
       wy=[int(i) for i in win_y.split()]
       data=data[wy[0]:wy[1],:]
       yref = yref-wy[0]
     except:
       print 'No subset in y direction'


     try:
       min
     except:
       min=np.nanmin(data)
     
     try:
       max
     except:
       max=np.nanmax(data)

     if flip_lr=='yes':
         data=np.fliplr(data)
         xref=np.shape(data)[1]-xref-1 
     if flip_ud=='yes':
         data=np.flipud(data)
         yref=np.shape(data)[0]-yref-1
     try:
       demFile
      # amp,dem,demRsc = readfile.read_float32(demFile)
       if os.path.basename(demFile).split('.')[1]=='hgt':
           amp,dem,demRsc = readfile.read_float32(demFile)
       elif os.path.basename(demFile).split('.')[1]=='dem':
           dem,demRsc = readfile.read_dem(demFile)

       try:
         win_x
         wx=[int(i) for i in win_x.split()]
         dem=dem[:,wx[0]:wx[1]]
         
       except:
         print ''
       try:
         win_y
         wy=[int(i) for i in win_y.split()]
         dem=dem[wy[0]:wy[1],:]
         
       except:
         print ''

       if flip_lr=='yes':
          dem=np.fliplr(dem)
       if flip_ud=='yes':
          dem=np.flipud(dem)

       cmap_dem=plt.get_cmap('gray')

       if disp_geo in ('yes','Yes','Y','y','YES') and geocoord in ('yes','Yes','Y','y','YES'):
          print 'display geo'
#          from mpl_toolkits.basemap import Basemap
     #     m = Basemap(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,resolution='f', area_thresh=1., projection='cyl')
      #    m.imshow(ut.hillshade(dem,50.0), interpolation='nearest', origin='upper')
        #  m.drawcoastlines(color='w',linewidth=0.8)
        #  m.drawmapboundary() # draw a line around the map region
        #  m.drawrivers()
        #  m.drawparallels(numpy.arange(int(d1.min()), int(d1.max()), 1),linewidth=0.2,labels=[1,0,0,0])
        #  m.drawmeridians(numpy.arange(int(d0.min()), int(d0.max()), 1),linewidth=0.2,labels=[0,0,0,1])
       else:
          print 'Not GEO'
          plt.imshow(ut.hillshade(dem,50.0),cmap=cmap_dem)
     except:
       print 'No DEM file'



     plt.imshow(data,cmap=ccmap, vmin=min, vmax=max)
     plt.colorbar()

     if k[0]=='velocity':
        plt.title('Velocity (m/yr)',fontsize=font_size)
        figName='velocity.pdf'
     elif k[0]=='temporal_coherence':
        plt.title('Temporal coherence',fontsize=font_size)
        figName='temporal_coherence.pdf'
     elif k[0]=='dem':
        plt.title('DEM error',fontsize=font_size)
        figName='DEM_error.pdf'
     elif k[0]=='rmse':
        plt.title('RMSE (m/yr)',fontsize=font_size)
        figName='rmse.pdf'
     elif k[0]=='mask':
        plt.title('Pixels with no valid value.',fontsize=font_size)
        figName='mask.pdf'
     if showRef=='yes':
        try: 
          refPoint=ref_color+ref_symbol
          plt.plot(xref,yref,refPoint,ms=ref_size)
        except:
          print 'No reference point'

     plt.xlim(0,np.shape(data)[1])
     plt.ylim(np.shape(data)[0],0)
     if saveFig=='yes':
        plt.savefig(figName)
     plt.show()
     
    # plt.savefig('fig.pdf')
 
    # fig = plt.figure()
    # ax.imshow(data,vmin=min, vmax=max)
    # ax.xaxis.label.set_fontsize(40)
####################################################################
####################################################################  

  if 'timeseries' in k and allData2display=='yes':
    
   if rewrapping=='yes':
    print 'rewrapping' 
    dateList=h5file['timeseries'].keys()
    nfigs=figs_rows*figs_cols
    ligram = len(dateList)
    range2phase=4*np.pi/float(h5file['timeseries'].attrs['WAVELENGTH'])
#    range2phase=4*np.pi/0.056
    print 'number of timeseries epochs to display:'+ str(ligram)
    kk=int(ligram/nfigs)+1
    ii=0
    for j in range(1,kk):
       fig = plt.figure(j)
       ii=(j-1)*nfigs+1
       for i in range(ii,ii+nfigs):
            ax = fig.add_subplot(figs_rows,figs_cols,i-ii+1) 
            dset=h5file['timeseries'].get(dateList[i-1])
            data = dset[0:dset.shape[0],0:dset.shape[1]]
            data=range2phase*data
          #  data=np.angle(np.exp(1j*data))
            data=rewrap(data)
            ax.imshow(data,cmap=ccmap)
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            if title=='out':
               ax.set_title(dateList[i-1],fontsize=font_size)
            elif title=='in':
               add_inner_title(ax, dateList[i-1], loc=1)
       fig.subplots_adjust(wspace=Wspace,hspace=Hspace)
       figName=k[0]+'_'+str(j)+'.pdf'
       if saveFig in ['yes','Yes','y','YES']:   
            plt.savefig(figName)


    fig = plt.figure(kk)
    ii=(kk-1)*nfigs+1
    for i in range(ii,ligram+1):
            ax = fig.add_subplot(figs_rows,figs_cols,i-ii+1)
            dset=h5file['timeseries'].get(dateList[i-1])
            data = dset[0:dset.shape[0],0:dset.shape[1]]
            data=range2phase*data
          #  data=np.angle(np.exp(1j*data))
            data=rewrap(data)
            ax.imshow(data,cmap=ccmap)
            ax.xaxis.label.set_fontsize(20)
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            if title=='out':
               ax.set_title(dateList[i-1],fontsize=font_size)
            elif title =='in':
               add_inner_title(ax, dateList[i-1], loc=1)
    fig.subplots_adjust(wspace=Wspace,hspace=Hspace)
    figName=k[0]+'_'+str(kk)+'.pdf'
    if saveFig in ['yes','Yes','y','YES']:
            plt.savefig(figName)
    
    plt.show()
   
   else:
    print 'No rewrapping'
    dateList=h5file['timeseries'].keys()
    nfigs=figs_rows*figs_cols
    ligram = len(dateList)
    print 'number of timeseries epochs to display:'+ str(ligram)
    kk=int(ligram/nfigs)+1
    ii=0
    for j in range(1,kk):
       fig = plt.figure(j)
       ii=(j-1)*nfigs+1
       for i in range(ii,ii+nfigs):

            ax = fig.add_subplot(figs_rows,figs_cols,i-ii+1)
            data=h5file['timeseries'].get(dateList[i-1])
            try:
              im=ax.imshow(data,cmap=ccmap,vmin=min,vmax=max)
           # print 'here'
            except:
              im=ax.imshow(data,cmap=ccmap)
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            if title=='out':
               ax.set_title(dateList[i-1],fontsize=font_size)
            elif title=='in':
               add_inner_title(ax, dateList[i-1], loc=1)
       fig.subplots_adjust(wspace=Wspace,hspace=Hspace)
 
    fig = plt.figure(kk)
    ii=(kk-1)*nfigs+1
    for i in range(ii,ligram+1):

            ax = fig.add_subplot(figs_rows,figs_cols,i-ii+1)
            data=h5file['timeseries'].get(dateList[i-1])
            try:
               im=ax.imshow(data,cmap=ccmap,vmin=min,vmax=max)
            except:
               im=ax.imshow(data,cmap=ccmap)
            ax.xaxis.label.set_fontsize(20)
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            if title=='out':
               ax.set_title(dateList[i-1],fontsize=font_size)
            if title=='in':
               add_inner_title(ax, dateList[i-1], loc=1)

    fig.subplots_adjust(wspace=Wspace,hspace=Hspace)
    plt.show()   


####################################################################
####################################################################  
            
  elif 'timeseries' in k and allData2display=='no':
            
    dateList=h5file['timeseries'].keys()
    try: 
      epoch_number
    except:
      epoch_number=dateList.index(epoch_date)
    range2phase=4*np.pi/float(h5file['timeseries'].attrs['WAVELENGTH'])
  #  range2phase=4*np.pi/0.056
    dset=h5file['timeseries'].get(dateList[epoch_number])
    data = dset[0:dset.shape[0],0:dset.shape[1]]
    if rewrapping=='yes':
       data=range2phase*data
      # data=np.angle(np.exp(1j*data))
       data=rewrap(data)
   
    try:
       min
    except: 
       min=np.nanmin(data)
     
    try:
       max
    except:
       max=np.nanmax(data)


    plt.imshow(data,cmap=ccmap,vmin=min,vmax=max)    
    plt.colorbar()
    plt.show()

################################################################
################################################################

  if k[0]in('interferograms','coherence','wrapped') and allData2display=='yes':

   if k[0] in ('coherence','wrapped'):
    rewrapping='no'
#    color_map = 'gray'
#    ccmap=plt.get_cmap(color_map)
   if rewrapping=='yes':

    ifgramList=h5file[k[0]].keys()
    nfigs=figs_rows*figs_cols
    ligram = len(ifgramList)
    print 'number of '+k[0]+' to display:'+ str(ligram)
    kk=int(ligram/nfigs)+1
    ii=0
    for j in range(1,kk):
       fig = plt.figure(j)
       ii=(j-1)*nfigs+1
       for i in range(ii,ii+nfigs):
            ax = fig.add_subplot(figs_rows,figs_cols,i-ii+1)
            dset = h5file[k[0]][ifgramList[i-1]].get(ifgramList[i-1])
            data = dset[0:dset.shape[0],0:dset.shape[1]] 
            data=np.angle(np.exp(1j*data))         
            ax.imshow(data,cmap=ccmap)
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            if title=='out':
               ax.set_title(h5file[k[0]][ifgramList[i-1]].attrs['DATE12'],fontsize=font_size)
            elif title=='in':
               add_inner_title(ax, h5file[k[0]][ifgramList[i-1]].attrs['DATE12'], loc=1) 
       fig.subplots_adjust(wspace=Wspace,hspace=Hspace)

    fig = plt.figure(kk)
    ii=(kk-1)*nfigs+1
    for i in range(ii,ligram+1):
            ax = fig.add_subplot(figs_rows,figs_cols,i-ii+1)
            dset = h5file[k[0]][ifgramList[i-1]].get(ifgramList[i-1])
            data = dset[0:dset.shape[0],0:dset.shape[1]]
            data=np.angle(np.exp(1j*data))
            ax.imshow(data,cmap=ccmap)
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            if title=='out':
               ax.set_title(h5file[k[0]][ifgramList[i-1]].attrs['DATE12'],fontsize=font_size)
            elif title=='in':
               add_inner_title(ax, h5file[k[0]][ifgramList[i-1]].attrs['DATE12'], loc=1)

    fig.subplots_adjust(wspace=Wspace,hspace=Hspace)
    plt.show()    
 
   else:
    ifgramList=h5file[k[0]].keys()
    nfigs=figs_rows*figs_cols
    ligram = len(ifgramList)
    print 'number of '+k[0]+' to display:'+ str(ligram)
    kk=int(ligram/nfigs)+1
    ii=0
    for j in range(1,kk):
       fig = plt.figure(j)
       ii=(j-1)*nfigs+1
       for i in range(ii,ii+nfigs):
            ax = fig.add_subplot(figs_rows,figs_cols,i-ii+1)
            print 'loading '+ifgramList[i-1]
            dset = h5file[k[0]][ifgramList[i-1]].get(ifgramList[i-1])
            data = dset[0:dset.shape[0],0:dset.shape[1]]
            try:
               ax.imshow(data,vmin=min,vmax=max,cmap=ccmap)
            except:
               ax.imshow(data,cmap=ccmap)
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            if title=='out':
               ax.set_title(h5file[k[0]][ifgramList[i-1]].attrs['DATE12'],fontsize=font_size)
            elif title=='in':
               add_inner_title(ax, h5file[k[0]][ifgramList[i-1]].attrs['DATE12'], loc=1)
       fig.subplots_adjust(wspace=Wspace,hspace=Hspace)
    fig = plt.figure(kk)
    ii=(kk-1)*nfigs+1
    for i in range(ii,ligram+1):
            ax = fig.add_subplot(figs_rows,figs_cols,i-ii+1)
            print 'loading '+ifgramList[i-1]
            dset = h5file[k[0]][ifgramList[i-1]].get(ifgramList[i-1]) 
            data = dset[0:dset.shape[0],0:dset.shape[1]]
            #data = h5file[k[0]][ifgramList[i-1]].get(ifgramList[i-1])

            try:
               ax.imshow(data,vmin=min,vmax=max,cmap=ccmap)
            except:
               ax.imshow(data,cmap=ccmap)
            
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            if title=='out':
               ax.set_title(h5file[k[0]][ifgramList[i-1]].attrs['DATE12'],fontsize=font_size)
            elif title=='in':
               add_inner_title(ax, h5file[k[0]][ifgramList[i-1]].attrs['DATE12'], loc=1)

    fig.subplots_adjust(wspace=Wspace,hspace=Hspace)
    plt.show()
  
  ####################################################################
####################################################################  

  elif k[0]in('interferograms','coherence','wrapped') and allData2display=='no':

    if k[0] in ('coherence','wrapped'):
      rewrapping=='no'

    ifgramList=h5file[k[0]].keys()
    try:
      epoch_number
    except:
      for i in range(len(ifgramList)):
        if epoch_date in ifgramList[i]:
           epoch_number = i
    dset = h5file[k[0]][ifgramList[epoch_number]].get(ifgramList[epoch_number])
    data = dset[0:dset.shape[0],0:dset.shape[1]]
    if rewrapping=='yes':
       data=np.angle(np.exp(1j*data))

    if dip_opposite in('yes','Yes','Y','y','YES'):
       data=-1*data

    #DEM basemap
    try:
       demFile
       if os.path.basename(demFile).split('.')[1]=='hgt':
           amp,dem,demRsc = readfile.read_float32(demFile)
       elif os.path.basename(demFile).split('.')[1]=='dem':
           dem,demRsc = readfile.read_dem(demFile)

       try:
         win_x
         wx=[int(i) for i in win_x.split()]
         dem=dem[:,wx[0]:wx[1]]
       except:
         print ''
       try:
         win_y
         wy=[int(i) for i in win_y.split()]
         dem=dem[wy[0]:wy[1],:]
       except:
         print ''

       if flip_lr=='yes':
          dem=np.fliplr(dem)
       if flip_ud=='yes':
          dem=np.flipud(dem)

       cmap_dem=plt.get_cmap('gray')

       if disp_geo in ('yes','Yes','Y','y','YES') and geocoord in ('yes','Yes','Y','y','YES'):
          print 'display geo'
       else:
          print 'Not GEO'
          plt.imshow(ut.hillshade(dem,50.0),cmap=cmap_dem)
    except:
       print 'No DEM file'

    try:
        plt.imshow(data,cmap=ccmap,vmin=min,vmax=max)
    except:
        plt.imshow(data,cmap=ccmap)
    
    plt.colorbar()
#    plt.title(h5file[k[0]][ifgramList[epoch_number]].attrs['DATE12'],fontsize=font_size)
    plt.title(ifgramList[epoch_number],fontsize=font_size)
    plt.show()
    
################################################################
################################################################

  h5file.close()
if __name__ == '__main__':

  main(sys.argv[1:])


