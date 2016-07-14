#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Aug 2015: Add 'contour','dispFig' options
#                   Finish 'saveFig' option
# Yunjun, Sep 2015: merge all casees into 'plot one' and 'plot all'
#                   Add 'sub_lat/sub_lon' option
#                   Enable '-G' option to display lon/lat
# Yunjun, Oct 2015: Add support of ROI_PAC products, modifiedy from
#                       basic_viewer.py written by Scott
#                   Add extend to colorbar, show value in status bar

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import getopt


############################################################
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

############################################################
def rewrap(unw):
  rewrapped = unw - np.round(unw/(2*np.pi)) * 2*np.pi
  return rewrapped

############################################################
def yymmdd2yyyymmdd(date):
  if date[0] == '9':  date = '19'+date
  else:               date = '20'+date
  return date

############################################################
def Usage():
  print '''
****************************************************************
****************************************************************
  Display PySAR / ROI_PAC products:

  -f: file to display, including:
          PySAR HDF5 files: velocity.h5, timeseries.h5, LoadedData.h5, ...
          ROI_PAC    files: .unw .cor .int .hgt .dem .trans .mli
  -e: display a epoch            of time-series/interferograms (if not specified then all epochs are diplayed).
  -d: display a specific date(s) of time-series/interferograms (if not specified then all epochs are diplayed)
  -m: minimum bound of the colorscale (default is the minimum value of the data set, if set, -w will be no)
  -M: Maximum bound of the colorscale (default is the maximum value of the data set, if set, -w will be no)

  -w: rewrap data to display the time-series epochs and interferograms: yes or no (default is yes)
  -O: Opposite sign. if yes then data is multiplied by -1 and then displayed (default is no)
  -v: flip left-right: yes or no (deafult is no)
  -u: flip up-down:    yes or no (default is no)
  -x: subset in x direction 
  -y: subset in y direction
  -l: subset in latitude
  -L: subset in longitude

  -r: row    number of figures in each window (used only when display the time-series or interferograms. Default: 5)
  -p: column number of figures in each window (used only when display the time-series or interferograms. Default: 8)
  -i: width  space between subplots (default 0.1)
  -j: height space between subplots (default 0.1)
  -s: font size (default: 12 for plot_one, 8 for plot_all)
  -c: colormaps of matplotlib found in (http://matplotlib.org/examples/pylab_examples/show_colormaps.html)
      (default is jet). some options are: seismic, bwr, spectral, jet, ... 
  -t: title of time-series epochs or interferograms. 
      options are 'in' and 'out'. (default is out)

  -G: display geo coordinates for geocoded files: yes or no (default is yes)
  -R: display reference point: yes or no (default is yes)
  -a: color  of marker for the reference point (k,g,r,b,...)
  -b: symbol of marker for the reference point (s,o,^,v,p,x,'*'...)
  -k: size   of marker for the reference point ()

  -P: display the figure (yes or no, default is yes)
  -S: save    the figure (yes or no, default is no)
  -o: output figure name

  -D: dem file
  -C: display Contour when dem file available: no, yes, or only 
      'only' - contour (default)
      'no'   - DEM
      'yes'  - DEM and contour
  -V: contour step (default is 200 meters)
  -g: contour smooth ( Sigma of Gaussian Filter, default is 3.0; Set to 0 for no smoothing) 

  Usage:  
           view.py file.h5
           view.py -f file.h5 
           view.py -f file.h5 -D DEM.dem -C yes
           view.py -f file.h5 -S yes -P no
           view.py -f file.h5 -m minValue -M maxValue -r 5 -p 5
           view.py -f file.h5 -l Latsub -L Lonsub

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  Example:
           view.py velocity.h5
           view.py -f velocity.h5 -m -0.02 -M 0.02 -l yes -c bwr
           view.py -f velocity.h5 -D radar_8rlks.hgt -l yes -m -0.01 -M 0.01 -O yes -y '1 900' -S yes -c hsv -a r -b ^ -k 5
           view.py -f velocity.h5 -D SanAndreas.dem -C only -S yes -P no
           view.py filt_060924-070927.unw
           view.py SanAndreas.dem

           view.py -f timeseries.h5
           view.py -f timeseries.h5 -d 20030502
           view.py -f timeseries.h5 -e 5 
           view.py -f timeseries.h5 -r 5 -p 8 -s 8 -i 0.1 -j 0.1 -w no
           view.py -f timeseries.h5 -d 20030502 -w no -D SanAndreas.dem -C yes -g no
           view.py -f timeseries.h5 -S yes -P no (save figure and not display it)

           view.py -f LoadedData.h5
           view.py -f LoadedData.h5 -e 5 -w no
           view.py -f LoadedData.h5 -d 070927-100217 -w no
           view.py -f LoadedData.h5 -t in
           view.py -f LoadedData.h5 -e 5 -w no -D SanAndreas.dem -C only
           view.py -f Coherence.h5  -e 5
           view.py -f Wrapped.h5    -e 5

   Display a subset:
           view.py -f velocity.h5 -x 100:600     -y 200:800
           view.py -f velocity.h5 -l 31.05:31.10 -L 130.05:130.10
           view.py -f timeseries.h5 -d 20100102      -x 100:600 -y 200:800
           view.py -f LoadedData.h5 -d 070927-100217 -x 100:600 -y 200:800

   Not display reference point:
           view.py -f velocity.h5 -R no -O yes

****************************************************************
****************************************************************           
  '''

############################################################
def main(argv):

  #################  default values  ################
  flip_lr='no'
  flip_ud='no'
  disp_geo = 'yes'
  #font_size=8
  color_map='jet'
  figs_rows=5
  figs_cols=8
  rewrapping='yes'
  allData2display='yes'
  Wspace = 0.1
  Hspace = 0.1
  title = 'out'
  showRef = 'yes'
  ref_color='k'
  ref_symbol='s'
  ref_size =10
  dip_opposite = 'no'
  saveFig='no'
  dispFig='yes'
  dispContour='only'
  contour_step=200
  contour_sigma=3.0
  fig_dpi=300

  ###################  get input args  ###############
  if len(sys.argv)>2:
     try:
        opts, args = getopt.getopt(argv,"h:D:O:G:S:f:m:M:v:u:s:c:e:d:r:p:w:i:j:t:R:a:b:k:x:y:C:V:P:o:g:l:L:")
     except getopt.GetoptError:
        Usage() ; sys.exit(1)
     if opts==[]: Usage() ; sys.exit(1)

     for opt,arg in opts:
        if opt in ("-h","--help"):
           Usage() ; sys.exit()
        elif opt == '-f': File = arg
        elif opt == '-D': demFile=arg
        elif opt == '-w': rewrapping = arg
        elif opt == '-m': min = float(arg);         rewrapping='no'
        elif opt == '-M': max = float(arg);         rewrapping='no'
        elif opt == '-v': flip_lr = arg
        elif opt == '-u': flip_ud = arg
        elif opt == '-s': font_size = int(arg)
        elif opt == '-c': color_map = arg
        elif opt == '-e': epoch_number = int(arg);  allData2display='no'
        elif opt == '-d': epoch_date = arg;         allData2display='no'
        elif opt == '-r': figs_rows = int(arg)
        elif opt == '-p': figs_cols = int(arg)
        elif opt == '-i': Wspace = float(arg)
        elif opt == '-j': Hspace = float(arg)
        elif opt == '-t': title = arg
        elif opt == '-R': showRef = arg
        elif opt == '-a': ref_color = arg
        elif opt == '-b': ref_symbol = arg
        elif opt == '-k': ref_size=int(arg)
        elif opt == '-x': win_x = [int(i) for i in arg.split(':')];      win_x.sort()
        elif opt == '-y': win_y = [int(i) for i in arg.split(':')];      win_y.sort()
        elif opt == '-G': disp_geo = arg
        elif opt == '-O': dip_opposite=arg
        elif opt == '-S': saveFig=arg 
        elif opt == '-C': dispContour=arg
        elif opt == '-V': contour_step=float(arg)
        elif opt == '-P': dispFig=arg
        elif opt == '-o': figName=arg
        elif opt == '-g': contour_sigma = float(arg)
        elif opt == '-l': win_lat = [float(i) for i in arg.split(':')];  win_lat.sort()
        elif opt == '-L': win_lon = [float(i) for i in arg.split(':')];  win_lon.sort()

  elif len(sys.argv)==2:
     if argv[0]=='-h':              Usage(); sys.exit(1)
     elif os.path.isfile(argv[0]):  File = argv[0]
     else:
        print 'Input file does not existed: '+argv[0];  sys.exit(1)
  elif len(sys.argv)<2:             Usage(); sys.exit(1)

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

     from matplotlib.colors import LinearSegmentedColormap
     ccmap = LinearSegmentedColormap('BlueRed1', cdict1)
        
     ################################################
  else:  ccmap=plt.get_cmap(color_map)


  ##################################################
  ext = os.path.splitext(File)[1]
  if ext == '.h5':
     import h5py
     h5file=h5py.File(File,'r')
     k=h5file.keys()
     if   'interferograms' in k: k[0] = 'interferograms'
     elif 'coherence'      in k: k[0] = 'coherence'
     elif 'timeseries'     in k: k[0] = 'timeseries'
     print 'Input: '+str(k)
     if k[0] in ('dem','velocity','mask','temporal_coherence','rmse'):
        allData2display = 'no'

  elif ext in ['.unw','.int','.cor','.hgt','.dem','.trans','.mli','.slc']:
     import pysar._readfile as readfile
     allData2display = 'no'
     k = [ext]
  else: 
     print 'File extension not recogized: '+ext
     print 'Support file format:\n\
                PySAR HDF5 files: velocity.h5, timeseries.h5, LoadedData.h5, ...\n\
                ROI_PAC    files: .unw .cor .int .hgt .dem .trans .mli'
     sys.exit(1)


####################################################################
########################## Plot One ################################

  if allData2display == 'no':
    try:    font_size
    except: font_size=12

    ################# File Reading ##################
    ##### PySAR HDF5
    if k[0] in ('dem','velocity','mask','temporal_coherence','rmse'):
       atr  = h5file[k[0]].attrs
       dset = h5file[k[0]].get(k[0])
       data = dset[0:dset.shape[0],0:dset.shape[1]]
       # rewrapping
       if rewrapping in ('yes','Yes','Y','y','YES'):
          print 'Rewrapping disabled for '+k[0]

    elif k[0] == 'timeseries':
       dateList=h5file[k[0]].keys()
       try:
          epoch_number
       except:
          try:
             epoch_date
             if len(epoch_date)==6:  epoch_date=yymmdd2yyyymmdd(epoch_date)
             epoch_number=dateList.index(epoch_date)
          except:
             print 'Unrecognized epoch input!';  sys.exit(1)
       print 'Displaying date: '+dateList[epoch_number]
       atr  = h5file[k[0]].attrs
       dset = h5file[k[0]].get(dateList[epoch_number])
       data = dset[0:dset.shape[0],0:dset.shape[1]]
       # rewrapping
       if rewrapping in ('yes','Yes','y','Y','YES'):
          print 'Rewrapping. Set min/max to -pi/pi. Showing phase'
          range2phase=4*np.pi/float(atr['WAVELENGTH'])         #double-way, 2*2*pi/lamda
          data=range2phase*data
          data=rewrap(data)
          min = -np.pi
          max = np.pi
          figUnit='(radian)'
       elif rewrapping in ('no','No','N','n','NO'):
          print 'No rewrapping. Showing displacement.'
          figUnit='(m)'

    elif k[0] in ('interferograms','coherence','wrapped'):
       ifgramList=h5file[k[0]].keys()
       try:
          epoch_number
       except:
          for i in range(len(ifgramList)):
             if epoch_date in ifgramList[i]:   epoch_number = i
       print 'Displaying: '+ifgramList[epoch_number]
       atr  = h5file[k[0]][ifgramList[epoch_number]].attrs
       dset = h5file[k[0]][ifgramList[epoch_number]].get(ifgramList[epoch_number])
       data = dset[0:dset.shape[0],0:dset.shape[1]]

       # rewrapping
       if k[0] in ('coherence','wrapped') and rewrapping in ('yes','Yes','y','Y','YES'):
          print 'No rewrapping for coherence/wrapped files, set to "no"'
          rewrapping='no'
       if rewrapping in ('yes','Yes','y','Y','YES'):
          print 'Rewrapping. Set min/max to -pi/pi.'
          data = np.angle(np.exp(1j*data))
          min  = -np.pi
          max  = np.pi

    ##### ROI_PAC product
    elif k[0] in ['.slc','.mli']:
       data,p,atr = readfile.read_complex64(File)
       data= np.nanlog10(data)
       figUnit = '(dB)'
    elif k[0] == '.int':
       a,data,atr = readfile.read_complex64(File)
       min = -np.pi
       max =  np.pi
       rewrapping = 'no'
       figUnit = '(radian)'
    elif k[0] in ['.unw', '.cor', '.hgt', '.trans']:
       a,data,atr = readfile.read_float32(File)
       if   k[0] == '.unw':   figUnit = '(radian)'
       elif k[0] == '.hgt':   figUnit = '(m)'
       if k[0] in ['.cor','.hgt','.trans']: rewrapping='no'
       if rewrapping in ('yes','Yes','y','Y','True','true'):
          print 'Rewrapping. Set min/max to -pi/pi.'
          data = rewrap(data)
          min = -np.pi
          max =  np.pi
    elif k[0] == '.dem':
       data,atr = readfile.read_dem(File)
       figUnit = '(m)'


    ############## Data Option ##################
    # Opposite Sign
    if dip_opposite in ('yes','Yes','Y','y','YES'):
       data=-1*data
    # Subset
    try:      # y/latitude direction
      win_lat
      try:
        atr['Y_FIRST']
        win_y=[0]*2
        win_y[0]=int((win_lat[1]-float(atr['Y_FIRST']))/float(atr['Y_STEP']))
        win_y[1]=int((win_lat[0]-float(atr['Y_FIRST']))/float(atr['Y_STEP']))
        if win_y[0]<0: win_y[0]=0; print 'input latitude > max latitude! Set to max'
        print 'subset in latitude  - '+str(win_lat[0])+':'+str(win_lat[1])
      except:  print 'Non-geocoded file, cannot use LatLon option';   Usage(); sys.exit(1)
    except:
      try:
        win_y
        print 'subset in y direction - '+str(win_y[0])+':'+str(win_y[1])
      except: win_y = [0,int(atr['FILE_LENGTH'])]
    try:      # x/longitude direction
      win_lon
      try:
        atr['X_FIRST']
        win_x=[0]*2
        win_x[0]=int((win_lon[0]-float(atr['X_FIRST']))/float(atr['X_STEP']))
        win_x[1]=int((win_lon[1]-float(atr['X_FIRST']))/float(atr['X_STEP']))
        if win_x[0]<0: win_x[0]=0; print 'input longitude > max longitude! Set to max'
        print 'subset in longitude - '+str(win_lon[0])+':'+str(win_lon[1])
      except:  print 'Non-geocoded file, cannot use LatLon option';   Usage(); sys.exit(1)
    except:
      try:
        win_x
        print 'subset in x direction - '+str(win_x[0])+':'+str(win_x[1])
      except: win_x = [0,int(atr['WIDTH'])]

    data = data[win_y[0]:win_y[1],win_x[0]:win_x[1]]

    # Reference Point
    try:
       xref=atr['ref_x']-win_x[0]
       yref=atr['ref_y']-win_y[0]
    except:  print 'No reference point'
    try:
       xref=xref-atr['subset_x0']
       yref=yref-atr['subset_y0']
    except:  print 'No subset'
    # Geo coordinate
    try:
       lon_step = float(atr['X_STEP'])
       lat_step = float(atr['Y_STEP'])
       lon_unit = atr['Y_UNIT']
       lat_unit = atr['X_UNIT']
       ullon     = float(atr['X_FIRST'])+win_x[0]*lon_step
       ullat     = float(atr['Y_FIRST'])+win_y[0]*lat_step
       llcrnrlon = ullon
       llcrnrlat = ullat+lat_step*data.shape[0]
       urcrnrlon = ullon+lon_step*data.shape[1]
       urcrnrlat = ullat
       geocoord='yes'
       print 'Input file is Geocoded'
    except:  geocoord='no'
    # Flip
    if flip_lr in ('yes','Yes','Y','y','YES'):  data=np.fliplr(data);  xref=np.shape(data)[1]-xref-1 
    if flip_ud in ('yes','Yes','Y','y','YES'):  data=np.flipud(data);  yref=np.shape(data)[0]-yref-1

    # Colorbar Extend
    data_min = np.nanmin(data)
    data_max = np.nanmax(data)
    try:    min
    except: min = data_min
    try:    max
    except: max = data_max
    if   min <= data_min and max >= data_max: cb_extend='neither'
    elif min >  data_min and max >= data_max: cb_extend='min'
    elif min <= data_min and max <  data_max: cb_extend='max'
    else:                                     cb_extend='both'

    ############## DEM Option ##################
    try:
       demFile
       print 'Show topography'
       import pysar._readfile as readfile
       if   os.path.basename(demFile).split('.')[1]=='hgt':  amp,dem,demRsc = readfile.read_float32(demFile)
       elif os.path.basename(demFile).split('.')[1]=='dem':      dem,demRsc = readfile.read_dem(demFile)

       # Subset
       dem = dem[win_y[0]:win_y[1],win_x[0]:win_x[1]]
       # Flip
       if flip_lr in ('yes','Yes','Y','y','YES'):  dem=np.fliplr(dem)
       if flip_ud in ('yes','Yes','Y','y','YES'):  dem=np.flipud(dem)

       # DEM data preparation
       if dispContour in ('no','No','n','N','NO','yes','Yes','y','Y','YES'):           #DEM basemap
          print 'plot DEM as basemap'
          cmap_dem=plt.get_cmap('gray')
          import pysar._pysar_utilities as ut
       if dispContour in ('only','Only','o','O','ONLY','yes','Yes','y','Y','YES'):     #contour
          print 'plot contour'
          #if smoothContour in ('yes','Yes','y','Y','YES'):
          import scipy.ndimage as ndimage
          dem=ndimage.gaussian_filter(dem,sigma=contour_sigma,order=0)
          contour_sequence=np.arange(-6000,9000,contour_step)
    except:  print 'No DEM file'

    ############## Data Plot and Output  ################
    # Figure Title
    if   k[0]=='velocity':                    figTitle = 'Velocity (m/yr)'
    elif k[0]=='temporal_coherence':          figTitle = 'Temporal coherence'
    elif k[0]=='dem':                         figTitle = 'DEM error'
    elif k[0]=='rmse':                        figTitle = 'RMSE (m/yr)'
    elif k[0]=='mask':                        figTitle = 'Pixels with no valid value.'
    elif k[0]=='coherence':                   figTitle = ifgramList[epoch_number]
    elif k[0]in('interferograms','wrapped'):  figTitle = ifgramList[epoch_number]+' (radian)'
    elif k[0]=='timeseries':
       try:    master_date=atr['ref_date']
       except: master_date=atr['DATE']
       if len(master_date)==6:     master_date=yymmdd2yyyymmdd(master_date)
       if dip_opposite in ('yes','Yes','Y','y','YES'): figTitle = dateList[epoch_number]+'_'+master_date+' '+figUnit
       else:                                           figTitle = master_date+'_'+dateList[epoch_number]+' '+figUnit
    elif k[0] in ['.unw','.cor','.hgt','.dem','.trans','.mli','.slc']:
       try:    figTitle = File+' '+figUnit
       except: figTitle = File

    # Plot in Geo-coordinate: plot in map
    if geocoord == 'yes' and disp_geo in ('yes','Yes','y','Y','YES'):
       print 'display Lat/Lon'
       fig = plt.figure()
       ax = fig.add_axes([0.1,0.1,0.8,0.8])
       try: plt.title(figTitle,fontsize=font_size)
       except: pass

       # Map - DEM - Data
       from mpl_toolkits.basemap import Basemap
       m = Basemap(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                   resolution='l', area_thresh=1., projection='cyl',suppress_ticks=False,ax=ax)
       try:
          demFile
          if dispContour in ('no','No','n','N','NO','yes','Yes','y','Y','YES'):
             m.imshow(ut.hillshade(np.flipud(dem),50.0),cmap=cmap_dem)
          if dispContour in ('only','Only','o','O','ONLY','yes','Yes','y','Y','YES'):
             import numpy.matlib
             c_x = np.linspace(llcrnrlon,urcrnrlon,num=dem.shape[1],endpoint='FALSE').reshape(1,dem.shape[1])
             c_xx= np.matlib.repmat(c_x,dem.shape[0],1)
             c_y = np.linspace(llcrnrlat,urcrnrlat,num=dem.shape[0],endpoint='FALSE').reshape(dem.shape[0],1)
             c_yy= np.matlib.repmat(c_y,1,dem.shape[1])
             m.contour(c_xx,c_yy,np.flipud(dem),contour_sequence,origin='lower',colors='black',alpha=0.5,latlon='FALSE')
       except:  pass
       try:     im = m.imshow(np.flipud(data),cmap=ccmap,vmin=min,vmax=max)
       except:  im = m.imshow(np.flipud(data),cmap=ccmap)

       # Reference Point
       if showRef in ('yes','Yes','Y','y','YES'):
          try:
             refPoint=ref_color+ref_symbol
             ref_lon = llcrnrlon + xref*lon_step
             ref_lat = urcrnrlat + yref*lat_step
             plt.plot(ref_lon,ref_lat,refPoint,ms=ref_size)
          except:  print 'No reference point'

       # Colorbar
       from mpl_toolkits.axes_grid1 import make_axes_locatable
       divider = make_axes_locatable(ax)
       cax = divider.append_axes("right",size="5%", pad=0.30)
       plt.colorbar(im,cax=cax,extend=cb_extend)
       #plt.colorbar(im,cax=cax)

       # Status bar
       def format_coord(x,y):
         col = int((x-ullon)/lon_step+0.5)
         row = int((y-ullat)/lat_step+0.5)
         if col>=0 and col<=data.shape[1] and row >=0 and row<=data.shape[0]:
            z = data[row,col]
            return 'lat=%.4f,  lon=%.4f,  value=%.4f'%(x,y,z)
         else:
            return 'lat=%.4f,  lon=%.4f'%(x,y)
       ax.format_coord = format_coord


    # Plot in x/y coordinate: row and column
    else:
       fig = plt.figure()
       ax = fig.add_axes([0.1,0.1,0.8,0.8])
       try: plt.title(figTitle,fontsize=font_size)
       except: pass

       # Plot
       try:
          demFile
          if dispContour in ('no','No','n','N','NO','yes','Yes','y','Y','YES'):
             ax.imshow(ut.hillshade(dem,50.0),cmap=cmap_dem)
          if dispContour in ('only','Only','o','O','ONLY','yes','Yes','y','Y','YES'):
             ax.contour(dem,contour_sequence,origin='lower',colors='black',alpha=0.5)
       except:  pass
       try:     im = ax.imshow(data,cmap=ccmap, vmin=min, vmax=max)
       except:  im = ax.imshow(data,cmap=ccmap)
       cbar = plt.colorbar(im,extend=cb_extend)
       #cbar.set_label('m/yr')

       # Reference Point
       if showRef in ('yes','Yes','Y','y','YES'):
          try:
             refPoint=ref_color+ref_symbol
             ax.plot(xref,yref,refPoint,ms=ref_size)
          except:  print 'No reference point'

       plt.xlim(0,np.shape(data)[1])
       plt.ylim(  np.shape(data)[0],0)

       # Status bar
       def format_coord(x,y):
         col = int(x+0.5)
         row = int(y+0.5)
         if col>=0 and col<=data.shape[1] and row >=0 and row<=data.shape[0]:
            z = data[row,col]
            return 'x=%.4f,  y=%.4f,  value=%.4f'%(x,y,z)
         else:
            return 'x=%.4f,  y=%.4f'%(x,y)
       ax.format_coord = format_coord

    # Save Figure
    if saveFig in ('yes','Yes','Y','y','YES'):
       try:  figName
       except:
          if   k[0]=='velocity':            figName='velocity.pdf'
          elif k[0]=='temporal_coherence':  figName='temporal_coherence.pdf'
          elif k[0]=='dem':                 figName='DEM_error.pdf'
          elif k[0]=='rmse':                figName='rmse.pdf'
          elif k[0]=='mask':                figName='mask.pdf'
          elif k[0]=='timeseries':
             figName=os.path.basename(File).split('.')[0]+'_'+dateList[epoch_number]+'.pdf'
          elif k[0] in ('interferograms','coherence','wrapped'):
             figName=ifgramList[epoch_number]+'.pdf'
          elif k[0] in ['.unw','.cor','.hgt','.dem','.trans','.mli','.slc']:
             figName=File+'.pdf'
       plt.savefig(figName,dpi=fig_dpi)
       print 'Saved figure to '+figName

    # Show Figure
    if dispFig in ('yes','Yes','Y','y','YES'):  
       plt.show()
    
####################################################################
########################## Plot All ################################  

  elif allData2display == 'yes':
    try:    font_size
    except: font_size=8

    if k[0] == 'timeseries':
       atr = h5file[k[0]].attrs
       # rewrapping
       if rewrapping in ('yes','Yes','Y','y','YES'):
          print 'Rewrapping. Set min/max to -pi/pi. Showing phase.' 
          range2phase=4*np.pi/float(h5file['timeseries'].attrs['WAVELENGTH'])
          min=-np.pi
          max=np.pi
       else:  print 'No rewrapping. Showing displacement.'

    elif k[0] in ('interferograms','coherence','wrapped'):
       atr = h5file[k[0]][h5file[k[0]].keys()[0]].attrs
       # rewrapping
       if k[0] in ('coherence','wrapped') and rewrapping in ('yes','Yes','y','Y','YES'):
          print 'No rewrapping for coherence/wrapped files, set to "no"'
          rewrapping='no'
       if rewrapping in ('yes','Yes','Y','y','YES'):
          print 'Rewrapping. Set min/max to -pi/pi.'
          min=-np.pi
          max=np.pi
       else:  print 'No rewrapping'

    ### Subset Option ###
    try:      # y/latitude direction
      win_lat
      try:
        atr['Y_FIRST']
        win_y=[0]*2
        win_y[0]=int((win_lat[1]-float(atr['Y_FIRST']))/float(atr['Y_STEP']))
        win_y[1]=int((win_lat[0]-float(atr['Y_FIRST']))/float(atr['Y_STEP']))
        if win_y[0]<0: win_y[0]=0; print 'input latitude > max latitude! Set to max'
        print 'subset in latitude  - '+str(win_lat[0])+':'+str(win_lat[1])
      except:  print 'Non-geocoded file, cannot use LatLon option';   Usage(); sys.exit(1)
    except:
      try:
        win_y
        print 'subset in y direction - '+str(win_y[0])+':'+str(win_y[1])
      except: win_y = [0,int(atr['FILE_LENGTH'])]
    try:      # x/longitude direction
      win_lon
      try:
        atr['X_FIRST']
        win_x=[0]*2
        win_x[0]=int((win_lon[0]-float(atr['X_FIRST']))/float(atr['X_STEP']))
        win_x[1]=int((win_lon[1]-float(atr['X_FIRST']))/float(atr['X_STEP']))
        if win_x[0]<0: win_x[0]=0; print 'input longitude > max longitude! Set to max'
        print 'subset in longitude - '+str(win_lon[0])+':'+str(win_lon[1])
      except:  print 'Non-geocoded file, cannot use LatLon option';   Usage(); sys.exit(1)
    except:
      try:
        win_x
        print 'subset in x direction - '+str(win_x[0])+':'+str(win_x[1])
      except: win_x = [0,int(atr['WIDTH'])]

    # Figure Name
    if saveFig in ('yes','Yes','Y','y','YES'):
       try:
          figName
          figNameBase = os.path.basename(figName).split('.')[0]
          figNameExt  = os.path.basename(figName).split('.')[1]
       except:
          figNameBase = os.path.basename(File).split('.')[0]
          figNameExt  = '.pdf'

    ################## DEM Options ####################
    try:
       demFile
       print 'Show topography'
       import pysar._readfile as readfile
       if   os.path.basename(demFile).split('.')[1]=='hgt':  amp,dem,demRsc = readfile.read_float32(demFile)
       elif os.path.basename(demFile).split('.')[1]=='dem':      dem,demRsc = readfile.read_dem(demFile)

       # Subset
       dem = dem[win_y[0]:win_y[1],win_x[0]:win_x[1]]

       if dispContour in ('no','No','n','N','NO','yes','Yes','y','Y','YES'):           #DEM basemap
          print 'plot DEM as basemap'
          cmap_dem=plt.get_cmap('gray')
          import pysar._pysar_utilities as ut
          hillshade_dem=ut.hillshade(dem,50.0)
       if dispContour in ('only','Only','o','O','ONLY','yes','Yes','y','Y','YES'):     #contour
          print 'plot contour'
          #if smoothContour in ('yes','Yes','y','Y','YES'):
          import scipy.ndimage as ndimage
          dem=ndimage.gaussian_filter(dem,sigma=contour_sigma,order=0)
          contour_sequence=np.arange(-6000,9000,contour_step)
    except:  print 'No DEM file'

    ################## Plot Loop ####################
    ifgramList=h5file[k[0]].keys()
    nfigs   = figs_rows*figs_cols
    lifgram = len(ifgramList)
    print 'number of  epochs/interferograms to display:'+ str(lifgram)
    kk=int(lifgram/nfigs)+1
    ii=0

    # plot (1,end-1) figures
    for j in range(1,kk):
       fig = plt.figure(j)
       ii=(j-1)*nfigs+1
       for i in range(ii,ii+nfigs):
           print 'loading '+ifgramList[i-1]
           ax = fig.add_subplot(figs_rows,figs_cols,i-ii+1) 

           # Data option
           if k[0] == 'timeseries':
              figTitle = ifgramList[i-1]
              dset = h5file[k[0]].get(ifgramList[i-1])
              data = dset[0:dset.shape[0],0:dset.shape[1]]
              if rewrapping in ('yes','Yes','Y','y','YES'):
                 data=range2phase*data
                 data=rewrap(data)
           elif k[0] in ('interferograms','coherence','wrapped'):
              figTitle = str(i)+' : '+h5file[k[0]][ifgramList[i-1]].attrs['DATE12']
              dset = h5file[k[0]][ifgramList[i-1]].get(ifgramList[i-1])
              data = dset[0:dset.shape[0],0:dset.shape[1]]
              if rewrapping in ('yes','Yes','Y','y','YES'):
                 data=np.angle(np.exp(1j*data))

           data = data[win_y[0]:win_y[1],win_x[0]:win_x[1]]

           # Plot
           try:
              demFile
              if dispContour in ('no','No','n','N','NO','yes','Yes','y','Y','YES'):
                 plt.imshow(hillshade_dem,cmap=cmap_dem)
              if dispContour in ('only','Only','o','O','ONLY','yes','Yes','y','Y','YES'):
                 plt.contour(dem,contour_sequence,origin='lower',colors='black',alpha=0.5)
           except:  pass
           try:     ax.imshow(data,cmap=ccmap,vmin=min,vmax=max)
           except:  ax.imshow(data,cmap=ccmap)

           ax.set_yticklabels([])
           ax.set_xticklabels([])
           ax.set_xticks([])
           ax.set_yticks([])
           if   title=='out':  ax.set_title(figTitle,fontsize=font_size)
           elif title=='in':   add_inner_title(ax, figTitle, loc=1)
       fig.subplots_adjust(wspace=Wspace,hspace=Hspace)
       if saveFig in ('yes','Yes','Y','y','YES'):   
           figName=figNameBase+'_'+str(j)+figNameExt
           plt.savefig(figName,dpi=fig_dpi)

    # plot the last figure
    fig = plt.figure(kk)
    ii=(kk-1)*nfigs+1
    for i in range(ii,lifgram+1):
           print 'loading '+ifgramList[i-1]
           ax = fig.add_subplot(figs_rows,figs_cols,i-ii+1)

           # Data option
           if k[0] == 'timeseries':
              figTitle = ifgramList[i-1]
              dset = h5file[k[0]].get(ifgramList[i-1])
              data = dset[0:dset.shape[0],0:dset.shape[1]]
              if rewrapping in ('yes','Yes','Y','y','YES'):
                 data=range2phase*data
                 data=rewrap(data)
           elif k[0] in ('interferograms','coherence','wrapped'):
              figTitle = str(i)+' : '+h5file[k[0]][ifgramList[i-1]].attrs['DATE12']
              dset = h5file[k[0]][ifgramList[i-1]].get(ifgramList[i-1])
              data = dset[0:dset.shape[0],0:dset.shape[1]]
              if rewrapping in ('yes','Yes','Y','y','YES'):
                 data=np.angle(np.exp(1j*data))

           data = data[win_y[0]:win_y[1],win_x[0]:win_x[1]]

           # Plot
           try:
              demFile
              if dispContour in ('no','No','n','N','NO','yes','Yes','y','Y','YES'):
                 plt.imshow(hillshade_dem,cmap=cmap_dem)
              if dispContour in ('only','Only','o','O','ONLY','yes','Yes','y','Y','YES'):
                 plt.contour(dem,contour_sequence,origin='lower',colors='black',alpha=0.5)
           except:  pass

           try:     ax.imshow(data,cmap=ccmap,vmin=min,vmax=max)
           except:  ax.imshow(data,cmap=ccmap)

           ax.xaxis.label.set_fontsize(20)
           ax.set_yticklabels([])
           ax.set_xticklabels([])
           ax.set_xticks([])
           ax.set_yticks([])
           if    title=='out':  ax.set_title(figTitle,fontsize=font_size)
           elif title =='in':   add_inner_title(ax, figTitle, loc=1)
    fig.subplots_adjust(wspace=Wspace,hspace=Hspace)
    if saveFig in ('yes','Yes','Y','y','YES'):
       figName=figNameBase+'_'+str(kk)+figNameExt
       plt.savefig(figName,dpi=fig_dpi)
       print 'Saved figure to '+figNameBase+'_*'+figNameExt

    if dispFig in ('yes','Yes','Y','y','YES'):
       plt.show()
   
####################################################################
####################################################################  

  try: h5file.close()
  except: pass

if __name__ == '__main__':
  main(sys.argv[1:])


