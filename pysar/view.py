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
# Yunjun, Nov 2015: Add data range output
# Yunjun, Dec 2015: Add double-date -d option for timeseries file
#                   Add long options
# Yunjun, Jan 2016: Support multiple epoch display for -d -e option
#                   Change -R option from show reference point to set reference date
#                   Add -t template option 
# Yunjun, May 2016: Add unit_and_scale(), -u option
#                   Change seldom used/bool judge option from - to -- option
#                   Use pysar.subset, readfile.read and box option
# Yunjun, Jun 2016: Use multilook from pysar.multi_looking for multiple display
#                   Add --point/line option to plot points and lines
#                   Add orbit_direction(), default flip based on asc_desc
#                   Simplified code for multiple plots
# Yunjun, Jul 2016: add --mask input option


import sys
import os
import getopt

import h5py
import numpy as np
import matplotlib.pyplot as plt

import pysar._readfile as readfile
import pysar._pysar_utilities as ut
import pysar.subset as subset
import pysar._datetime as ptime


###################################  Sub Function  ######################################
##################################################
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

##################################################
def rewrap(data,atr):
  if atr['UNIT'] == 'm':
      range2phase = -4*np.pi/float(atr['WAVELENGTH'])         #double-way, -2*2*pi/lamda
      data = range2phase*data
      rewrapped = data - np.round(data/(2*np.pi)) * 2*np.pi
  elif atr['UNIT'] == 'radian':
      rewrapped = data - np.round(data/(2*np.pi)) * 2*np.pi
  else: print 'Can not rewrap file in unit: '+atr['UNIT']
  #print 'rewrapping to -pi/pi'
  return rewrapped

##################################################
def unit_and_scale(data_unit,display_unit):
  ## Calculate the scale factor based on data file's unit and display unit
  ## Default data file units in PySAR are:
  ## m, m/yr, radian, 1

  display_scale = 1
  #if not len(data_unit.split('/')) == len(display_unit.split('/')):
  #    print 'Data unit and display unit are not compatible. use original unit.'
  #    return display_scale, data_unit

  display_unit = display_unit.split('/')
  data_unit    =    data_unit.split('/')

  if   display_unit[0] == 'mm': display_scale = display_scale*1000.0
  elif display_unit[0] == 'cm': display_scale = display_scale*100.0
  elif display_unit[0] == 'm' : pass
  elif display_unit[0] == 'km': display_scale = display_scale*0.001
  else: print 'Unrecognized length unit: '+display_unit[0]

  try:
      if   display_unit[1] in ['y','yr','year'  ]: display_unit[1] = 'yr' ; pass
      elif display_unit[1] in ['m','mon','month']: display_unit[1] = 'mon'; display_scale = display_scale/12.0
      elif display_unit[1] in ['d','day'        ]: display_unit[1] = 'day'; display_scale = display_scale/365.25
      else: print 'Unrecognized time unit: '+display_unit[1]
  except: pass

  ######
  try:
      display_unit = display_unit[0]+'/'+display_unit[1]
  except:
      if len(data_unit) == 1:  display_unit = display_unit[0]
      else:                    display_unit = display_unit[0]+'/yr'

  return display_unit, display_scale

##################################################
def unit_type(unit_in):
  unit = unit_in.split('/')
  if   len(unit) == 1:
      if   unit[0] in ['radian']:           type = 'angle'
      elif unit[0] in ['km','m','cm','mm']: type = 'length'
      elif unit[0] in ['1','dB']:           type = '1'
      else: print 'Unrecognized unit type: '+unit_in; type = 'unknown'
  elif len(unit) == 2:
      if unit[0] in ['km','m','cm','mm'] and unit[1] in ['yr','mon','day']:
          type = 'velocity'
      else: print 'Unrecognized unit type: '+unit_in; type = 'unknown'
  else: print 'Unrecognized unit type: '+unit_in; type = 'unknown'

  return type

##################################################
def orbit_direction(atr):
  try:
      atr['ORBIT_DIRECTION']
      if   atr['ORBIT_DIRECTION'].lower() in ['ascending', 'ascend' ]:   return  'ascending'
      elif atr['ORBIT_DIRECTION'].lower() in ['descending','descend']:   return 'descending'
      else: print 'Unrecognized attribute: ORBIT_DIRECTION = '+atr['ORBIT_DIRECTION']; sys.exit(1)
  except: pass

  try: heading = float(atr['HEADING'])
  except: pass
  try: heading = float(atr['HEADING_DEG'])
  except: pass
  try:
      heading
      if abs(heading) < 90: return  'ascending'
      else:                 return 'descending'
  except: print 'Cannot found HEADING or HEADING_DEG attribute.'; sys.exit(1)


##################  Usage  #######################
def Usage():
  print '''
*****************************************************************************************

  Display PySAR / ROI_PAC products:


  -f            : file to display, including:
                  PySAR HDF5 files: velocity.h5, timeseries.h5, LoadedData.h5, ...
                  ROI_PAC    files: .unw .cor .int .hgt .dem .trans .mli
  -d            : display a specific date(s)     of time-series/interferograms (if not specified then all epochs are diplayed)
  -e            : display a epoch (start from 1) of time-series/interferograms (if not specified then all epochs are diplayed).
  -E/--exclude  : exclude epoch list for timeseries/interferograms.
                  Set '-E --' to disable the exclude list in template, i.e.
                      view.py -f timeseries.h5 -t KyushuT73F2980AlosD.template -E --
  -m            : minimum bound of the colorscale (default is the minimum value of the data set, if set, -w will be no)
  -M            : Maximum bound of the colorscale (default is the maximum value of the data set, if set, -w will be no)
  -t            : template file, i.e.
                  pysar.view.row     = 5
                  pysar.view.column  = 20
                  pysar.view.min     = -7
                  pysar.view.max     = 7
  --mask        : mask file for display (useful when displaying un-masked files)

  DEM:
  -D               : dem file (show shaded relief by default)
  --dem-contour    : show DEM contour
  --dem-noshade    : do not show DEM shaded relief
  --contour-step   : contour step                      (default is 200 meters)
  --contour-smooth : contour smooth ( Sigma of Gaussian Filter, default is 3.0; Set to 0 for no smoothing) 

  Data Option:
  --wrap        : rewrap data to display the time-series epochs and interferograms
  --displacement: show displacement, instead of phase
  --opposite    : opposite sign - multiply data by -1
  --fliplr      : flip left-right
  --flipud      : flip up-down
  -x            : subset in x direction 
  -y            : subset in y direction
  -l            : subset in latitude
  -L            : subset in longitude
  --no-multilook: do not multilook for big data (by default, multilook applied for big data display)
                  use this option when high quality figure is needed.

  Point and Line:
  --point       : point coordinate in x,y,x,y,...
  --line        : line end coordinate x,y,x,y;x,y,x,y;...
                  p_yx = '1468,1002,1498,1024,1354,1140,1394,1174'
                  l_yx = '1468,1002,1498,1024;1354,1140,1394,1174,1560,1155'
                  view.py -f mask_all.h5 --point=p_yx --line=l_yx

  Figure Setting:
  --figsize     : figure size in inches (width, length), i.e. '15,10'
  -r            : row    number of figures in each window (used only when display the time-series or interferograms. Default: 5)
  -p            : column number of figures in each window (used only when display the time-series or interferograms. Default: 8)
  -i            : width  space between subplots (default 0.1)
  -j            : height space between subplots (default 0.1)
  -s            : font size (default: 12 for plot_one, 8 for plot_all)
  -c            : colormaps of matplotlib found in (http://matplotlib.org/examples/pylab_examples/show_colormaps.html)
                  (default is jet). some options are: seismic, bwr, spectral, jet, ... 
  --noaxis      : turn off axis display of figure 
  -T            : title of time-series epochs or interferograms. 
                  options are 'in' and 'out'. (default is out)
  -u/--unit     : unit for display
                  displacement: mm, cm, m (default)
                  velocity    : m/day, m/mon, cm/mon, cm/yr, m/yr (default)

  --display-radar : display in radar coordinates for geocoded files. 
                    By default, for geocoded file display in geo coordinate, otherwise display in radar coordinate

  Reference (in time and space):
  --noreference   : do not show reference point. By default, it will show if there is reference point in attributes.
  --ref-epoch     : reference date / epoch for timeseries / interferograms and wrapped file
  --ref-color     : color  of marker for the reference point (k,g,r,b,...)
  --ref-symbol    : symbol of marker for the reference point (s,o,^,v,p,x,'*'...)
  --ref-size      : size   of marker for the reference point (10 by default)

  Output and Display:
  --save        : save                    the figure
  --nodisplay   : save and do not display the figure
  -o/--output   : output figure name, use input file name by default
  --dpi         : save and change dpi number for output file (150 by default, 300 for normal print quality, 600 for high quality)

  Usage:  
           view.py file
           view.py -f file -m minValue -M maxValue -r 5 -p 5
           view.py -f file -t template_file
           view.py -f file -D DEM.dem
           view.py -f file --save --nodisplay
           view.py -f file -l Latsub -L Lonsub

  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  Example:
           view.py velocity.h5
           view.py filt_060924-070927.unw
           view.py SanAndreas.dem
           view.py -f filt_060924-070927.unw --displacement  --save
           view.py -f velocity.h5 -t ShikokuT417F650_690AlosA.template -u cm/yr
           view.py -f velocity.h5 -m -0.02 -M 0.02 -c bwr --fliplr
           view.py -f velocity.h5 --ref-color=r --ref-symbol=^ --ref-size=5

           view.py -f timeseries.h5
           view.py -f timeseries.h5 -d 20030502
           view.py -f timeseries.h5 -d 080411 -R 110118
           view.py -f timeseries.h5 -e 5 
           view.py -f timeseries.h5 -r 5 -p 8 -i 0.1 -j 0.1 --wrap

           view.py -f LoadedData.h5
           view.py -f LoadedData.h5 -d 070927-100217
           view.py -f LoadedData.h5 -d geo_filt_070927-100217-sim_HDR_4rlks_c10.unw
           view.py -f LoadedData.h5 -T in
           view.py -f Coherence.h5  -e 5
           view.py -f Wrapped.h5    -e 5

   Showing DEM:
           view.py -f velocity.h5 -D SanAndreas.dem
           view.py -f velocity.h5 -D SanAndreas.dem --dem-contour
           view.py -f velocity.h5 -D SanAndreas.dem --dem-contour --dem-noshade

   Display in subset:
           view.py -f velocity.h5 -x 100:600     -y 200:800
           view.py -f velocity.h5 -l 31.05:31.10 -L 130.05:130.10
           view.py -f timeseries.h5 -d 20100102      -x 100:600 -y 200:800
           view.py -f LoadedData.h5 -d 070927-100217 -x 100:600 -y 200:800

   Masking:
           view.py -f Seeded_LoadedData.h5 -d 931018-950809 --mask Mask_tempCoh.h5

   Showing reference:
           view.py -f velocity.h5 --noreference

   Save and Output:
           view.py -f velocity.h5 --save
           view.py -f velocity.h5 -o velocity.pdf
           view.py -f velocity.h5 --nodisplay

*****************************************************************************************
  '''


#########################################################################################
##################################  Main Function  ######################################
def main(argv):

  #################  default values  ################
  contour_step  = 200.0
  contour_sigma = 3.0
  demContour    = 'no'
  demShade      = 'yes'
  disp_axis     = 'yes'
  disp_geo      = 'yes'
  dispDisplacement = 'no'
  dispFig        = 'yes'
  dispOne        = 'yes'
  dispOpposite   = 'no'
  fig_dpi     = 100            # 150 for display, 300 for print
  figNameExt  = '.png'         # emf, eps, pdf, png, ps, raw, rgba, svg, svgz.
  #fig_rows   = 5
  #fig_cols   = 8

  #fig_size   = [15.0,8.0]     # in inches; [15.0,8.0] for 13 inch Mac;
  #fig_size   = [30.0,16.0]     # in inches; [25,15] for 28 inch Monitor;
  #flip_lr    = 'no'
  #flip_ud    = 'no'
  #font_size = 8
  Hspace     = 0.1
  Wspace     = 0.1
  masking    = 'no'
  multilook  = 'yes'
  ref_color  = 'k'
  ref_symbol = 's'
  ref_size   = 10
  rewrapping = 'no'
  saveFig    = 'no'
  showRef    = 'yes'
  title      = 'out'

  ###################  Read Input Args  ###############
  if len(sys.argv)>2:
     try:
        opts, args = getopt.getopt(argv,'c:d:D:e:E:f:h:i:j:l:L:m:M:o:p:r:s:t:T:u:x:y:',\
                                       ['help','wrap','displacement','opposite','fliplr','flipud','save','unit=',\
                                        'scale=','nodisplay','noreference','figsize=','dem-contour','dem-noshade',\
                                        'contour-step=','contour-smooth=','ref-epoch=','ref-color=','ref-symbol=',\
                                        'ref-size=','display-radar','title=','dpi=','output=','exclude=','noaxis',\
                                        'point=','line=','no-multilook','mask='])

     except getopt.GetoptError:
        print 'Error in reading input options!';  Usage() ; sys.exit(1)
     if opts==[]: Usage() ; sys.exit(1)

     for opt,arg in opts:
        if opt in ("-h","--help"):    Usage() ; sys.exit()
        elif opt == '-f': File       = arg
        elif opt == '-c': color_map  = arg
        elif opt == '-D': demFile    = arg
        elif opt == '-d': epoch_date   = [i        for i in arg.split(',')];   epoch_date.sort()
        elif opt == '-e': epoch_number = [int(i)-1 for i in arg.split(',')];   epoch_number.sort()   ## start from 1
        elif opt == '-i': Wspace     = float(arg)
        elif opt == '-j': Hspace     = float(arg)
        elif opt == '-m': disp_min = float(arg);         rewrapping='no'
        elif opt == '-M': disp_max = float(arg);         rewrapping='no'
        elif opt == '-r': fig_rows  = int(arg)
        elif opt == '-p': fig_cols  = int(arg)
        elif opt == '-s': font_size  = int(arg)
        elif opt == '-t': templateFile = arg
        elif opt == '-l': win_lat = [float(i) for i in arg.split(':')];  win_lat.sort()
        elif opt == '-L': win_lon = [float(i) for i in arg.split(':')];  win_lon.sort()
        elif opt == '-x': win_x = [int(i) for i in arg.split(':')];      win_x.sort()
        elif opt == '-y': win_y = [int(i) for i in arg.split(':')];      win_y.sort()

        elif opt in ['-E','--exclude']: exclude_epoch = [i for i in arg.split(',')];   exclude_epoch.sort()
        elif opt in ['-o','--output'] : figName      = arg;    saveFig = 'yes'
        elif opt in ['-T','--title']  : title        = arg
        elif opt in ['-u','--unit']   : disp_unit    = arg.lower()
        elif opt == '--contour-step'  : contour_step = float(arg)
        elif opt == '--contour-smooth': contour_sigma = float(arg)
        elif opt == '--dem-contour'   : demContour   = 'yes'
        elif opt == '--dem-noshade'   : demShade     = 'no'
        elif opt == '--displacement'  : dispDisplacement = 'yes';  rewrapping = 'no'
        elif opt == '--display-radar' : disp_geo = 'no'
        elif opt == '--dpi'           : fig_dpi  = int(arg);   saveFig = 'yes'
        elif opt == '--figsize'       : fig_size = [float(i) for i in arg.split(',')][0:2]
        elif opt == '--fliplr'        : flip_lr      = 'yes'
        elif opt == '--flipud'        : flip_ud      = 'yes'
        elif opt == '--mask'          : maskFile     = arg
        elif opt == '--noaxis'        : disp_axis    = 'no'
        elif opt == '--nodisplay'     : dispFig      = 'no';       saveFig = 'yes'
        elif opt == '--noreference'   : showRef      = 'no'
        elif opt == '--opposite'      : dispOpposite = 'yes'
        elif opt == '--ref-epoch'     : ref_epoch    = arg
        elif opt == '--ref-color'     : ref_color    = arg
        elif opt == '--ref-symbol'    : ref_symbol   = arg
        elif opt == '--ref-size'      : ref_size     = int(arg)
        elif opt == '--save'          : saveFig      = 'yes'
        elif opt == '--scale'         : disp_scale   = float(arg)
        elif opt == '--wrap'          : rewrapping   = 'yes'
        elif opt == '--point'         : point_yx     = [i for i in arg.split(',')];
        elif opt == '--line'          : line_yx    = [i for i in arg.split(';')];
        elif opt == '--no-multilook'  : multilook    = 'no'

  elif len(sys.argv)==2:
     if argv[0] in ['-h','--help']:              Usage(); sys.exit(1)
     elif os.path.isfile(argv[0]):  File = argv[0]
     else:    print 'Input file does not existed: '+argv[0];  sys.exit(1)
  elif len(sys.argv)<2:             Usage(); sys.exit(1)

  ##### Read File Info / Attributes
  try: atr = readfile.read_attributes(File)
  except: print 'Can not read file: '+File; sys.exit(1)
  ext = os.path.splitext(File)[1].lower()
  print '\n******************** Display ********************'
  print 'Input file is '+atr['PROCESSOR']+' '+atr['FILE_TYPE']+': '+File
  k = atr['FILE_TYPE']

  ##################  Color Map  ######################
  try: color_map
  except:
      if k in ['coherence','temporal_coherence','.cor','.dem','.hgt']:
            color_map = 'gray'
      else: color_map = 'jet'

  if color_map == 'hsv':
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
  else:  ccmap=plt.get_cmap(color_map)

  ##### Check subset range
  width  = int(atr['WIDTH'])
  length = int(atr['FILE_LENGTH'])
  print 'file size: '+str(length)+', '+str(width)

  try: win_y = subset.coord_geo2radar(win_lat,atr,'latitude')
  except:
      try:    win_y
      except: win_y = [0,length]
  try: win_x = subset.coord_geo2radar(win_lon,atr,'longitude')
  except:
      try:    win_x
      except: win_x = [0,width]
  
  win_y,win_x = subset.check_subset_range(win_y,win_x,atr)
  box = (win_x[0],win_y[0],win_x[1],win_y[1])
  if win_y[1]-win_y[0] == length and win_x[1]-win_x[0] == width:
        subsetData = 'no'
  else: subsetData = 'yes'

  ## Geo coordinate
  try:
      lon_step = float(atr['X_STEP'])
      lat_step = float(atr['Y_STEP'])
      lon_unit = atr['Y_UNIT']
      lat_unit = atr['X_UNIT']
      ullon     = float(atr['X_FIRST'])+win_x[0]*lon_step
      ullat     = float(atr['Y_FIRST'])+win_y[0]*lat_step
      llcrnrlon = ullon
      llcrnrlat = ullat+lat_step*(win_y[1]-win_y[0])
      urcrnrlon = ullon+lon_step*(win_x[1]-win_x[0])
      urcrnrlat = ullat
      geocoord='yes'
      print 'Input file is Geocoded'
      geo_box = (ullon,ullat,urcrnrlon,llcrnrlat)
  except:  geocoord='no'

  ##### Template File
  try:
      templateFile
      templateContents = readfile.read_template(templateFile)
      print 'reading pysar.view.* option from template'
  except: pass
  try:        fig_rows
  except:
      try:    fig_rows = int(templateContents['pysar.view.row'])
      except: fig_rows = 5
  try:        fig_cols
  except:
      try:    fig_cols = int(templateContents['pysar.view.column'])
      except: fig_cols = 8
  try:        exclude_epoch
  except:
      try:    exclude_epoch = templateContents['pysar.drop.date'].replace(' ','').split(',')
      except: pass
  try:        disp_unit
  except:
      try:    disp_unit = templateContents['pysar.view.unit']
      except: pass
  try: saveFig = templateContents['pysar.view.save']
  except: pass
  try: font_size = int(templateContents['pysar.view.fontSize'])
  except: pass


  ##### show displacement instead of phase
  if k in ['interferograms','.unw'] and dispDisplacement == 'yes':
      print 'show displacement'
      phase2range = -float(atr['WAVELENGTH']) / (4*np.pi)
      try:    disp_unit
      except: disp_unit = 'm'
  else: dispDisplacement = 'no'

  ## rewrapping
  if rewrapping == 'yes':
      if k in ['velocity','coherence','wrapped','temporal_coherence','mask','rmse','dem',\
               '.dem','.hgt','.slc','.mli','.trans','.cor']:
          rewrapping == 'no'
          print 'Rewrapping is disabled for '+k
      else:
          disp_min  = -np.pi
          disp_max  =  np.pi
          disp_unit = 'radian'
          print 'rewrap to -pi/pi'

  ## Scale Factor
  try:        disp_scale
  except:
      try:    disp_unit,disp_scale = unit_and_scale(atr['UNIT'],disp_unit)
      except: disp_scale = 1; disp_unit = atr['UNIT']
  print 'display in unit: '+disp_unit

  ## Flip
  try:
      orb_dir = orbit_direction(atr)
      print orb_dir+' orbit'
      ## flip by default if in radar coord
      try: flip_lr
      except:
          if orb_dir == 'descending' and geocoord == 'no': flip_lr = 'yes'
          else:                                            flip_lr = 'no'
      try: flip_ud
      except:
          if orb_dir == 'ascending'  and geocoord == 'no': flip_ud = 'yes'
          else:                                            flip_ud = 'no'
  except:
      flip_lr = 'no'
      flip_ud = 'no'

  if flip_lr == 'yes':  print 'flip left and right'
  if flip_ud == 'yes':  print 'flip up   and down'

  ## Display Min / Max
  try:
      disp_min
      disp_max
  except:
      try:
          templateFile
          if   unit_type(disp_unit) == 'velocity':
              try:
                  lim = templateContents['pysar.view.velocityLim']
                  disp_min,disp_max = sorted([float(i) for i in lim.split(':')])
              except: pass
          elif unit_type(disp_unit) == 'length':
              try:
                  lim = templateContents['pysar.view.displacementLim']
                  disp_min,disp_max = sorted([float(i) for i in lim.split(':')])
              except: pass
          else: pass
      except:
          if k in ['coherence','temporal_coherence','.cor']:
              disp_min = 0
              disp_max = 1
          elif k in ['wrapped','.int']:
              disp_min = -np.pi
              disp_max =  np.pi
          else: pass

  ##### Input File and Date List  -  One / Multiple Display
  if k in ['interferograms','coherence','wrapped','timeseries']:
      h5file = h5py.File(File,'r')
      epochList = h5file[k].keys()
      epochList = sorted(epochList)

      ## exclude epoch
      try:
          exclude_epoch
          exclude_epoch = ut.yymmdd(exclude_epoch)
          print 'exclude dates below:'
          epochList2 = []
          for epoch in epochList: epochList2.append(epoch)
          for epoch_ex in exclude_epoch:
              for epoch in epochList2:
                  if epoch_ex in epoch:
                      epochList.remove(epoch)
                      print epoch
          if len(epochList2) == len(epochList): del exclude_epoch
      except: pass

      ## convert epoch_date to epoch_number
      try:
          epoch_date
          epoch_date = ut.yymmdd(epoch_date)
          epoch_number=[]
          for i in range(len(epoch_date)):
              for j in range(len(epochList)):
                  if epoch_date[i] in epochList[j]:  epoch_number.append(j)
          epoch_number.sort()
      except: pass

      ## judge dispOne or not
      try:
          epoch_number
          if len(epoch_number)==1:
              dispOne = 'yes'
              epoch_number = epoch_number[0]
          else: dispOne = 'no'
      except:   dispOne = 'no';   epoch_number = range(0,len(epochList))
      if dispOne == 'no':
          if len(epoch_number)==0:   print '0 epoch found!';  sys.exit(1)
          else:  print 'number of epochs to display: '+str(len(epoch_number))

      ## convert reference epoch date to reference epoch number
      try:
          ref_epoch
          for i in range(len(epochList)):
              if ref_epoch in epochList[i]:  ref_epoch_number = i
          print 'Reference date: '+epochList[ref_epoch_number]
      except: pass
  else:
      try: del ref_epoch
      except: pass

  ##### Read Input Mask File
  try:
      maskFile
      msk,msk_atr = readfile.read(maskFile)
      msk = msk[win_y[0]:win_y[1],win_x[0]:win_x[1]]
      ndx = msk == 0
      print 'masking data with: '+maskFile
      masking = 'yes'
  except:
      masking = 'no'


  ####################################################################
  ########################## Display One #############################

  if dispOne == 'yes':
    ##### Setting for One Display #####
    try:    font_size
    #except: font_size=24
    except: font_size=16
    try:    fig_size
    except: fig_size   = [12.5,8.0]
    #try:    disp_axis
    #except: disp_axis = 'yes'

    ################# Data Reading ##################
    ##### Multiple Datasets File
    if k == 'timeseries':
        ## read data for display
        dset = h5file[k].get(epochList[epoch_number])
        data = dset[win_y[0]:win_y[1],win_x[0]:win_x[1]]
        ## reference date
        try:     ref_date = atr['ref_date']
        except:
            try: ref_date = ptime.yyyymmdd(atr['DATE'])
            except: pass
        try:
            ref_epoch_number
            if not epochList[ref_epoch_number] == ref_date:
                ref_date = epochList[ref_epoch_number]
                ref_dset = h5file[k].get(epochList[ref_epoch_number])
                ref_data = ref_dset[win_y[0]:win_y[1],win_x[0]:win_x[1]]
                data = data - ref_data
                del ref_data
            else: print 'input reference epoch is the same as current one, no reference change.'
        except: pass

    elif k in ('interferograms','coherence','wrapped'):
        print 'Displaying: '+epochList[epoch_number]
        dset = h5file[k][epochList[epoch_number]].get(epochList[epoch_number])
        data = dset[win_y[0]:win_y[1],win_x[0]:win_x[1]]

    ##### Single Dataset File
    else:
        data,atr = readfile.read(File,box)

    ############## Data Option ##################
    ## mask
    if masking == 'yes':  data[ndx] = np.nan

    ## show displacement instead of phase
    if dispDisplacement == 'yes':
        data = data*phase2range

    ## rewrapping
    if rewrapping == 'yes':
        data = rewrap(data,atr)

    ## Opposite Sign
    if dispOpposite == 'yes':
        print 'show opposite'
        data=-1*data

    ## Scale Factor
    if not disp_scale == 1:
        data = data * disp_scale

    ## Reference Point
    try:
       yref = (float(atr['ref_lat']) - ullat)/lat_step
       xref = (float(atr['ref_lon']) - ullon)/lon_step
    except:
       try:
           xref = int(atr['ref_x']) - win_x[0]
           yref = int(atr['ref_y']) - win_y[0]
       except:  pass

    ## Min and Max for Colorbar Extend
    data_min = np.nanmin(data)
    data_max = np.nanmax(data)
    try:    disp_min
    except: disp_min = data_min
    try:    disp_max
    except: disp_max = data_max
    if   disp_min <= data_min and disp_max >= data_max: cb_extend='neither'
    elif disp_min >  data_min and disp_max >= data_max: cb_extend='min'
    elif disp_min <= data_min and disp_max <  data_max: cb_extend='max'
    else:                                               cb_extend='both'
    print 'data    range: '+str(data_min)+' - '+str(data_max)
    print 'display range: '+str(disp_min)+' - '+str(disp_max)

    ############## Read DEM ##################
    try:
       demFile
       demRsc = readfile.read_attributes(demFile)
       print 'Show topography'

       ##### Read DEM
       if demRsc['WIDTH'] == width and demRsc['FILE_LENGTH'] == length:
           dem,demRsc = readfile.read(demFile,box)
       ##### Support Different Resolution/Size DEM
       elif subsetData == 'no':
           dem,demRsc = readfile.read(demFile)
       else:
           try:
               geo_box
               dem_win_y = subset.coord_geo2radar([geo_box[1],geo_box[3]],demRsc,'latitude')
               dem_win_x = subset.coord_geo2radar([geo_box[0],geo_box[2]],demRsc,'longitude')
               dem_win_y,dem_win_x = subset.check_subset_range(dem_win_y,dem_win_x,demRsc)
               dem_box = (dem_win_x[0],dem_win_y[0],dem_win_x[1],dem_win_y[1])
               dem,demRsc = readfile.read(demFile,dem_box)
           except: print 'Can not use different size DEM file in radar coordinate.'; sys.exit

       ##### DEM extension
       if demShade == 'yes':           #DEM basemap
          print 'show shaded relief DEM'
       if demContour == 'yes':     #contour
          print 'show contour: step = '+str(contour_step)+' m'
          import scipy.ndimage as ndimage
          dem_contour=ndimage.gaussian_filter(dem,sigma=contour_sigma,order=0)
          contour_sequence=np.arange(-6000,9000,contour_step)
    except: pass

    ##################### Display #####################
    fig = plt.figure(figsize=fig_size)
    ax = fig.add_axes([0.1,0.1,0.8,0.8])

    ## Title
    if k in ['coherence','interferograms','wrapped']:
        figTitle = epochList[epoch_number]
        if 'unwCor' in File: figTitle += '_unwCor'
    elif k == 'timeseries':
        date12 = epochList[epoch_number]
        try:
            ref_date
            if dispOpposite == 'yes': date12 = date12+'_'+ref_date
            else:                     date12 = ref_date+'_'+date12
        except: pass
        print 'Displaying '+date12
        try:    processMark = '_ts'+os.path.basename(File).split('timeseries')[1].split('.h5')[0]
        except: processMark = '_'+os.path.basename(File).split('.h5')[0]
        figTitle = date12+processMark
    else:  figTitle = File
    if rewrapping == 'yes':  figTitle += '_wrap'
    if subsetData == 'yes':  figTitle += '_sub'
    plt.title(figTitle,fontsize=font_size)

    ##### Plot in Geo-coordinate: plot in map
    if geocoord == 'yes' and disp_geo == 'yes':
       print 'plot in Lat/Lon'

       ## Map Setup
       from mpl_toolkits.basemap import Basemap
       m = Basemap(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                   resolution='l', area_thresh=1., projection='cyl',suppress_ticks=False,ax=ax)

       ## Plot DEM
       try:
          demFile
          if demShade == 'yes':
             m.imshow(ut.hillshade(dem,50.0),origin='upper', cmap='gray')
          if demContour == 'yes':
             import numpy.matlib
             c_x = np.linspace(llcrnrlon,urcrnrlon,num=dem.shape[1],endpoint='FALSE').reshape(1,dem.shape[1])
             c_xx= np.matlib.repmat(c_x,dem.shape[0],1)
             c_y = np.linspace(urcrnrlat,llcrnrlat,num=dem.shape[0],endpoint='FALSE').reshape(dem.shape[0],1)
             c_yy= np.matlib.repmat(c_y,1,dem.shape[1])
             m.contour(c_xx,c_yy,dem_contour,contour_sequence,origin='upper',colors='black',alpha=0.5,latlon='FALSE')
       except:  pass

       ## Plot Data
       try:     im = m.imshow(data,cmap=ccmap,origin='upper',vmin=disp_min,vmax=disp_max)
       except:  im = m.imshow(data,cmap=ccmap,origin='upper')

       # Reference Point
       if showRef == 'yes':
          try:
             refPoint=ref_color+ref_symbol
             ref_lon = llcrnrlon + xref*lon_step
             ref_lat = urcrnrlat + yref*lat_step
             plt.plot(ref_lon,ref_lat,refPoint,ms=ref_size)
          except:  pass

       # Colorbar
       from mpl_toolkits.axes_grid1 import make_axes_locatable
       divider = make_axes_locatable(ax)
       cax = divider.append_axes("right",size="5%", pad=0.30)
       cbar = plt.colorbar(im,cax=cax,extend=cb_extend)
       cbar.set_label(disp_unit)
       #plt.colorbar(im,cax=cax)

       # Status bar
       def format_coord(x,y):
         col = int((x-ullon)/lon_step+0.5)
         row = int((y-ullat)/lat_step+0.5)
         if col>=0 and col<=data.shape[1] and row >=0 and row<=data.shape[0]:
            z = data[row,col]
            try:
                h = dem[row,col]
                return 'lon=%.4f,  lat=%.4f,  elev=%.1f m,  value=%.4f'%(x,y,h,z)
            except:
                return 'lon=%.4f,  lat=%.4f,  value=%.4f'%(x,y,z)
         else:
            return 'lon=%.4f,  lat=%.4f'%(x,y)
       ax.format_coord = format_coord


    ##### Plot in x/y coordinate: row and column
    else:
       print 'plot in Y/X'

       ## Plot DEM
       try:
          demFile
          if demShade == 'yes':
             ax.imshow(ut.hillshade(dem,50.0), cmap='gray')
          if demContour == 'yes':
             ax.contour(dem_contour,contour_sequence,origin='lower',colors='black',alpha=0.5)
       except:  pass

       ## Plot Data
       try:     im = ax.imshow(data,cmap=ccmap, vmin=disp_min, vmax=disp_max)
       except:  im = ax.imshow(data,cmap=ccmap)

       ## Colorbar
       cbar = plt.colorbar(im,extend=cb_extend)
       cbar.set_label(disp_unit)

       ## Reference Point
       if showRef == 'yes':
          try:
             refPoint=ref_color+ref_symbol
             ax.plot(xref,yref,refPoint,ms=ref_size)
          except:  pass

       plt.xlim(0,np.shape(data)[1])
       plt.ylim(  np.shape(data)[0],0)

       ##### Plot Points and Lines
       try:
           point_yx
           point_num = len(point_yx)/2*2
           point_yx = point_yx[0:point_num]
           point_x = point_yx[1::2]
           point_y = point_yx[::2]
           plt.plot(point_x,point_y,'ro')
           print 'plot points'
       except: pass

       try:
           line_yx
           for i in range(0,len(line_yx)):
               line_xx = line_yx[i].split(',')[1::2]
               line_yy = line_yx[i].split(',')[::2]
               plt.plot(line_xx,line_yy,'r',lw=2)
           print 'plot lines'
       except: pass

       ## Status bar
       def format_coord(x,y):
         col = int(x+0.5)
         row = int(y+0.5)
         if col>=0 and col<=data.shape[1] and row >=0 and row<=data.shape[0]:
            z = data[row,col]
            try:
                h = dem[row,col]
                return 'x=%.4f,  y=%.4f,  elev=%.1f m,  value=%.4f'%(x,y,h,z)
            except:
                return 'x=%.4f,  y=%.4f,  value=%.4f'%(x,y,z)
         else:
            return 'x=%.4f,  y=%.4f'%(x,y)
       ax.format_coord = format_coord

    ##### Figure Setting
    ## Flip
    if flip_lr == 'yes':  fig.gca().invert_xaxis()
    if flip_ud == 'yes':  fig.gca().invert_yaxis()
    ## Turn off axis
    if disp_axis == 'no': ax.axis('off')

    ##### Save Figure
    if saveFig == 'yes':
       try:
          figName
       except:
          figNameBase = figTitle
          figName = figNameBase+figNameExt

       plt.savefig(figName,bbox_inches='tight',transparent=True,dpi=fig_dpi)
       print 'Saved figure to '+figName

    ###### Show Figure
    if dispFig == 'yes':
       plt.show()
    

  ####################################################################
  ########################## Display Multiple ########################  

  elif dispOne == 'no':
    ##### Setting for Multiple Display #####
    print 'row    number: '+str(fig_rows)
    print 'column number: '+str(fig_cols)
    try:    font_size
    except: font_size=12
    try:    fig_size
    except: fig_size   = [30.0,16.0]
    #try:    disp_axis
    #except: disp_axis = 'no'

    ## Figure Name
    if saveFig == 'yes':
        try:
            figName
            figNameExt  = os.path.splitext(figName)[1].lower()
            figNameBase = os.path.basename(figName).split(figNameExt)[0]
        except:
            figNameBase = os.path.basename(File).split(ext)[0]
        try:
            exclude_epoch
            figNameBase += '_ex'
        except: pass
        if rewrapping == 'yes': figNameBase += '_wrap'
        if subsetData == 'yes': figNameBase += '_sub'

    ## Reference date for timeseries
    if k == 'timeseries':
        try:     ref_date = atr['ref_date']
        except:
            try: ref_date = ptime.yyyymmdd(atr['DATE'])
            except: pass
        try:
            ref_epoch_number
            if not epochList[ref_epoch_number] == ref_date:
                ref_date = epochList[ref_epoch_number]
                ref_dset = h5file[k].get(epochList[ref_epoch_number])
                ref_data = ref_dset[win_y[0]:win_y[1],win_x[0]:win_x[1]]
                figNameBase = figNameBase+'_ref'+ref_date
            else: print 'input reference epoch is the same as current one, no reference change.'
        except: pass

    ##### Multilook if too many subplots in one figure
    ## for less memory and faster speed
    nfigs = fig_rows*fig_cols                                  ## number of subplots per figure
    lks = 1
    if multilook == 'yes':
        win_size = (win_y[1]-win_x[0])*(win_x[1]-win_x[0])
        if   win_size * nfigs > (8e6*160):   lks=16;       ## 2k * 2k image with 120 subplots
        elif win_size * nfigs > (4e6*80) :   lks=8;       ## 2k * 2k image with 80  subplots
        elif win_size * nfigs > (4e6*20) :   lks=4;       ## 2k * 2k image with 40  subplots
        elif win_size * nfigs > (1e6*20) :   lks=2;       ## 2k * 2k image with 40  subplots
    if lks > 1:
        from pysar.multi_looking import multilook         
        print 'number of data points per figure: '+'%.1E' %(win_size*nfigs)+', multilook with a factor of '+str(lks)

    ################## DEM Options ####################
    try:
        ## Read DEM
        dem,demRsc = readfile.read(demFile,box)
        print 'Show topography'
        if lks > 1:   dem = multilook(dem,lks,lks)

        ## DEM Extension
        if demShade == 'yes':           #DEM basemap
            print 'plot DEM as basemap'
            cmap_dem=plt.get_cmap('gray')
            hillshade_dem=ut.hillshade(dem,50.0)
        if demContour == 'yes':     #contour
            print 'plot contour: step = '+str(contour_step)+' m'
            import scipy.ndimage as ndimage
            dem=ndimage.gaussian_filter(dem,sigma=contour_sigma,order=0)
            contour_sequence=np.arange(-6000,9000,contour_step)
    except:  pass


    ################## Plot Loop ####################
    #epochList=h5file[k].keys()
    nepoch = len(epoch_number)
    fig_num = int(float(nepoch)/float(nfigs) - 1e-5) + 1
    print 'figure number: '+str(fig_num)

    ## Find min and value for all data, reference for better min/max setting next time
    all_data_min=0
    all_data_max=0

    ##### Loop - Figures
    for j in range(1,fig_num+1):
        print '----------------------------------------'
        print 'figure '+str(j)
        fig = plt.figure(j,figsize=fig_size)
        fig_data_min=0
        fig_data_max=0

        ## Starting Epoch and Ending Epoch in this figure
        i_start = (j-1)*nfigs
        i_end   = min([nepoch,i_start+nfigs])

        ##### Loop - Subplots
        for i in range(i_start,i_end):
            epoch = epochList[epoch_number[i]]
            print 'loading '+epoch+'  '+str(i+1)
            ax = fig.add_subplot(fig_rows,fig_cols,i-i_start+1) 

            ##### Data Reading
            if k == 'timeseries':
                figTitle = epoch[0:4]+'-'+epoch[4:6]+'-'+epoch[6:8]
                dset = h5file[k].get(epoch)
                data = dset[win_y[0]:win_y[1],win_x[0]:win_x[1]]
                try: data = data - ref_data
                except: pass
            elif k in ('interferograms','coherence','wrapped'):
                if   nfigs > 100:   figTitle = str(epoch_number[i]+1)
                elif nfigs > 50 :   figTitle = str(epoch_number[i]+1)+'\n' +h5file[k][epoch].attrs['DATE12']
                else:               figTitle = str(epoch_number[i]+1)+' : '+h5file[k][epoch].attrs['DATE12']
                dset = h5file[k][epoch].get(epoch)
                data = dset[win_y[0]:win_y[1],win_x[0]:win_x[1]]

            ##### Data Option
            ## mask file
            if masking          == 'yes':   data[ndx] = np.nan
            ## multilooking
            if lks              >  1    :   data = multilook(data,lks,lks)
            ## show displacement instead of phase
            if dispDisplacement == 'yes':   data = data*phase2range
            ## rewrapping
            if rewrapping       == 'yes':   data = rewrap(data,atr)
            ## Opposite Sign
            if dispOpposite     == 'yes':   data = -1*data
            ## Scale Factor
            if not disp_scale   == 1:       data = data*disp_scale
            ## Data Min/Max
            fig_data_min = np.nanmin([fig_data_min,np.nanmin(data)])
            fig_data_max = np.nanmax([fig_data_max,np.nanmax(data)])

            # Plot
            try:
                demFile
                if demShade   == 'yes':  plt.imshow(hillshade_dem,cmap='gray')
                if demContour == 'yes':  plt.contour(dem,contour_sequence,origin='lower',colors='black',alpha=0.5)
            except:  pass
            try:     im = ax.imshow(data,cmap=ccmap,vmin=disp_min,vmax=disp_max)
            except:  im = ax.imshow(data,cmap=ccmap)
            #del data

            ax.set_yticklabels([])
            ax.set_xticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            if   title == 'out':  ax.set_title(figTitle,fontsize=font_size)
            elif title == 'in':   add_inner_title(ax, figTitle, loc=1)

            ## Flip
            if flip_lr == 'yes':  fig.gca().invert_xaxis()
            if flip_ud == 'yes':  fig.gca().invert_yaxis()
            ## Turn off axis
            if disp_axis == 'no': ax.axis('off')

        ## Min and Max for this figure
        all_data_min = np.nanmin([all_data_min,fig_data_min])
        all_data_max = np.nanmax([all_data_max,fig_data_max])
        print 'data    range: '+str(fig_data_min)+' - '+str(fig_data_max)
        try:  print 'display range: '+str(disp_min)+' - '+str(disp_max)
        except: pass

        ##### Colorbar
        try:
            disp_min
            disp_max
            if   disp_min <= fig_data_min and disp_max >= fig_data_max: cb_extend='neither'
            elif disp_min >  fig_data_min and disp_max >= fig_data_max: cb_extend='min'
            elif disp_min <= fig_data_min and disp_max <  fig_data_max: cb_extend='max'
            else:                                                       cb_extend='both'
            fig.subplots_adjust(wspace=Wspace,hspace=Hspace, right=0.9)
            cbar_ax = fig.add_axes([0.91, 0.15, 0.015, 0.7])
            cbar = fig.colorbar(im, cax=cbar_ax, extend=cb_extend)
            cbar.set_label(disp_unit)
        except: print 'Different color scale for each subplot!'

        ## Save Figure
        if saveFig == 'yes':
           if fig_num > 1:  figName = figNameBase+'_'+str(j)+figNameExt
           else:            figName = figNameBase+figNameExt

           plt.savefig(figName,bbox_inches='tight',transparent=True,dpi=fig_dpi)
           print 'saved figure to '+figName
           if dispFig == 'no':  fig.clf()


    print '----------------------------------------'
    print 'all data range: '+str(all_data_min)+' - '+str(all_data_max)
    try: print 'display  range: '+str(disp_min)+' - '+str(disp_max)
    except: pass

    ## Display Figure
    if dispFig == 'yes':    plt.show()
   
  ##########################################################
  try: h5file.close()
  except: pass


#########################################################################################
#########################################################################################
if __name__ == '__main__':
  main(sys.argv[1:])


