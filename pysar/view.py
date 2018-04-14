#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2017, Zhang Yunjun, Heresh Fattahi          #
# Author:  Zhang Yunjun, Heresh Fattahi                    #
############################################################
#
# Recommended usage:
#   import pysar.view as pv
#

import os
import sys
import argparse
from datetime import datetime as dt

import h5py
import numpy as np
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LinearSegmentedColormap, LightSource
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import cm, pyproj

import pysar
import pysar.utils.datetime as ptime
import pysar.utils.readfile as readfile
import pysar.utils.utils as ut
import pysar.utils.plot as pp

import pysar.mask as mask
import pysar.subset as subset
from pysar.multilook import multilook_matrix

from pysar.utils.readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file
from pysar.utils.plot import Basemap2

fig = None

##################################################################################################
def auto_figure_title(fname, epoch=[], inps_dict=None):
    '''Get auto figure title from meta dict and input options
    Inputs:
        fname - string, input file name
        epoch - list of string, optional, epoch to read for multi dataset/group files
        inps_dict - dict, optional, processing attributes, including:
                    ref_date
                    pix_box
                    wrap
                    disp_scale
                    opposite
    Output:
        fig_title - string, output figure title
    Example:
        'geo_velocity.h5' = auto_figure_title('geo_velocity.h5', None, vars(inps))
        '101020-110220_ECMWF_demErr_quadratic' = auto_figure_title('timeseries_ECMWF_demErr_quadratic.h5', '110220')
    '''
    atr = readfile.read_attribute(fname)
    k = atr['FILE_TYPE']
    width = int(atr['WIDTH'])
    length = int(atr['LENGTH'])

    #if not epoch and k in multi_group_hdf5_file+multi_dataset_hdf5_file:
    #    print "No date/date12 input.\nIt's required for "+k+" file\nReturn None"
    #    return None
    if not epoch:
        epoch = []

    if len(epoch)==1 and k in multi_group_hdf5_file:
        fig_title = epoch[0]
        if 'unwCor' in fname:
            fig_title += '_unwCor'

    elif len(epoch)==1 and k in multi_dataset_hdf5_file+['GIANT_TS']:
        if inps_dict['REF_DATE']:
            ref_date = inps_dict['REF_DATE']
        else:
            try:
                ref_date = atr['REF_DATE']
            except:
                ref_date = None
        if not ref_date:
            fig_title = epoch[0]
        else:
            fig_title = ptime.yymmdd(ref_date)+'-'+ptime.yymmdd(epoch[0])

        try:
            ext = os.path.splitext(fname)[1]
            processMark = os.path.basename(fname).split('timeseries')[1].split(ext)[0]
            fig_title += processMark
        except: pass
    else:
        fig_title = os.path.splitext(os.path.basename(fname))[0]

    # mark - subset
    try:
        pix_box = inps_dict['pix_box']
        if (pix_box[2]-pix_box[0])*(pix_box[3]-pix_box[1]) < width*length:
            fig_title += '_sub'
    except: pass

    # mark - rewrapping
    try:
        rewrapping = inps_dict['wrap']
        if rewrapping:
            fig_title += '_wrap'
    except: pass

    # mark - scale
    try:
        scaling = inps_dict['disp_scale']
        if not scaling == 1.0:
            fig_title += '_scale'+str(scaling)
    except: pass

    # mark - opposite
    try:
        disp_opposite = inps_dict['opposite']
        if disp_opposite:
            fig_title += '_oppo'
    except: pass

    return fig_title


##################################################################################################
def check_multilook_input(pixel_box, row_num, col_num):
    # Estimate multilook_num
    box_size = (pixel_box[2]-pixel_box[0])*(pixel_box[3]-pixel_box[1])
    pixel_num_per_figure = box_size*row_num*col_num
    if   pixel_num_per_figure > (8e6*160):   multilook_num=16;      ## 2k * 2k image with 120 subplots
    elif pixel_num_per_figure > (4e6*80) :   multilook_num=8;       ## 2k * 2k image with 80  subplots
    elif pixel_num_per_figure > (4e6*20) :   multilook_num=4;       ## 2k * 2k image with 40  subplots
    elif pixel_num_per_figure > (1e6*20) :   multilook_num=2;       ## 2k * 2k image with 40  subplots
    else: multilook_num=1
    # Update multilook based on multilook_num
    if multilook_num > 1:
        multilook = True
        print('number of data points per figure: '+'%.1E' %(pixel_num_per_figure))
        print('multilook with a factor of '+str(multilook_num)+' for display')
    else:
        multilook = False
    return multilook, multilook_num


##################################################################################################
def get_epoch_full_list_from_input(all_epoch_list, epoch_input_list=[], epoch_num_input_list=[]):
    '''Read/Get input epoch list from input epoch and epoch_num'''
    epoch_input_list = sorted(epoch_input_list)
    epoch_num_input_list = sorted([i-1 for i in epoch_num_input_list])

    # Default value
    epoch_list = []
    epoch_num_list = []

    # Update epoch_num_list from epoch_input_list
    if epoch_input_list:
        for epoch_input in epoch_input_list:
            epoch_list_tmp = [s for s in all_epoch_list if epoch_input.lower() in s.lower()]
            epoch_list += list(set(epoch_list_tmp) - set(epoch_list))
        epoch_num_list += [all_epoch_list.index(s) for s in epoch_list]

    # Update epoch_num_list from epoch_num_input_list
    epoch_num_list += list(set(epoch_num_input_list) - set(epoch_num_list))

    # Convert epoch_num_list to epoch_list
    epoch_list = [all_epoch_list[i] for i in epoch_num_list]

    return epoch_list, epoch_num_list


##################################################################################################
def scale_data4disp_unit_and_rewrap(data, atr, disp_unit=None, rewrapping=False):
    '''Scale 2D matrix value according to display unit and re-wrapping flag
    Disable rewrapping option 1) for specific data types, which rewrapping has no physical meaning;
                              2) if disp_unit exists and != 'radian'; priority: disp_unit > rewrapping 
    Inputs:
        data - 2D np.array
        atr  - dict, including the following attributes:
               UNIT
               FILE_TYPE
               WAVELENGTH
        disp_unit  - string, optional
        rewrapping - bool, optional
    Outputs:
        data
        disp_unit
        rewrapping
    '''

    # Check re-wrap's conflict with disp_unit
    k = atr['FILE_TYPE']
    if not disp_unit and rewrapping:
        if k not in ['coherence','temporal_coherence','mask', 'dem', '.dem','.hgt',\
                     '.slc','.mli','.trans','.cor']:
            disp_unit = 'radian'
        else:
            rewrapping = False
            print('WARNING: re-wrap is disabled for '+k)
    elif disp_unit not in ['radian'] and rewrapping:
        print('WARNING: re-wrap is disabled because display unit is not radian.')
        rewrapping = False

    # Default unit
    if not disp_unit:
        if k in ['.mli','.slc','.amp']:
            disp_unit = 'dB'
        else:
            disp_unit = atr['UNIT']

    # Data Operation - Scale to display unit
    if not disp_unit == atr['UNIT']:
        data, disp_unit = scale_data2disp_unit(data, atr, disp_unit)        
    print('display in unit: '+disp_unit)

    # Data Operation - Rewrapping
    if rewrapping and 'radian' in disp_unit:
        print('re-wrapping data to [-pi, pi]')
        data -= np.round(data/(2*np.pi)) * (2*np.pi)

    return data, disp_unit, rewrapping


def scale_data2disp_unit(matrix, atr_dict, disp_unit):
    '''Scale data based on data unit and display unit
    Inputs:
        matrix    : 2D np.array
        atr_dict  : dictionary, meta data
        disp_unit : str, display unit
    Outputs:
        matrix    : 2D np.array, data after scaling
        disp_unit : str, display unit
    Default data file units in PySAR are:  m, m/yr, radian, 1
    '''

    # Initial 
    scale = 1.0
    data_unit = atr_dict['UNIT'].lower().split('/')
    disp_unit = disp_unit.lower().split('/')

    # if data and display unit is the same
    if disp_unit == data_unit:
        return matrix, disp_unit

    # Calculate scaling factor  - 1
    # phase unit - length / angle 
    if data_unit[0].endswith('m'):
        if   disp_unit[0] == 'mm': scale *= 1000.0
        elif disp_unit[0] == 'cm': scale *= 100.0
        elif disp_unit[0] == 'dm': scale *= 10.0
        elif disp_unit[0] == 'km': scale *= 1/1000.0
        elif disp_unit[0] in ['radians','radian','rad','r']:
            range2phase = -(4*np.pi) / float(atr_dict['WAVELENGTH'])
            scale *= range2phase
        else:
            print('Unrecognized display phase/length unit: '+disp_unit[0])
            return
        
        if   data_unit[0] == 'mm': scale *= 0.001
        elif data_unit[0] == 'cm': scale *= 0.01
        elif data_unit[0] == 'dm': scale *= 0.1
        elif data_unit[0] == 'km': scale *= 1000.
        
    elif data_unit[0] == 'radian':
        phase2range = -float(atr_dict['WAVELENGTH']) / (4*np.pi)
        if   disp_unit[0] == 'mm': scale *= phase2range * 1000.0
        elif disp_unit[0] == 'cm': scale *= phase2range * 100.0
        elif disp_unit[0] == 'dm': scale *= phase2range * 10.0
        elif disp_unit[0] == 'km': scale *= phase2range * 1/1000.0
        elif disp_unit[0] in ['radians','radian','rad','r']:
            pass
        else:
            print('Unrecognized phase/length unit: '+disp_unit[0])
            return

    # amplitude/coherence unit - 1
    elif data_unit[0] == '1':
        if disp_unit[0] == 'db':
            ind = np.nonzero(matrix)
            matrix[ind] = 10*np.log10(np.absolute(matrix[ind]))
            disp_unit[0] = 'dB'
        else:
            try:
                scale /= float(disp_unit[0])
            except:
                print('Un-scalable display unit: '+disp_unit[0])
    else:
        print('Un-scalable data unit: '+data_unit)

    # Calculate scaling factor  - 2
    if len(data_unit)==2:
        try:
            disp_unit[1]
            if   disp_unit[1] in ['y','yr','year'  ]: disp_unit[1] = 'yr'
            elif disp_unit[1] in ['m','mon','month']: disp_unit[1] = 'mon'; scale *= 12.0
            elif disp_unit[1] in ['d','day'        ]: disp_unit[1] = 'day'; scale *= 365.25
            else: print('Unrecognized time unit for display: '+disp_unit[1])
        except:
            disp_unit.append('yr')
        disp_unit = disp_unit[0]+'/'+disp_unit[1]
    else:
        disp_unit = disp_unit[0]

    # Scale input matrix
    matrix *= scale

    return matrix, disp_unit


##################################################################################################
def update_plot_inps_with_display_setting_file(inps, disp_set_file):
    '''Update inps using values from display setting file'''
    disp_set_dict = readfile.read_template(disp_set_file)
    if not inps.disp_unit and 'plot.displayUnit' in disp_set_dict.keys():
        inps.disp_unit = disp_set_dict['plot.displayUnit']
    if not inps.disp_min and 'plot.displayMin' in disp_set_dict.keys():
        inps.disp_min = float(disp_set_dict['plot.displayMin'])
    if not inps.disp_max and 'plot.displayMax' in disp_set_dict.keys():
        inps.disp_max = float(disp_set_dict['plot.displayMax'])
     
    if not inps.colormap and 'plot.colormap' in disp_set_dict.keys():
        inps.colormap = disp_set_dict['plot.colormap']

    if not inps.subset_lat and 'plot.subset.lalo' in disp_set_dict.keys():
        inps.subset_lat = [float(n) for n in disp_set_dict['plot.subset.lalo'].replace(',',' ').split()[0:2]]
    if not inps.subset_lon and 'plot.subset.lalo' in disp_set_dict.keys():
        inps.subset_lon = [float(n) for n in disp_set_dict['plot.subset.lalo'].replace(',',' ').split()[2:4]]
    if not inps.seed_lalo and 'plot.seed.lalo' in disp_set_dict.keys():
        inps.seed_lalo = [float(n) for n in disp_set_dict['plot.referenceLalo'].replace(',',' ').split()]
    
    return inps

def update_plot_inps_with_meta_dict(inps, meta_dict):
    k = meta_dict['FILE_TYPE']
    width = int(float(meta_dict['WIDTH']))
    length = int(float(meta_dict['LENGTH']))

    # default mask file:
    if not inps.mask_file and k in ['velocity','timeseries']:
        if os.path.basename(meta_dict['FILE_PATH']).startswith('geo_'):
            inps.mask_file = 'geo_maskTempCoh.h5'
        else:
            inps.mask_file = 'maskTempCoh.h5'
        if not os.path.isfile(inps.mask_file):
            inps.mask_file = None

    # Subset
    ## Convert subset input into bounding box in radar / geo coordinate
    ## geo_box = None if atr is not geocoded. 
    inps.pix_box, inps.geo_box = subset.subset_input_dict2box(vars(inps), meta_dict)
    inps.pix_box = subset.check_box_within_data_coverage(inps.pix_box, meta_dict)
    inps.geo_box = subset.box_pixel2geo(inps.pix_box, meta_dict)
    # Out message
    data_box = (0,0,width,length)
    print('data   coverage in y/x: '+str(data_box))
    print('subset coverage in y/x: '+str(inps.pix_box))
    print('data   coverage in lat/lon: '+str(subset.box_pixel2geo(data_box, meta_dict)))
    print('subset coverage in lat/lon: '+str(inps.geo_box))
    print('------------------------------------------------------------------------')
    
    # Multilook
    # if too many subplots in one figure for less memory and faster speed
    if inps.multilook_num > 1:
        inps.multilook = True

    # Colormap
    inps.colormap = pp.check_colormap_input(meta_dict, inps.colormap)

    # Seed Point
    # Convert seed_lalo if existed, to seed_yx, and use seed_yx for the following
    # seed_yx is referenced to input data coverage, not subseted area for display
    if inps.seed_lalo and inps.geo_box:
        inps.seed_yx = [subset.coord_geo2radar(inps.seed_lalo[0], meta_dict, 'lat'), \
                        subset.coord_geo2radar(inps.seed_lalo[1], meta_dict, 'lon')]
        print('input reference point in lat/lon: '+str(inps.seed_lalo))
        print('input reference point in y  /x  : '+str(inps.seed_yx))

    # Unit and Wrap
    # Check re-wrap's conflict with disp_unit.  Priority: disp_unit > wrap > disp_unit(None)
    if not inps.disp_unit and inps.wrap:
        if k not in ['coherence','temporal_coherence','mask', 'dem', '.dem','.hgt',\
                     '.slc','.mli','.trans','.cor']:
            inps.disp_min  = -np.pi
            inps.disp_max  =  np.pi
            inps.disp_unit = 'radian'
            print('re-wrap data to [-pi, pi] for display')
        else:
            inps.wrap = False
            print('WARNING: re-wrap is disabled for '+k)
    elif inps.disp_unit not in ['radian'] and inps.wrap:
        print('WARNING: re-wrap is disabled because display unit is not radian.')
        inps.wrap = False

    # Default unit for amplitude image
    if not inps.disp_unit and k in ['.mli','.slc','.amp']:
        inps.disp_unit = 'dB'

    # Min / Max - Display
    if not inps.disp_min and not inps.disp_max:
        if k in ['coherence','temporal_coherence','.cor']:
            inps.disp_min = 0.0
            inps.disp_max = 1.0
        elif k in ['wrapped','.int']:
            inps.disp_min = -np.pi
            inps.disp_max =  np.pi

    # Transparency - Alpha
    if not inps.transparency:
        ## Auto adjust transparency value when showing shaded relief DEM
        if inps.dem_file and inps.disp_dem_shade:
            inps.transparency = 0.8
        else:
            inps.transparency = 1.0    

    # Flip Left-Right / Up-Down
    if not inps.flip_lr and not inps.flip_ud:
        inps.flip_lr, inps.flip_ud = pp.auto_flip_direction(meta_dict)

    # Figure Title
    if not inps.fig_title:
        try:    inps.fig_title = auto_figure_title(meta_dict['FILE_PATH'], inps.epoch, vars(inps))
        except: inps.fig_title = os.path.splitext(os.path.basename(meta_dict['FILE_PATH']))[0]
    print('figure title: '+inps.fig_title)

    # Figure output file name
    if not inps.outfile:
        inps.outfile = inps.fig_title+inps.fig_ext

    return inps


##################################################################################################
def update_matrix_with_plot_inps(data, meta_dict, inps):

    # Seed Point
    # If value of new seed point is not nan, re-seed the data and update inps.seed_yx/lalo
    # Otherwise, try to read seed info from atrributes into inps.seed_yx/lalo
    if inps.seed_yx and ('REF_Y' not in meta_dict.keys() or \
                         inps.seed_yx != [int(meta_dict['REF_Y']), int(meta_dict['REF_X'])]):
        inps.seed_value = data[inps.seed_yx[0]-inps.pix_box[1], inps.seed_yx[1]-inps.pix_box[0]]
        if not np.isnan(inps.seed_value):
            data -= inps.seed_value
            if meta_dict['FILE_TYPE'] in multi_group_hdf5_file+multi_dataset_hdf5_file:
                print('set reference point to: '+str(inps.seed_yx))
            if inps.geo_box:
                inps.seed_lalo = [subset.coord_radar2geo(inps.seed_yx[0], meta_dict, 'y'), \
                                  subset.coord_radar2geo(inps.seed_yx[1], meta_dict, 'x')]
            else:
                inps.seed_lalo = None
        else:
            print('WARNING: input reference point has nan value, continue with original reference info')
            inps.seed_yx = None
    else:
        if 'REF_Y' in meta_dict.keys():
            inps.seed_yx = [int(meta_dict['REF_Y']), int(meta_dict['REF_X'])]
        else:
            inps.seed_yx = None

        if 'REF_LAT' in meta_dict.keys():
            inps.seed_lalo = [float(meta_dict['REF_LAT']), float(meta_dict['REF_LON'])]
        elif inps.seed_yx and inps.geo_box:
            inps.seed_lalo = [subset.coord_radar2geo(inps.seed_yx[0], meta_dict, 'y'), \
                              subset.coord_radar2geo(inps.seed_yx[1], meta_dict, 'x')]
        else:
            inps.seed_lalo = None

    # Multilook
    if inps.multilook and inps.multilook_num > 1:
        data = multilook_matrix(data, inps.multilook_num, inps.multilook_num)

    # Convert data to display unit
    if not inps.disp_unit:
        inps.disp_unit = meta_dict['UNIT']
    if not inps.disp_unit == meta_dict['UNIT']:
        data, inps.disp_unit = scale_data2disp_unit(data, meta_dict, inps.disp_unit)        
    #print 'display in unit: '+inps.disp_unit

    # Re-wrap
    if inps.wrap and 'radian' in inps.disp_unit:
        #print 're-wrapping data to [-pi, pi]'
        data -= np.round(data/(2*np.pi)) * (2*np.pi)

    # 1.4 Scale 
    if not inps.disp_scale == 1.0:
        print('scaling data by a factor of '+str(inps.disp_scale))
        data *= inps.disp_scale

    # 1.5 Opposite
    if inps.opposite:
        print('show opposite')
        data *= -1
    
    return data, inps


##################################################################################################
def plot_matrix(ax, data, meta_dict, inps=None):
    '''Plot 2D matrix 
    
    Inputs:
        ax   : matplot.pyplot axes object
        data : 2D np.array, 
        meta_dict : dictionary, attributes of data
        inps : Namespace, optional, input options for display
    
    Outputs:
        ax  : matplot.pyplot axes object
    
    Example:
        import matplotlib.pyplot as plt
        import pysar._readfile as readfile
        import pysar.view as pv
         
        data, atr = readfile.read('velocity.h5')
        fig = plt.figure()
        ax = fig.add_axes([0.1,0.1,0.8,0.8])
        ax = pv.plot_matrix(ax, data, atr)
        plt.show()
    '''

    #----------------------- 0. Initial a inps Namespace if no inps input --------------------#
    if not inps:
        inps = cmdLineParse([''])
        inps = update_plot_inps_with_meta_dict(inps, meta_dict)

    #----------------------- 1.1 Update plot inps with metadata dict -------------------------#
    #inps = update_plot_inps_with_meta_dict(inps, meta_dict)

    #----------------------- 1.2 Update plot inps/data with data matrix ----------------------#
    data, inps = update_matrix_with_plot_inps(data, meta_dict, inps)
    print('data    unit: '+meta_dict['UNIT'])
    print('display unit: '+inps.disp_unit)

    # 1.6 Min / Max - Data/Display
    data_min = np.nanmin(data)
    data_max = np.nanmax(data)
    if inps.disp_min is None:  inps.disp_min = data_min
    if inps.disp_max is None:  inps.disp_max = data_max
    print('data    range: %f - %f' % (data_min, data_max))
    print('display range: %f - %f' % (inps.disp_min, inps.disp_max))

    # 1.7 DEM
    if inps.dem_file:
        dem_meta_dict = readfile.read_attribute(inps.dem_file)
        print('reading DEM: '+os.path.basename(inps.dem_file)+' ...')
        if inps.geo_box:
            # Support DEM with different Resolution and Coverage 
            inps.dem_pix_box = subset.box_geo2pixel(inps.geo_box, dem_meta_dict)
        else:
            inps.dem_pix_box = inps.pix_box
        dem, dem_meta_dict = readfile.read(inps.dem_file, inps.dem_pix_box, epoch='height')

        # If data is too large, do not show DEM contour
        if inps.geo_box:
            lat_length = abs(inps.geo_box[1]-inps.geo_box[3])
            lon_length = abs(inps.geo_box[2]-inps.geo_box[0])
            if max(lat_length, lon_length) > 1.0:
                inps.disp_dem_contour = False
                print('area is too large (lat or lon > 1 deg), turn off the DEM contour display')

    print('display data in transparency: '+str(inps.transparency))


    #-------------------- 2.1 Plot in Geo-coordinate using Basemap --------------------------------#
    if inps.geo_box and inps.fig_coord=='geo':
        print('plot in Lat/Lon coordinate')
        # Map Setup
        print('map projection: '+inps.map_projection)
        print('boundary database resolution: '+inps.resolution)
        if inps.map_projection in ['cyl','merc','mill','cea','gall']:
            m = Basemap2(llcrnrlon=inps.geo_box[0], llcrnrlat=inps.geo_box[3],\
                        urcrnrlon=inps.geo_box[2], urcrnrlat=inps.geo_box[1],\
                        projection=inps.map_projection,\
                        resolution=inps.resolution, area_thresh=1., suppress_ticks=False, ax=ax)
        elif inps.map_projection in ['ortho']:
            m = Basemap2(lon_0=(inps.geo_box[0]+inps.geo_box[2])/2.0,\
                        lat_0=(inps.geo_box[3]+inps.geo_box[1])/2.0,\
                        projection=inps.map_projection,\
                        resolution=inps.resolution, area_thresh=1., suppress_ticks=False, ax=ax)
        else:
            m = Basemap2(lon_0=(inps.geo_box[0]+inps.geo_box[2])/2.0,\
                        lat_0=(inps.geo_box[3]+inps.geo_box[1])/2.0,\
                        llcrnrlon=inps.geo_box[0], llcrnrlat=inps.geo_box[3],\
                        urcrnrlon=inps.geo_box[2], urcrnrlat=inps.geo_box[1],\
                        projection=inps.map_projection,\
                        resolution=inps.resolution, area_thresh=1., suppress_ticks=False, ax=ax)

        # Draw coastline
        if inps.coastline:
            print('draw coast line')
            m.drawcoastlines()

        # Plot DEM
        if inps.dem_file:
            print('plotting DEM background ...')
            m = pp.plot_dem_lalo(m, dem, inps.geo_box, vars(inps))

        # Plot Data
        print('plotting Data ...')
        im = m.imshow(data, cmap=inps.colormap, origin='upper', vmin=inps.disp_min, vmax=inps.disp_max,\
                      alpha=inps.transparency, interpolation='nearest', animated=inps.animation)

        # Scale Bar
        if inps.disp_scalebar:
            print('plot scale bar')
            if not inps.scalebar:  inps.scalebar=[999,999,999]
            # Default Distance - 20% of data width
            if inps.scalebar[0] == 999.0:
                gc = pyproj.Geod(a=m.rmajor, b=m.rminor) 
                az12, az21, wid_dist = gc.inv(inps.geo_box[0], inps.geo_box[3], inps.geo_box[2], inps.geo_box[3])
                inps.scalebar[0] = ut.round_to_1(wid_dist*0.2)
            # Default center - Lower Left Corner
            if inps.scalebar[1] == 999.0:  inps.scalebar[1] = inps.geo_box[3]+0.1*(inps.geo_box[1]-inps.geo_box[3])
            if inps.scalebar[2] == 999.0:  inps.scalebar[2] = inps.geo_box[0]+0.2*(inps.geo_box[2]-inps.geo_box[0])
            # Draw scale bar
            m.drawscale(inps.scalebar[1], inps.scalebar[2], inps.scalebar[0], ax=ax,\
                        font_size=inps.font_size, color=inps.font_color)

        # Lat Lon labels
        if inps.lalo_label:
            print('plot lat/lon labels')
            m.draw_lalo_label(inps.geo_box, ax=ax, lalo_step=inps.lalo_step, font_size=inps.font_size, color=inps.font_color)
        else:
            ax.tick_params(labelsize=inps.font_size, colors=inps.font_color)
        
        # Plot Seed Point
        if inps.disp_seed and inps.seed_lalo:
            ax.plot(inps.seed_lalo[1], inps.seed_lalo[0], inps.seed_color+inps.seed_symbol, ms=inps.seed_size)
            print('plot reference point')

        # Status bar
        def format_coord(x,y):
            col = subset.coord_geo2radar(x, meta_dict, 'lon') - inps.pix_box[0]
            row = subset.coord_geo2radar(y, meta_dict, 'lat') - inps.pix_box[1]
            if 0 <= col < data.shape[1] and 0 <= row < data.shape[0]:
                z = data[row, col]
                if inps.dem_file:
                    dem_col = subset.coord_geo2radar(x, dem_meta_dict, 'lon') - inps.dem_pix_box[0]
                    dem_row = subset.coord_geo2radar(y, dem_meta_dict, 'lat') - inps.dem_pix_box[1]
                    h = dem[dem_row, dem_col]
                    return 'lon=%.4f, lat=%.4f, value=%.4f, elev=%.1f m, x=%.1f, y=%.1f'\
                           %(x,y,z,h,col+inps.pix_box[0],row+inps.pix_box[1])
                else:
                    return 'lon=%.4f, lat=%.4f, value=%.4f, x=%.1f, y=%.1f'\
                           %(x,y,z,col+inps.pix_box[0],row+inps.pix_box[0])
            else:
                return 'lon=%.4f, lat=%.4f'%(x,y)
        ax.format_coord = format_coord
     
    #-------------------- 2.2 Plot in Y/X-coordinate ------------------------------------------------#
    else:
        print('plotting in Y/X coordinate ...')
        
        # Plot DEM
        if inps.dem_file:
            print('plotting DEM background ...')
            ax = pp.plot_dem_yx(ax, dem, vars(inps))

        # Plot Data
        print('plotting Data ...')
        im = ax.imshow(data, cmap=inps.colormap, vmin=inps.disp_min, vmax=inps.disp_max,\
                       alpha=inps.transparency, interpolation='nearest')

        ax.tick_params(labelsize=inps.font_size)

        # Plot Seed Point
        if inps.disp_seed and inps.seed_yx:
            ax.plot(inps.seed_yx[1]-inps.pix_box[0], inps.seed_yx[0]-inps.pix_box[1],\
                    inps.seed_color+inps.seed_symbol, ms=inps.seed_size)
            print('plot reference point')

        ax.set_xlim(-0.5, np.shape(data)[1]-0.5)
        ax.set_ylim(np.shape(data)[0]-0.5, -0.5)

        # Status bar
        def format_coord(x,y):
            col = int(x)
            row = int(y)
            if 0 <= col < data.shape[1] and 0 <= row < data.shape[0]:
                z = data[row,col]
                if inps.dem_file:
                    h = dem[row,col]
                    return 'x=%.4f,  y=%.4f,  elev=%.1f m,  value=%.4f'%(x,y,h,z)
                else:
                    return 'x=%.4f,  y=%.4f,  value=%.4f'%(x,y,z)
            else:
                return 'x=%.4f,  y=%.4f'%(x,y)
        ax.format_coord = format_coord


    #-------------------- 3 Figure Setting --------------------------------------------------------#
    # 3.1 Colorbar
    if inps.disp_cbar:
        # Colorbar Extend
        if not inps.cbar_ext:
            if   inps.disp_min <= data_min and inps.disp_max >= data_max: inps.cbar_ext='neither'
            elif inps.disp_min >  data_min and inps.disp_max >= data_max: inps.cbar_ext='min'
            elif inps.disp_min <= data_min and inps.disp_max <  data_max: inps.cbar_ext='max'
            else:  inps.cbar_ext='both'
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", "3%", pad="3%")
        if inps.wrap and 'radian' in inps.disp_unit:
            cbar = plt.colorbar(im, cax=cax, ticks=[-np.pi, 0, np.pi])
            cbar.ax.set_yticklabels([r'-$\pi$', '0', r'$\pi$'])
        else:
            cbar = plt.colorbar(im, cax=cax, extend=inps.cbar_ext)
        if inps.cbar_nbins:
            cbar.locator = ticker.MaxNLocator(nbins=inps.cbar_nbins)
            cbar.update_ticks()
        cbar.ax.tick_params(labelsize=inps.font_size, colors=inps.font_color)
        if not inps.cbar_label:
            cbar.set_label(inps.disp_unit, fontsize=inps.font_size, color=inps.font_color)
        else:
            cbar.set_label(inps.cbar_label, fontsize=inps.font_size, color=inps.font_color)

    # 3.2 Title
    if inps.disp_title:
        ax.set_title(inps.fig_title, fontsize=inps.font_size, color=inps.font_color)

    # 3.3 Flip Left-Right / Up-Down
    if inps.flip_lr:
        print('flip figure left and right')
        ax.invert_xaxis()
    if inps.flip_ud:
        print('flip figure up and down')
        ax.invert_yaxis()

    # 3.4 Turn off axis
    if not inps.disp_axis:
        ax.axis('off')
        print('turn off axis display')

    # 3.5 Turn off tick label
    if not inps.disp_tick:
        #ax.set_xticklabels([])
        #ax.set_yticklabels([])
        print("NO TICKS")
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])

    # Figure Output
    if inps.save_fig:
        print('save figure to '+inps.outfile)
        plt.savefig(inps.outfile, bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)

    return ax, inps


##################################################################################################
EXAMPLE='''example:
  view.py SanAndreas.dem
  view.py velocity.h5 -u cm -m -2 -M 2 -c bwr --mask Mask_tempCoh.h5 -d SanAndreas.dem

  view.py timeseries.h5 
  view.py unwrapIfgram.h5 070927-100217
  view.py Wrapped.h5    -n 5
  view.py geomap_4rlks.trans range

  # Display in subset:
  view.py velocity.h5 -x 100 600     -y 200 800
  view.py velocity.h5 -l 31.05 31.10 -L 130.05 130.10

  # Exclude Dates:
  view.py timeseries.h5 -ex drop_date.txt

  # Reference:
  view.py velocity.h5   --ref-yx   210 566
  view.py timeseries.h5 --ref-date 20101120

  # Save and Output:
  view.py velocity.h5 --save
  view.py velocity.h5 -o velocity.pdf
  view.py velocity.h5 --nodisplay
'''

PLOT_TEMPLATE='''Plot Setting:
  plot.name          = 'Yunjun et al., 2016, AGU, Fig 4f'
  plot.type          = LOS_VELOCITY
  plot.startDate     = 
  plot.endDate       = 
  plot.displayUnit   = cm/yr
  plot.displayMin    = -2
  plot.displayMax    = 2
  plot.colormap      = jet
  plot.subset.lalo   = 33.05:33.15, 131.15:131.27
  plot.seed.lalo = 33.0651, 131.2076
'''


def cmdLineParse(argv):
    parser = argparse.ArgumentParser(description='Plot InSAR Product in 2D',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    ##### Input 
    infile = parser.add_argument_group('Input File', 'File/Dataset to display')
    infile.add_argument('file', help='file for display')
    infile.add_argument('epoch', nargs='*', help='optional - date/epoch(s) to display')
    infile.add_argument('-n','--epoch-num', dest='epoch_num', metavar='NUM', type=int, nargs='*', default=[],\
                        help='optional - order number of date/epoch(s) to display')
    infile.add_argument('--exclude','--ex', dest='exclude_epoch', metavar='EPOCH', nargs='*',\
                        help='dates will not be displayed')
    infile.add_argument('--mask', dest='mask_file', metavar='FILE',\
                        help='mask file for display')
    infile.add_argument('--zero-mask', dest='zero_mask', action='store_true', help='mask pixels with zero value.')

    ##### Output
    outfile = parser.add_argument_group('Output', 'Save figure and write to file(s)')
    outfile.add_argument('--save', dest='save_fig', action='store_true',\
                         help='save the figure')
    outfile.add_argument('--nodisplay', dest='disp_fig', action='store_false',\
                         help='save and do not display the figure')
    outfile.add_argument('-o','--outfile',\
                         help="save the figure with assigned filename.\n"
                              "By default, it's calculated based on the input file name.")

    ###### Data Display Option
    disp = parser.add_argument_group('Display Options', 'Options to adjust the dataset display')
    disp.add_argument('-m', dest='disp_min', type=float, help='minimum value of color scale')
    disp.add_argument('-M', dest='disp_max', type=float, help='maximum value of color scale')
    disp.add_argument('-u','--unit', dest='disp_unit', metavar='UNIT',\
                      help='unit for display.  Its priority > wrap')
    disp.add_argument('--scale', dest='disp_scale', metavar='NUM', type=float, default=1.0,\
                      help='display data in a scaled range. \n'
                           'Equivelant to data*input_scale')
    disp.add_argument('-c','--colormap', dest='colormap',\
                      help='colormap used for display, i.e. jet, RdBu, hsv, jet_r etc.\n'
                           'Support colormaps in Matplotlib - http://matplotlib.org/users/colormaps.html')
    disp.add_argument('--projection', dest='map_projection', default='cyl',\
                      help='map projection when plotting in geo-coordinate. \n'
                           'Reference - http://matplotlib.org/basemap/users/mapsetup.html\n\n')

    disp.add_argument('--wrap', action='store_true',\
                      help='re-wrap data to display data in fringes.')
    disp.add_argument('--opposite', action='store_true',\
                      help='display in opposite sign, equivalent to multiply data by -1.')
    disp.add_argument('--flip-lr', dest='flip_lr', action='store_true', help='flip left-right')
    disp.add_argument('--flip-ud', dest='flip_ud', action='store_true', help='flip up-down')
    disp.add_argument('--multilook-num', dest='multilook_num', type=int, default=1, \
                      help='multilook data in X and Y direction with a factor for display')
    disp.add_argument('--nomultilook', '--no-multilook', dest='multilook', action='store_false',\
                      help='do not multilook, for high quality display. \n'
                           'If multilook and multilook_num=1, multilook_num will be estimated automatically.\n'
                           'Useful when displaying big datasets.')
    disp.add_argument('--alpha', dest='transparency', type=float,\
                      help='Data transparency. \n'
                           '0.0 - fully transparent, 1.0 - no transparency.')
    disp.add_argument('--plot-setting', dest='disp_setting_file',\
                      help='Template file with plot setting.\n'+PLOT_TEMPLATE)

    ##### DEM
    dem = parser.add_argument_group('DEM','display topography in the background')
    dem.add_argument('-d','--dem', dest='dem_file', metavar='DEM_FILE',\
                     help='DEM file to show topography as background')
    dem.add_argument('--dem-noshade', dest='disp_dem_shade', action='store_false',\
                     help='do not show DEM shaded relief')
    dem.add_argument('--dem-nocontour', dest='disp_dem_contour', action='store_false',\
                     help='do not show DEM contour lines')
    dem.add_argument('--contour-smooth', dest='dem_contour_smooth', type=float, default=3.0,\
                     help='Background topography contour smooth factor - sigma of Gaussian filter. \n'
                          'Default is 3.0; set to 0.0 for no smoothing.')
    dem.add_argument('--contour-step', dest='dem_contour_step', metavar='NUM', type=float, default=200.0,\
                     help='Background topography contour step in meters. \n'
                          'Default is 200 meters.')

    ###### Subset
    subset = parser.add_argument_group('Subset','Display dataset in subset range')
    subset.add_argument('-x', dest='subset_x', type=int, nargs=2, metavar='X', \
                        help='subset display in x/cross-track/range direction')
    subset.add_argument('-y', dest='subset_y', type=int, nargs=2, metavar='Y', \
                        help='subset display in y/along-track/azimuth direction')
    subset.add_argument('-l','--lat', dest='subset_lat', type=float, nargs=2, metavar='LAT', \
                        help='subset display in latitude')
    subset.add_argument('-L','--lon', dest='subset_lon', type=float, nargs=2, metavar='LON', \
                        help='subset display in longitude')
    #subset.add_argument('--pixel-box', dest='pix_box', type=tuple,\
    #                    help='subset display in box define in pixel coord (x_start, y_start, x_end, y_end).\n'
    #                         'i.e. (100, 500, 1100, 2500)')
    #subset.add_argument('--geo-box', dest='geo_box', type=tuple,\
    #                    help='subset display in box define in geo coord (UL_lon, UL_lat, LR_lon, LR_lat).\n'
    #                         'i.e. (130.2, 33.8, 131.2, 31.8)')

    ##### Reference
    ref = parser.add_argument_group('Reference','Show / Modify reference in time and space for display')
    ref.add_argument('--ref-date', dest='ref_date', metavar='DATE', \
                     help='Change reference date for display')
    ref.add_argument('--ref-lalo', dest='seed_lalo', metavar=('LAT','LON'), type=float, nargs=2,\
                     help='Change referene point LAT LON for display')
    ref.add_argument('--ref-yx', dest='seed_yx', metavar=('Y','X'), type=int, nargs=2,\
                     help='Change referene point Y X for display')
    ref.add_argument('--noreference', dest='disp_seed', action='store_false', help='do not show reference point')
    ref.add_argument('--ref-color', dest='seed_color', metavar='COLOR', default='k',\
                     help='marker color of reference point')
    ref.add_argument('--ref-symbol', dest='seed_symbol', metavar='SYMBOL', default='s',\
                     help='marker symbol of reference point')
    ref.add_argument('--ref-size', dest='seed_size', metavar='SIZE_NUM', type=int, default=10,\
                     help='marker size of reference point, default: 10')

    ##### Vectors
    #vec = parser.add_argument_group('Vectors','Plot vector geometry')
    #vec.add_argument('--point-yx', dest='point_yx', type=int, nargs='')

    ##### Figure 
    fig = parser.add_argument_group('Figure','Figure settings for display')
    fig.add_argument('-s','--fontsize', dest='font_size', type=int, help='font size')
    fig.add_argument('--fontcolor', dest='font_color', default='k', help='font color')
    fig.add_argument('--dpi', dest='fig_dpi', metavar='DPI', type=int, default=150,\
                     help='DPI - dot per inch - for display/write')
    fig.add_argument('-r','--row', dest='fig_row_num', type=int, default=1, help='subplot number in row')
    fig.add_argument('-p','--col', dest='fig_col_num', type=int, default=1, help='subplot number in column')
    fig.add_argument('--noaxis', dest='disp_axis', action='store_false', help='do not display axis')
    fig.add_argument('--nocbar','--nocolorbar', dest='disp_cbar', action='store_false', help='do not display colorbar')
    fig.add_argument('--cbar-nbins', dest='cbar_nbins', type=int, help='number of bins for colorbar')
    fig.add_argument('--cbar-ext', dest='cbar_ext', default=None, choices={'neither','min','max','both',None},\
                     help='Extend setting of colorbar; based on data stat by default.')
    fig.add_argument('--cbar-label', dest='cbar_label', default=None, help='colorbar label')
    fig.add_argument('--notitle', dest='disp_title', action='store_false', help='do not display title')
    fig.add_argument('--notick', dest='disp_tick', action='store_false', help='do not display tick in x/y axis')
    fig.add_argument('--title-in', dest='fig_title_in', action='store_true', help='draw title in/out of axes')
    fig.add_argument('--figtitle', dest='fig_title', help='Title shown in the figure.')
    fig.add_argument('--figsize', dest='fig_size', metavar=('WID','LEN'), type=float, nargs=2,\
                      help='figure size in inches - width and length')
    fig.add_argument('--figext', dest='fig_ext',\
                     default='.png', choices=['.emf','.eps','.pdf','.png','.ps','.raw','.rgba','.svg','.svgz'],\
                     help='File extension for figure output file')
    fig.add_argument('--fignum', dest='fig_num', type=int, help='number of figure windows')
    fig.add_argument('--wspace', dest='fig_wid_space', type=float, default=0.05,\
                     help='width space between subplots in inches')
    fig.add_argument('--hspace', dest='fig_hei_space', type=float, default=0.05,\
                     help='height space between subplots in inches')
    fig.add_argument('--coord', dest='fig_coord', choices=['radar','geo'], default='geo',\
                     help='Display in radar/geo coordination system, for geocoded file only.')
    fig.add_argument('--animation', action='store_true', help='enable animation mode')
    
    ##### Map
    map_group = parser.add_argument_group('Map', 'Map settings for display')
    map_group.add_argument('--coastline', action='store_true', help='Draw coastline.')
    map_group.add_argument('--resolution', default='c', choices={'c','l','i','h','f',None}, \
                           help='Resolution of boundary database to use.\n'+\
                                'c (crude, default), l (low), i (intermediate), h (high), f (full) or None.')
    map_group.add_argument('--lalo-label', dest='lalo_label', action='store_true',\
                           help='Show N, S, E, W tick label for plot in geo-coordinate.\n'
                                'Useful for final figure output.')
    map_group.add_argument('--lalo-step', dest='lalo_step', type=float, help='Lat/lon step for lalo-label option.')
    map_group.add_argument('--scalebar', nargs=3, metavar=('DISTANCE','LAT_C','LON_C'), type=float,\
                           help='set scale bar with DISTANCE in meters centered at [LAT_C, LON_C]\n'+\
                                'set to 999 to use automatic value, e.g.\n'+\
                                '--scalebar 2000 33.06 131.18\n'+\
                                '--scalebar 500  999   999\n'+\
                                '--scalebar 999  33.06 131.18')
    map_group.add_argument('--noscalebar', dest='disp_scalebar', action='store_false', help='do not display scale bar.')

    inps = parser.parse_args(argv)
    # If output flie name assigned or figure shown is turned off, turn on the figure save
    if inps.outfile or not inps.disp_fig:
        inps.save_fig = True
    if inps.coastline and inps.resolution in ['c','l']:
        inps.resolution = 'i'
    if inps.lalo_step:
        inps.lalo_label = True
    return inps


#########################################  Main Function  ########################################
def main(argv):
    global fig

    inps = cmdLineParse(argv)
    if not inps.disp_fig:
        plt.switch_backend('Agg')
    print('\n******************** Display ********************')
    print(inps.disp_tick)
    # File Basic Info
    ext = os.path.splitext(inps.file)[1].lower()
    fbase = os.path.splitext(os.path.basename(inps.file))[0]

    try:
        atr = readfile.read_attribute(inps.file)
    except:
        sys.exit('Can not read file: '+inps.file)

    print('Input file is '+atr['PROCESSOR']+' '+atr['FILE_TYPE']+': '+inps.file)
    k = atr['FILE_TYPE']
    width = int(float(atr['WIDTH']))
    length = int(float(atr['LENGTH']))
    print('file size in y/x: '+atr['LENGTH']+', '+atr['WIDTH'])

    #------------------------------ Epoch/Date Info -------------------------------------------#
    # Read "epoch list to display' and 'reference date' for multi-dataset files
    if k in multi_group_hdf5_file+multi_dataset_hdf5_file+['HDFEOS', 'GIANT_TS']:
        # Read Epoch List
        h5file = h5py.File(inps.file,'r')
        if k in ['HDFEOS']:
            epochList = h5file.attrs['DATE_TIMESERIES'].split()
        elif k in ['GIANT_TS']:
            epochList = [dt.fromordinal(int(i)).strftime('%Y%m%d') for i in h5file['dates'][:].tolist()]
            print('dates:')
            print(epochList)
        else:
            epochList = sorted(h5file[k].keys())
        #h5file.close()

        # Epochs to display
        inps.epoch = get_epoch_full_list_from_input(epochList, inps.epoch, inps.epoch_num)[0]
        # If no epoch info input, display all the epochs
        if not inps.epoch:
            inps.epoch = list(epochList)

        # Epochs to exclude
        if inps.exclude_epoch:
            if os.path.isfile(inps.exclude_epoch[0]):
                inps.exclude_epoch = ptime.read_date_list(inps.exclude_epoch[0])
            inps.exclude_epoch = get_epoch_full_list_from_input(epochList, inps.exclude_epoch)[0]
            inps.epoch = sorted(list(set(inps.epoch) - set(inps.exclude_epoch)))
        else:
            inps.exclude_epoch = []
        epochNum = len(inps.epoch)

        # Output message
        if k in multi_dataset_hdf5_file:
            print('num of dates in file %s : %d' % (inps.file, len(epochList)))
            print('dates to exclude (%d):' % (len(inps.exclude_epoch)));     print(inps.exclude_epoch)
            print('dates to display (%d):' % (len(inps.epoch)));             print(inps.epoch)
        else:
            print('num of dates in file   :'+str(len(epochList)))
            print('num of dates to exclude:'+str(len(inps.exclude_epoch)))
            print('num of dates to display:'+str(epochNum))

        # Reference Date
        if inps.ref_date:
            ref_date_list = get_epoch_full_list_from_input(epochList, [inps.ref_date])[0]
            if ref_date_list:
                inps.ref_date = ref_date_list[0]
                #if inps.ref_date == atr['REF_DATE']:
                #    inps.ref_date = None
                #    print 'Input date is already the reference date: '+atr['REF_DATE']
            else:
                print('input reference date is not included in file '+inps.file)
                print('input reference date: '+inps.ref_date)
                print('dates in file '+inps.file+':\n'+str(epochList))
                inps.ref_date = None

        print('------------------------------------------------------------------------')
        # Raise Error if no epoch left
        if epochNum == 0:
            raise Exception('Zero epoch found!')

    # for single-dataset file
    elif k.lower() in ['.trans','.utm_to_rdc']:
        epochList = ['range','azimuth']
    elif fbase.startswith('los'):
        epochList = ['incidenceAngle','headingAngle']
    else:
        epochList = ['']

    if not inps.epoch:
        inps.epoch = epochList
    else:
        inps.epoch = get_epoch_full_list_from_input(epochList, inps.epoch)[0]

    epochNum = len(inps.epoch)


    #------------------------------ Update Plot Inps with metadata dict ----------------#
    if inps.disp_setting_file:
        inps = update_plot_inps_with_display_setting_file(inps, inps.disp_setting_file)
    inps = update_plot_inps_with_meta_dict(inps, atr)

    ## Exit if 1) save and do not display figure, a.k.a. generate figure file,
    ##     and 2) figure file exists and newer than data file
    #if (inps.save_fig\
    #    and not inps.disp_fig\
    #    and not ut.update_file(inps.outfile, inps.file, check_readable=False)):
    #    return inps.outfile

    # Read mask file if inputed
    if inps.mask_file:
        try:
            atrMsk = readfile.read_attribute(inps.mask_file)
            if atrMsk['LENGTH'] == atr['LENGTH'] and atrMsk['WIDTH'] == atr['WIDTH']:
                msk = readfile.read(inps.mask_file, inps.pix_box, epoch='mask')[0]
                print('mask data with: '+os.path.basename(inps.mask_file))
            else:
                print('WARNING: input file has different size from mask file: %s. Continue without mask' % (inps.mask_file))
                inps.mask_file = None
        except:
            print('Can not open mask file: '+inps.mask_file)
            inps.mask_file = None
    elif k in ['HDFEOS']:
        inps.mask_file = inps.file
        h5msk = h5py.File(inps.file, 'r')
        msk = h5msk[k]['GRIDS']['timeseries']['quality'].get('mask')[:]
        h5msk.close()
        print('mask %s data with contained mask dataset.' % (k))
    elif inps.file.endswith('PARAMS.h5'):
        inps.mask_file = inps.file
        h5msk = h5py.File(inps.file, 'r')
        msk = h5msk['cmask'][:] == 1.
        h5msk.close()
        print('mask data with contained cmask dataset')

    ############################### Read Data and Display ###############################
    ##### Display One Dataset
    print("EPOCH NUM IS"+str(epochNum))
    if epochNum == 1:
        # Read Data
        print('reading data ...')
        data, atr = readfile.read(inps.file, inps.pix_box, inps.epoch[0])
        
        if inps.zero_mask:
            data[data==0] = np.nan
        if inps.ref_date:
            ref_data = readfile.read(inps.file, inps.pix_box, inps.ref_date)[0]
            data -= ref_data
            del ref_data
        # Mask Data
        if inps.mask_file:
            data = mask.mask_matrix(data, msk)
            del msk

        # Figure Setting 
        if not inps.font_size:
            inps.font_size = 16
        if not inps.fig_size:
            plot_shape = [data.shape[1]*1.25, data.shape[0]]
            fig_scale = min(pysar.figsize_single_min/min(plot_shape), pysar.figsize_single_max/max(plot_shape))
            inps.fig_size = [np.rint(i*fig_scale*2)/2 for i in plot_shape]
        print('create figure in size: '+str(inps.fig_size))
        fig = plt.figure(figsize=inps.fig_size)
        ax = fig.add_axes([0.1,0.1,0.8,0.8])

        # Plotting
        ax, inps = plot_matrix(ax, data, atr, inps)
        del data

        if inps.disp_fig:
            print('showing ...')
            plt.show()


    ##### Display Multiple Multiple Datasets
    else:
        #plt.switch_backend('Agg')   #to depress the warning from tight_layout.
        ##### Figure Setting 
        if not inps.font_size:
            inps.font_size = 12
        if not inps.fig_size:
            # Get screen size in inch
            #screen_dpi = plt.figure().dpi
            #import Tkinter as tk
            #root = tk.Tk()
            #screen_width_res = root.winfo_screenwidth() / screen_dpi
            #screen_height_res = root.winfo_screenheight() / screen_dpi
            inps.fig_size = pysar.figsize_multi
        print('create figure in size: '+str(inps.fig_size))

        # Figure number (<= 200 subplots per figure)
        if not inps.fig_num:
            inps.fig_num = 1
            while epochNum/float(inps.fig_num) > 200.0:
                inps.fig_num += 1

        # Row/Column number
        if inps.fig_row_num==1 and inps.fig_col_num==1:
            # calculate row and col number based on input info
            data_shape = [float(inps.pix_box[3]-inps.pix_box[1])*1.1, float(inps.pix_box[2]-inps.pix_box[0])]
            if not inps.disp_min and not inps.disp_max:
                fig_size4plot = inps.fig_size
            else:
                fig_size4plot = [inps.fig_size[0]*0.95, inps.fig_size[1]]

            inps.fig_row_num, inps.fig_col_num = pp.auto_row_col_num(epochNum, data_shape, fig_size4plot, inps.fig_num)

        inps.fig_num = np.ceil(float(epochNum) / float(inps.fig_row_num*inps.fig_col_num)).astype(int)
        print('dataset number: '+str(epochNum))
        print('row     number: '+str(inps.fig_row_num))
        print('column  number: '+str(inps.fig_col_num))
        print('figure  number: '+str(inps.fig_num))

        # Update multilook parameters with new num and col number
        if inps.multilook and inps.multilook_num == 1:
            inps.multilook, inps.multilook_num = check_multilook_input(inps.pix_box, inps.fig_row_num, inps.fig_col_num)

        # Output File Name
        if inps.outfile:
            inps.fig_ext  = os.path.splitext(inps.outfile)[1].lower()
            inps.outfile_base = os.path.basename(inps.outfile).split(inps.fig_ext)[0]
        else:
            inps.outfile_base = os.path.splitext(inps.file)[0]

            if (inps.pix_box[2]-inps.pix_box[0])*(inps.pix_box[3]-inps.pix_box[1]) < width*length:
                inps.outfile_base += '_sub'
            if inps.wrap:
                inps.outfile_base += '_wrap'
            if not inps.disp_scale == 1.0:
                inps.outfile_base += '_scale'+str(inps.disp_scale)
            if inps.opposite:
                inps.outfile_base += '_oppo'
            if inps.ref_date:
                inps.outfile_base += '_ref'+inps.ref_date
            if inps.exclude_epoch:
                inps.outfile_base += '_ex'

        ##### Aux Data
        # Reference date for timeseries
        if k == 'timeseries':
            if inps.ref_date:
                print('consider input reference date: '+inps.ref_date)
                ref_data = readfile.read(inps.file, inps.pix_box, inps.ref_date)[0]

        ref_yx = None
        if k in ['interferograms'] and 'REF_Y' in atr.keys():
            ref_yx = [int(atr[i]) for i in ['REF_Y','REF_X']]
            print('consider reference pixel of interferograms in y/x: %d/%d' % (ref_yx[0], ref_yx[1]))

        # Read DEM
        if inps.dem_file:
            print('reading DEM: '+os.path.basename(inps.dem_file)+' ...')
            dem, dem_meta_dict = readfile.read(inps.dem_file, inps.pix_box, epoch='height')
            if inps.multilook:
                dem = multilook_matrix(dem, inps.multilook_num, inps.multilook_num)

            # Shaded Relief and Contour
            if inps.disp_dem_shade: 
                print('show shaded relief DEM')
                ls = LightSource(azdeg=315, altdeg=45)
                dem_hillshade = ls.shade(dem, vert_exag=1.0, cmap=plt.cm.gray, vmin=-5000, vmax=np.nanmax(dem)+2000)
            if inps.disp_dem_contour:
                print('show contour: step = '+str(inps.dem_contour_step)+' m')
                dem_contour = ndimage.gaussian_filter(dem, sigma=inps.dem_contour_smooth, order=0)
                contour_sequence = np.arange(-6000, 9000, inps.dem_contour_step)

        ################## Plot Loop ####################

        ## Find min and value for all data, reference for better min/max setting next time
        all_data_min=0
        all_data_max=0

        #h5file = h5py.File(inps.file, 'r')
        # Check dropped interferograms
        drop_epoch_list = []
        if k in multi_group_hdf5_file and inps.disp_title:
            drop_epoch_list = sorted(list(set(inps.epoch) - \
                                          set(ut.check_drop_ifgram(h5file, print_msg=False))))
            print("mark interferograms with 'DROP_IFGRAM'='yes' in red colored title")

        ##### Loop 1 - Figures
        for j in range(1, inps.fig_num+1):
            # Output file name for current figure
            if inps.fig_num > 1:
                inps.outfile = inps.outfile_base+'_'+str(j)+inps.fig_ext
            else:
                inps.outfile = inps.outfile_base+inps.fig_ext
            fig_title = 'Figure '+str(j)+' - '+inps.outfile
            print('----------------------------------------')
            print(fig_title)
            # Open a new figure object
            fig = plt.figure(j, figsize=inps.fig_size)
            fig.canvas.set_window_title(fig_title)

            fig_data_min=0
            fig_data_max=0
            i_start = (j-1)*inps.fig_row_num*inps.fig_col_num
            i_end   = min([epochNum, i_start+inps.fig_row_num*inps.fig_col_num])
            ##### Loop 2 - Subplots
            prog_bar = ptime.progress_bar(maxValue=i_end-i_start, prefix='loading: ')
            for i in range(i_start, i_end):
                epoch = inps.epoch[i]
                ax = fig.add_subplot(inps.fig_row_num, inps.fig_col_num, i-i_start+1)
                prog_bar.update(i-i_start+1, suffix=str(i+1))

                # Read Data
                if k in multi_dataset_hdf5_file:
                    dset = h5file[k].get(epoch)
                    data = dset[inps.pix_box[1]:inps.pix_box[3], inps.pix_box[0]:inps.pix_box[2]]
                    if inps.ref_date:
                        data -= ref_data
                    try:    subplot_title = dt.strptime(epoch, '%Y%m%d').isoformat()[0:10]
                    except: subplot_title = epoch
                elif k in multi_group_hdf5_file:
                    if inps.fig_row_num*inps.fig_col_num > 100:
                        subplot_title = str(epochList.index(epoch)+1)
                    else:
                        subplot_title = str(epochList.index(epoch)+1)+'\n'+h5file[k][epoch].attrs['DATE12']
                    dset = h5file[k][epoch].get(epoch)
                    data = dset[inps.pix_box[1]:inps.pix_box[3], inps.pix_box[0]:inps.pix_box[2]]
                    if inps.zero_mask:
                        data[data==0] = np.nan
                    elif ref_yx:
                        data -= data[ref_yx[0], ref_yx[1]]
                elif k in ['HDFEOS']:
                    dset = h5file[k]['GRIDS']['timeseries']['observation'].get('displacement')[i,:,:]
                    data = dset[inps.pix_box[1]:inps.pix_box[3], inps.pix_box[0]:inps.pix_box[2]]
                    subplot_title = str(epoch)
                else:
                    data = readfile.read(inps.file, epoch=epoch)[0]
                    subplot_title = str(epoch)
                # mask
                if inps.mask_file:
                    data = mask.mask_matrix(data, msk)

                # Update data with plot inps
                data, inps = update_matrix_with_plot_inps(data, atr, inps)

                # Data Min/Max
                fig_data_min = np.nanmin([fig_data_min, np.nanmin(data)])
                fig_data_max = np.nanmax([fig_data_max, np.nanmax(data)])

                # Plot DEM
                if inps.dem_file and inps.disp_dem_shade:
                    ax.imshow(dem_hillshade, cmap='gray', interpolation='spline16')
                if inps.dem_file and inps.disp_dem_contour:
                    ax.contour(dem_contour, contour_sequence, origin='lower',colors='black',alpha=0.5)

                # Plot Data
                try:
                    im = ax.imshow(data, cmap=inps.colormap, interpolation='nearest', alpha=inps.transparency,\
                                   vmin=inps.disp_min, vmax=inps.disp_max)
                except:
                    im = ax.imshow(data, cmap=inps.colormap, interpolation='nearest', alpha=inps.transparency)

                ###### Subplot Setting
                # Tick and Label
                ax.set_yticklabels([])
                ax.set_xticklabels([])
                ax.set_xticks([])
                ax.set_yticks([])
                # Title
                if inps.disp_title:
                    if not inps.fig_title_in:
                        if epoch in drop_epoch_list:
                            ax.set_title(subplot_title, fontsize=inps.font_size, color='crimson', fontweight='bold')
                        else:
                            ax.set_title(subplot_title, fontsize=inps.font_size)
                    else:
                        pp.add_inner_title(ax, subplot_title, loc=1)   
                # Flip Left-Right / Up-Down
                if inps.flip_lr:        ax.invert_xaxis()
                if inps.flip_ud:        ax.invert_yaxis()
                # Turn off axis
                if not inps.disp_axis:
                    ax.axis('off')

            ##### Figure Setting - End of Loop 2
            prog_bar.close()
            fig.tight_layout()
            # Min and Max for this figure
            all_data_min = np.nanmin([all_data_min, fig_data_min])
            all_data_max = np.nanmax([all_data_max, fig_data_max])
            print('data    range: [%.2f, %.2f] %s' % (fig_data_min, fig_data_max, inps.disp_unit))
            try:  print('display range: [%.2f, %.2f] %s' % (disp_min, disp_max, inps.disp_unit))
            except: pass

            # Colorbar
            if not inps.disp_min and not inps.disp_max:
                print('Note: different color scale for EACH subplot!')
            else:
                if not inps.cbar_ext:
                    if   inps.disp_min <= fig_data_min and inps.disp_max >= fig_data_max: inps.cbar_ext='neither'
                    elif inps.disp_min >  fig_data_min and inps.disp_max >= fig_data_max: inps.cbar_ext='min'
                    elif inps.disp_min <= fig_data_min and inps.disp_max <  fig_data_max: inps.cbar_ext='max'
                    else:  inps.cbar_ext='both'
                print('show colorbar')
                #fig.subplots_adjust(wspace=inps.fig_wid_space, hspace=inps.fig_hei_space, right=0.965)
                fig.subplots_adjust(right=0.95)
                cax = fig.add_axes([0.96, 0.25, 0.01, 0.5])
                if inps.wrap and 'radian' in inps.disp_unit:
                    cbar = plt.colorbar(im, cax=cax, ticks=[-np.pi, 0, np.pi])
                    cbar.ax.set_yticklabels([r'-$\pi$', '0', r'$\pi$'])
                else:
                    cbar = plt.colorbar(im, cax=cax, extend=inps.cbar_ext)
                if inps.cbar_nbins:
                    cbar.locator = ticker.MaxNLocator(nbins=inps.cbar_nbins)
                    cbar.update_ticks()
                cbar.ax.tick_params(labelsize=inps.font_size, colors=inps.font_color)
                if not inps.cbar_label:
                    cbar.set_label(inps.disp_unit, fontsize=inps.font_size, color=inps.font_color)
                else:
                    cbar.set_label(inps.cbar_label, fontsize=inps.font_size, color=inps.font_color)

            # Save Figure
            if inps.save_fig:
                print('save figure to '+inps.outfile)
                fig.savefig(inps.outfile, bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)
                if not inps.disp_fig:
                    fig.clf()

        ##### End of Loop 1
        try: h5file.close()
        except: pass
        print('----------------------------------------')
        print('all data range: [%f, %f] %s' % (all_data_min, all_data_max, inps.disp_unit))
        if inps.disp_min and inps.disp_max:
            print('display  range: [%f, %f] %s' % (inps.disp_min, inps.disp_max, inps.disp_unit))

        # Display Figure
        if inps.disp_fig:
            print('showing ...')
            plt.show()


##################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])


