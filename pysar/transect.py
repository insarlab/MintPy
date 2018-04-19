#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os, sys, glob
import argparse
import h5py
import numpy as np
import scipy.ndimage
import scipy.io as sio
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from pysar.utils import readfile, utils as ut


#####################################################################
def get_scale_from_disp_unit(disp_unit, data_unit):
    # Initial 
    scale = 1.0
    data_unit = data_unit.lower().split('/')
    disp_unit = disp_unit.lower().split('/')

    # if data and display unit is the same
    if disp_unit == data_unit:
        return scale, disp_unit

    # Calculate scaling factor  - 1
    # phase unit - length
    if data_unit[0] == 'm':
        if   disp_unit[0] == 'mm': scale *= 1000.0
        elif disp_unit[0] == 'cm': scale *= 100.0
        elif disp_unit[0] == 'dm': scale *= 10.0
        elif disp_unit[0] == 'km': scale *= 1/1000.0
        else:
            print('Unrecognized display unit: '+disp_unit[0])
            return
    else:
        print('Unrecognized data unit: '+data_unit[0])
        return

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

    return scale, disp_unit


#####################################################################
def read_lonlat_file(lonlat_file):
    '''Read Start/End lat/lon from lonlat text file in gmt format.
    Inputs:
        lonlat_file : text file in gmt lonlat point file
    Outputs:
        start/end_lalo : list of 2 float
    '''
    fll = open(lonlat_file,'r')
    lines = fll.read().splitlines()
    [lon0,lat0] = [float(i) for i in lines[1].split()]
    [lon1,lat1] = [float(i) for i in lines[2].split()]
    fll.close()
    
    start_lalo = [lat0, lon0]
    end_lalo = [lat1, lon1]
    return start_lalo, end_lalo


#####################################################################
def manual_select_start_end_point(File):
    '''Manual Select Start/End Point in display figure.'''
    print('reading '+File+' ...')
    data, atr = readfile.read(File)
    print('displaying '+File+' ...')
    fig = plt.figure()
    ax=fig.add_subplot(111)
    ax.imshow(data)
    
    xc=[]
    yc=[]
    print('please click on start and end point of the desired profile')
    print('then close the figure to continue')
    def onclick(event):
        if event.button==1:
            xcc, ycc = int(event.xdata), int(event.ydata)
            xc.append(xcc);  yc.append(ycc)
            print('x = '+str(xcc)+'\ny = '+str(ycc))
            ax.plot(xcc,ycc,'ro')
    cid = fig.canvas.mpl_connect('button_release_event', onclick)
    plt.show();

    start_yx = [yc[0], xc[0]]
    end_yx = [yc[1], xc[1]]
    return start_yx, end_yx


#####################################################################
def transect_yx(z,atr,start_yx,end_yx,interpolation='nearest'):
    '''Extract 2D matrix (z) value along the line [x0,y0;x1,y1]
    Ref link: http://stackoverflow.com/questions/7878398/how-to-e
              xtract-an-arbitrary-line-of-values-from-a-numpy-array
    
    Inputs:
        z        - (np.array)   2D data matrix
        atr      - (dictionary) 2D data matrix attribute dictionary
        start_yx - (list) y,x coordinate of start point
        end_yx   - (list) y,x coordinate of end   point
        interpolation - sampling/interpolation method, including:
                'nearest'  - nearest neighbour, by default
                'cubic'    - cubic interpolation
                'bilinear' - bilinear interpolation
    
    Output:
        transect - N*2 matrix containing distance - 1st col - and its corresponding 
                   values - 2nd col - along the line, N is the number of points.
    
    Example:
        transect = transect_yx(dem,demRsc,[10,15],[100,115])
    '''

    ## Extract the line
    [y0,x0] = start_yx
    [y1,x1] = end_yx
    length = int(np.hypot(x1-x0, y1-y0))
    x, y = np.linspace(x0, x1, length), np.linspace(y0, y1, length)
    transect = np.zeros([length,2])

    ## Y - Extract the value along the line
    if   interpolation.lower() == 'nearest' : zi = z[np.rint(y).astype(np.int), np.rint(x).astype(np.int)]
    elif interpolation.lower() == 'cubic'   : zi = scipy.ndimage.map_coordinates(z,np.vstack((y,x)))
    elif interpolation.lower() == 'bilinear': zi = scipy.ndimage.map_coordinates(z,np.vstack((y,x)),order=2)
    else:
        print('Unrecognized interpolation method: '+interpolation)
        print('Continue with nearest ...')
        zi = z[np.rint(y).astype(np.int), np.rint(x).astype(np.int)]  # nearest neighbour
    transect[:,1] = zi

    ## X - Distance along the line
    earth_radius = 6371.0e3;    # in meter
    try:
        atr['X_FIRST']
        [lat0,lat1] = ut.coord_radar2geo([y0,y1],atr,'y')
        x_step = float(atr['X_STEP'])*np.pi/180.0*earth_radius*np.cos((lat0+lat1)/2*np.pi/180)
        y_step = float(atr['Y_STEP'])*np.pi/180.0*earth_radius
    except:
        x_step = ut.range_ground_resolution(atr)
        y_step = ut.azimuth_ground_resolution(atr)
    dis_x = (x-x0)*x_step
    dis_y = (y-y0)*y_step
    transect[:,0] = np.hypot(dis_x,dis_y)

    return transect

def transect_lalo(z,atr,start_lalo,end_lalo,interpolation='nearest'):
    '''Extract 2D matrix (z) value along the line [start_lalo, end_lalo]'''
    [y0,y1] = ut.coord_geo2radar([start_lalo[0],end_lalo[0]],atr,'lat')
    [x0,x1] = ut.coord_geo2radar([start_lalo[1],end_lalo[1]],atr,'lon')
    transect = transect_yx(z,atr,[y0,x0],[y1,x1],interpolation)
    return transect

def transect_list(fileList, inps):
    '''Get transection along input line from file list
    Inputs:
        fileList : list of str, path of files to get transect
        inps     : Namespace including the following items:
                   start/end_lalo
                   start/end_yx
                   interpolation
    Outputs:
        transectList : list of N*2 matrix containing distance and its value
        atrList      : list of attribute dictionary, for each input file
    '''

    transectList = []
    atrList      = []
    for File in fileList:
        print('reading '+File)
        data, atr = readfile.read(File)
        if inps.start_lalo and inps.end_lalo:
            transect = transect_lalo(data, atr, inps.start_lalo, inps.end_lalo, inps.interpolation)
        else:
            transect = transect_yx(data, atr, inps.start_yx, inps.end_yx, inps.interpolation)
        transectList.append(transect)
        atrList.append(atr)
    return transectList, atrList


#####################################################################
EXAMPLE='''example:
  transect.py velocity.h5 -s 5290 5579 -e 12177 482
  transect.py velocity.h5 --start-lalo 30.125 129.988 --end-lalo 30.250 130.116
  transect.py velocity.h5 --line-file  transect_lonlat.xy -d gsi10m.dem
  transect.py AlosA*/velocity.h5 AlosD*/velocity.h5 --line-file  transect_lonlat.xy -d gsi10m.dem
'''

def create_parser():
    parser = argparse.ArgumentParser(description='Generate transect/profile along a line',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs='+', help='input file to show transection')
    parser.add_argument('-m','--min', dest='disp_min', type=float, help='minimum value for data display')
    parser.add_argument('-M','--max', dest='disp_max', type=float, help='maximum value for data display')
    parser.add_argument('-u','--unit', dest='disp_unit', default='cm', help='unit for data display. Default: cm')
    parser.add_argument('--offset', dest='disp_offset', type=float, default=3.0, help='offset between each data profile')
    parser.add_argument('--interpolation', default='nearest', choices=['nearest','bilinear','cubic'],\
                        help='interpolation method while extacting profile along the line')

    # Start / End Point
    end = parser.add_argument_group('Start and End Point of Profile')
    end.add_argument('-s','--start-yx', dest='start_yx', type=int, nargs=2,\
                     help='start point of the profile in pixel number [y, x]')
    end.add_argument('-e','--end-yx', dest='end_yx', type=int, nargs=2,\
                     help='end   point of the profile in pixel number [y, x]')
    end.add_argument('--start-lalo', dest='start_lalo', type=float, nargs=2,\
                     help='start point of the profile in [lat, lon]')
    end.add_argument('--end-lalo', dest='end_lalo', type=float, nargs=2,\
                     help='end   point of the profile in [lat, lon]')
    end.add_argument('--line-file', dest='lola_file',\
                     help='file with start and end point info in lon lat, same as GMT format.\n'
                          'i.e. transect_lonlat.xy:\n'
                          '>\n'
                          '131.1663    33.1157\n'
                          '131.2621    33.0860')

    # DEM
    dem = parser.add_argument_group('DEM','display topography in the bottom')
    dem.add_argument('-d','--dem', help='DEM file')
    dem.add_argument('--dem-min', dest='dem_disp_min', type=float, help='min display value for DEM display, in km')
    dem.add_argument('--dem-max', dest='dem_disp_max', type=float, help='max display value for DEM display, in km')

    # Output
    outfile = parser.add_argument_group('Output', 'Save figure and write to file(s)')
    outfile.add_argument('--save', dest='save_fig', action='store_true',\
                         help='save the figure')
    outfile.add_argument('--nodisplay', dest='disp_fig', action='store_false',\
                         help='save and do not display the figure')
    outfile.add_argument('-o','--outfile',\
                                help="save the figure with assigned filename.\n"
                                     "By default, it's calculated based on inputs.")

    # Figure 
    fig = parser.add_argument_group('Figure','Figure settings for display')
    fig.add_argument('--dpi', dest='fig_dpi', type=int, default=300, help='DPI - dot per inch - for display/write')
    fig.add_argument('--figsize', dest='fig_size', type=float, nargs=2, default=[7.0, 6.0],\
                     help='figure size in inches - width and length')
    fig.add_argument('--figext', dest='outfile_ext',\
                     default='.png', choices=['.emf','.eps','.pdf','.png','.ps','.raw','.rgba','.svg','.svgz'],\
                     help='File extension for figure output file')
    fig.add_argument('--fontsize', dest='font_size', type=int, help='font size')
    fig.add_argument('--ms','--markersize', dest='marker_size', type=float, default=2.0,\
                     help='Point marker size. Default: 2.0')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    inps.file = ut.get_file_list(inps.file)
    if inps.outfile or not inps.disp_fig:
        inps.save_fig = True
    return inps


############################ Main ###################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    print('\n**************** Transect *********************')
    print('number of file: '+str(len(inps.file)))
    print(inps.file)
    

    ##### Start / End Point Input
    # 1. lonlat file
    if inps.lola_file:
        inps.start_lalo, inps.end_lalo = read_lonlat_file(inps.lola_file)
    # 2. Manually select in window display
    if (not inps.start_yx or not inps.end_yx) and (not inps.start_lalo or not inps.end_lalo):
        print('No input yx/lalo found.')
        print('Continue with manually select start and end point.')
        inps.start_yx, inps.end_yx = manual_select_start_end_point(inps.file[0])
    # Message
    if inps.start_lalo and inps.end_lalo:
        print('Start point:  lat = '+str(inps.start_lalo[0])+', lon = '+str(inps.start_lalo[1]))
        print('End   point:  lat = '+str(inps.end_lalo[0])+', lon = '+str(inps.end_lalo[1]))
    else:
        print('Start point:  y = '+str(inps.start_yx[0])+', x = '+str(inps.start_yx[1]))
        print('End   point:  y = '+str(inps.end_yx[0])+', x = '+str(inps.end_yx[1]))

    ##### Get Transection/Profiles Data
    print('extract transect from input files ...')
    transectList, atrList = transect_list(inps.file, inps)
    if inps.dem:
        demTransectList, demAtrList = transect_list([inps.dem], inps)
    print('profile length: '+str(transectList[0][-1,0]/1000.0)+' km')

    ##### Plot
    # Figure 1 - Profile line in the 1st input file
    print('plot profile line in the 1st input file')
    fig0 = plt.figure()
    ax0 = fig0.add_axes([0.1,0.1,0.8,0.8])
    data0, atr0 = readfile.read(inps.file[0])
    ax0.imshow(data0)
    if inps.start_lalo and inps.end_lalo:
        [y0,y1] = ut.coord_geo2radar([inps.start_lalo[0], inps.end_lalo[0]], atr0, 'lat')
        [x0,x1] = ut.coord_geo2radar([inps.start_lalo[1], inps.end_lalo[1]], atr0, 'lon')
    else:
        [y0,y1] = [inps.start_yx[0], inps.end_yx[0]]
        [x0,x1] = [inps.start_yx[1], inps.end_yx[1]]
    ax0.plot([x0, x1], [y0, y1], 'ro-')
    ax0.set_xlim(0,np.shape(data0)[1])
    ax0.set_ylim(np.shape(data0)[0],0)
    ax0.set_title('Transect Line in '+inps.file[0])

    # Status bar
    def format_coord(x,y):
        col = int(x)
        row = int(y)
        if 0 <= col < data0.shape[1] and 0 <= row < data0.shape[0]:
            z = data0[row,col]
            if 'X_FIRST' in atr0.keys():
                lat = ut.coord_radar2geo(row, atr0, 'row')
                lon = ut.coord_radar2geo(col, atr0, 'col')
                return 'lon=%.4f, lat=%.4f, x=%.0f,  y=%.0f,  value=%.4f' % (lon, lat, x,y,z)
            else:
                return 'x=%.0f,  y=%.0f,  value=%.4f'%(x,y,z)
        else:
            return 'x=%.0f,  y=%.0f'%(x,y)
    ax0.format_coord = format_coord


    # Figure 2 - Transections/Profiles
    print('plot profiles')
    fig,ax = plt.subplots(figsize = inps.fig_size)
    # Plot 2.1 - Input Files
    if not inps.disp_unit:
        inps.disp_unit = atrList[0]['UNIT']
    inps.disp_scale, inps.disp_unit = get_scale_from_disp_unit(inps.disp_unit, atrList[0]['UNIT'])

    value_min = 0
    value_max = 0
    for i in range(len(inps.file)):
        # Profile Color based on Asc/Desc direction
        if atrList[i]['ORBIT_DIRECTION'][0].upper() == 'A':
            p_color = 'crimson'
        else:
            p_color = 'royalblue'
        # Plot
        distance = transectList[i][:,0]/1000.0          # km
        value = transectList[i][:,1]*inps.disp_scale - inps.disp_offset*i
        ax.plot(distance, value, '.', color=p_color, markersize=inps.marker_size)
        # Y Stat
        value_min = np.nanmin([value_min, np.nanmin(value)])
        value_max = np.nanmax([value_max, np.nanmax(value)])

    # Y axis
    if not inps.disp_min:
        inps.disp_min = np.floor(value_min - (value_max-value_min)*1.2/len(inps.file))
    if not inps.disp_max:
        inps.disp_max = np.ceil(value_max)
    ax.set_ylim(inps.disp_min, inps.disp_max)
    ax.set_ylabel('Mean LOS Velocity ('+inps.disp_unit+')', fontsize=inps.font_size)
    # X axis
    ax.set_xlabel('Distance (km)', fontsize=inps.font_size)
    ax.tick_params(which='both', direction='out', labelsize=inps.font_size)

    # Plot 2.2 - DEM
    if inps.dem:
        ax2 = ax.twinx()
        distance = demTransectList[0][:,0]/1000.0    # km
        value    = demTransectList[0][:,1]/1000.0    # km
        ax2.fill_between(distance, 0, value, facecolor='gray')

        # Y axis - display DEM in the bottom
        value_min = np.nanmin(value)
        value_max = np.nanmax(value)
        if not inps.dem_disp_min:
            inps.dem_disp_min = np.floor(value_min*2.0)/2.0
        if not inps.dem_disp_max:
            inps.dem_disp_max = np.ceil((value_max + (value_max-value_min)*(len(inps.file)+0.0))*2.0)/2.0
        ax2.set_ylim(inps.dem_disp_min, inps.dem_disp_max)
        ## Show lower part of yaxis
        #dem_tick = ax2.yaxis.get_majorticklocs()
        #dem_tick = dem_tick[:len(dem_tick)/2]
        #ax2.set_yticks(dem_tick)
        ax2.set_ylabel('Elevation (km)', fontsize=inps.font_size)
        ax2.tick_params(which='both', direction='out', labelsize=inps.font_size)

    ## X axis - Shared
    distanceMax = np.nanmax(transectList[0][:,0]/1000.0)   # in km
    plt.xlim(0, distanceMax)
    plt.tight_layout()

    ##### Output
    if not inps.outfile:
        figBase = 'transect_x'+str(x0)+'y'+str(y0)+'_x'+str(x1)+'y'+str(y1)
    else:
        figBase, inps.outfile_ext = os.path.splitext(inps.outfile)
        if not inps.outfile_ext:
            inps.outfile_ext = '.png'
    if inps.save_fig:
        print('writing >>> '+figBase+inps.outfile_ext)
        fig0.savefig(inps.file[-1]+inps.outfile_ext, bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)
        fig.savefig(figBase+inps.outfile_ext, bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)

        print('writing >>> '+figBase+'.mat')
        transect_mat = {}
        for i in range(len(inps.file)):
            project = atrList[i]['PROJECT_NAME']
            transect_mat[project] = transectList[i]
        if inps.dem:
            transect_mat['elevation'] = demTransectList[0]
        sio.savemat(figBase+'.mat', {'transection':transect_mat})

    # Display
    if inps.disp_fig:
        print('showing ...')
        plt.show()
    return


#####################################################################
if __name__ == '__main__':
    main()

