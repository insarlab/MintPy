#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, Zhang Yunjun, 2013               #
############################################################


import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
import scipy.ndimage

from mintpy.utils import plot as pp, readfile, utils as ut

#####################################################################
EXAMPLE = """example:
  transect.py velocity.h5 -s 5290 5579 -e 12177 482
  transect.py velocity.h5 --start-lalo 30.125 129.988 --end-lalo 30.250 130.116
  transect.py velocity.h5 --line-file  transect_lonlat.xy --dem gsi10m.dem

  # profile from multiple files
  transect.py AlosA*/velocity.h5 AlosD*/velocity.h5 --line-file  transect_lonlat.xy --dem gsi10m.dem
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Generate transect/profile along a line',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs='+',
                        help='input file to show transection')
    parser.add_argument('--dset', dest='dset', help='Dataset name to read')
    parser.add_argument('-m', '--min', dest='disp_min',
                        type=float, help='minimum value for data display')
    parser.add_argument('-M', '--max', dest='disp_max',
                        type=float, help='maximum value for data display')
    parser.add_argument('-u', '--unit', dest='disp_unit',
                        default='cm', help='unit for data display. Default: cm')
    parser.add_argument('--offset', dest='disp_offset', type=float,
                        default=3.0, help='offset between each data profile')
    parser.add_argument('--interpolation', default='nearest', choices=['nearest', 'bilinear', 'cubic'],
                        help='interpolation method while extacting profile along the line')

    # Start / End Point
    end = parser.add_argument_group('Start and End Point of Profile')
    end.add_argument('-s', '--start-yx', dest='start_yx', metavar=('Y0', 'X0'), type=int, nargs=2,
                     help='start point of the profile in pixel number [y, x]')
    end.add_argument('-e', '--end-yx', dest='end_yx', metavar=('Y1', 'X1'), type=int, nargs=2,
                     help='end   point of the profile in pixel number [y, x]')
    end.add_argument('--start-lalo', dest='start_lalo', metavar=('LAT0', 'LON0'), type=float, nargs=2,
                     help='start point of the profile in [lat, lon]')
    end.add_argument('--end-lalo', dest='end_lalo', metavar=('LAT1', 'LON1'), type=float, nargs=2,
                     help='end   point of the profile in [lat, lon]')
    end.add_argument('--line-file', dest='lola_file',
                     help='file with start and end point info in lon lat, same as GMT format.\n'
                          'i.e. transect_lonlat.xy:\n'
                          '>\n'
                          '131.1663    33.1157\n'
                          '131.2621    33.0860')

    # DEM
    dem = parser.add_argument_group('DEM', 'display topography in the bottom')
    dem.add_argument('--dem', help='DEM file')
    dem.add_argument('--dem-min', dest='dem_disp_min', type=float,
                     help='min display value for DEM display, in km')
    dem.add_argument('--dem-max', dest='dem_disp_max', type=float,
                     help='max display value for DEM display, in km')

    # Output
    outfile = parser.add_argument_group('Output', 'Save figure and write to file(s)')
    outfile.add_argument('--save', dest='save_fig', action='store_true',
                         help='save the figure/data')
    outfile.add_argument('--nodisplay', dest='disp_fig', action='store_false',
                         help='save and do not display the figure')
    outfile.add_argument('-o', '--outfile',
                         help="save the figure with assigned filename.\n"
                         "By default, it's calculated based on inputs.")

    # Figure
    fig = parser.add_argument_group('Figure', 'Figure settings for display')
    fig.add_argument('--dpi', dest='fig_dpi', type=int, default=300,
                     help='DPI - dot per inch - for display/write')
    fig.add_argument('--figsize', dest='fig_size', type=float, nargs=2, default=[7.0, 6.0],
                     help='figure size in inches - width and length')
    fig.add_argument('--figext', dest='outfile_ext',
                     default='.png', choices=['.emf', '.eps', '.pdf', '.png', '.ps',
                                              '.raw', '.rgba', '.svg', '.svgz'],
                     help='File extension for figure output file')
    fig.add_argument('--fontsize', dest='font_size',
                     type=int, help='font size')
    fig.add_argument('--ms', '--markersize', dest='marker_size', type=float, default=2.0,
                     help='Point marker size. Default: 2.0')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    inps.file = ut.get_file_list(inps.file)
    if inps.outfile or not inps.disp_fig:
        inps.save_fig = True
    return inps


#####################################################################
def read_lonlat_file(lonlat_file):
    """Read Start/End lat/lon from lonlat text file in gmt format.
    Inputs:
        lonlat_file : text file in gmt lonlat point file
    Outputs:
        start/end_lalo : list of 2 float
    """
    fll = open(lonlat_file)
    lines = fll.read().splitlines()
    [lon0, lat0] = [float(i) for i in lines[1].split()]
    [lon1, lat1] = [float(i) for i in lines[2].split()]
    fll.close()

    start_lalo = [lat0, lon0]
    end_lalo = [lat1, lon1]
    return start_lalo, end_lalo


#####################################################################
def manual_select_start_end_point(File, dset=None):
    """Manual Select Start/End Point in display figure."""
    print('reading '+File+' ...')
    data, atr = readfile.read(File, datasetName=dset)
    print('displaying '+File+' ...')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(data)

    lines = []
    def draw_line():
        xy = plt.ginput(2)
        x = [p[0] for p in xy]
        y = [p[1] for p in xy]
        line = plt.plot(x,y)
        ax.figure.canvas.draw()
        return line

    line = draw_line()
    #xc = []
    #yc = []
    #print('1) Click on start and end point of the desired profile')
    #print('2) Close the figure to continue the profile plotting')
#
    #def onclick(event):
    #    if event.button == 1:
    #        xcc, ycc = int(event.xdata), int(event.ydata)
    #        xc.append(xcc)
    #        yc.append(ycc)
    #        print('({}, {})'.format(xcc, ycc))
    #        ax.plot(xcc, ycc, 'ro')
    #cid = fig.canvas.mpl_connect('button_release_event', onclick)
    plt.show()

    #start_yx = [yc[0], xc[0]]
    #end_yx = [yc[1], xc[1]]
    #return start_yx, end_yx
    return [0,0], [10,10]


#####################################################################
def transect_yx(z, atr, start_yx, end_yx, interpolation='nearest'):
    """Extract 2D matrix (z) value along the line [x0,y0;x1,y1]
    Ref link: http://stackoverflow.com/questions/7878398/how-to-e
              xtract-an-arbitrary-line-of-values-from-a-numpy-array

    Parameters: z : (np.array)   2D data matrix
                atr : (dict) 2D data matrix attribute dictionary
                start_yx : (list) y,x coordinate of start point
                end_yx : (list) y,x coordinate of end   point
                interpolation : str, sampling/interpolation method, including:
                    'nearest'  - nearest neighbour, by default
                    'cubic'    - cubic interpolation
                    'bilinear' - bilinear interpolation

    Returns:    transect: N*2 matrix containing distance and value
                    N is the number of points.
                    1st col - distance in m
                    2nd col - value

    Example: transect = transect_yx(dem,demRsc,[10,15],[100,115])
    """

    # Extract the line
    [y0, x0] = start_yx
    [y1, x1] = end_yx
    length = int(np.hypot(x1-x0, y1-y0))
    x, y = np.linspace(x0, x1, length), np.linspace(y0, y1, length)
    transect = np.zeros([length, 2])

    # Y - Extract the value along the line
    if interpolation.lower() == 'nearest':
        zi = z[np.rint(y).astype(np.int), np.rint(x).astype(np.int)]
    elif interpolation.lower() == 'cubic':
        zi = scipy.ndimage.map_coordinates(z, np.vstack((y, x)))
    elif interpolation.lower() == 'bilinear':
        zi = scipy.ndimage.map_coordinates(z, np.vstack((y, x)), order=2)
    else:
        print('Unrecognized interpolation method: '+interpolation)
        print('Continue with nearest ...')
        zi = z[np.rint(y).astype(np.int), np.rint(x).astype(np.int)]  # nearest neighbour
    transect[:, 1] = zi

    # X - Distance along the line
    earth_radius = 6371.0e3    # in meter
    coord = ut.coordinate(atr)
    try:
        atr['X_FIRST']
        [lat0, lat1] = coord.yx2lalo([y0, y1], coord_type='y')
        lat_c = (lat0 + lat1) / 2.
        x_step = float(atr['X_STEP']) * np.pi/180.0 * earth_radius * np.cos(lat_c * np.pi/180)
        y_step = float(atr['Y_STEP']) * np.pi/180.0 * earth_radius
    except:
        x_step = ut.range_ground_resolution(atr)
        y_step = ut.azimuth_ground_resolution(atr)
    dis_x = (x - x0) * x_step
    dis_y = (y - y0) * y_step
    transect[:, 0] = np.hypot(dis_x, dis_y)

    return transect


def transect_lalo(z, atr, start_lalo, end_lalo, interpolation='nearest'):
    """Extract 2D matrix (z) value along the line [start_lalo, end_lalo]"""
    coord = ut.coordinate(atr)
    [y0, y1], [x0, x1] = coord.lalo2yx([start_lalo[0], end_lalo[0]], [start_lalo[1], end_lalo[1]])
    transect = transect_yx(z, atr, [y0, x0], [y1, x1], interpolation)
    return transect


def transect_list(fileList, inps):
    """Get transection along input line from file list
    Inputs:
        fileList : list of str, path of files to get transect
        inps     : Namespace including the following items:
                   start/end_lalo
                   start/end_yx
                   interpolation
    Outputs:
        transectList : list of N*2 matrix containing distance and its value
        atrList      : list of attribute dictionary, for each input file
    """

    transectList = []
    atrList = []
    for File in fileList:
        print('reading '+File)
        data, atr = readfile.read(File, datasetName=inps.dset)
        if inps.start_lalo and inps.end_lalo:
            transect = transect_lalo(data, atr,
                                     inps.start_lalo, inps.end_lalo,
                                     inps.interpolation)
        else:
            transect = transect_yx(data, atr,
                                   inps.start_yx, inps.end_yx,
                                   inps.interpolation)
        transectList.append(transect)
        atrList.append(atr)
    return transectList, atrList


def get_start_end_points(inps):
    # 1. lonlat file
    if inps.lola_file:
        inps.start_lalo, inps.end_lalo = read_lonlat_file(inps.lola_file)

    # 2. Manually select in window display
    if (not inps.start_yx or not inps.end_yx) and (not inps.start_lalo or not inps.end_lalo):
        print('No input yx/lalo found.')
        print('Continue with manually select start and end point.')
        inps.start_yx, inps.end_yx = manual_select_start_end_point(inps.file[0], dset=inps.dset)

    # Message
    if inps.start_lalo and inps.end_lalo:
        print(f'Start point in lat/lon: {inps.start_lalo}')
        print(f'End   point in lat/lon: {inps.end_lalo}')
    else:
        print(f'Start point in y/x: {inps.start_yx}')
        print(f'End   point in y/x: {inps.end_yx}')
    return


def plot_transect_location(ax, inps):
    print('plot profile line in the 1st input file')
    data0, atr0 = readfile.read(inps.file[0], datasetName=inps.dset)
    ax.imshow(data0)

    coord = ut.coordinate(atr0)
    if inps.start_lalo and inps.end_lalo:
        [y0, y1], [x0, x1] = coord.lalo2yx([inps.start_lalo[0], inps.end_lalo[0]],
                                           [inps.start_lalo[1], inps.end_lalo[1]])
        inps.start_yx = [y0, x0]
        inps.end_yx = [y1, x1]

    ax.plot([inps.start_yx[1], inps.end_yx[1]],
            [inps.start_yx[0], inps.end_yx[0]], 'ro-')
    ax.set_xlim(0, np.shape(data0)[1])
    ax.set_ylim(np.shape(data0)[0], 0)
    ax.set_title('Transect Line in '+inps.file[0])

    # Status bar
    def format_coord(x, y):
        col = int(x)
        row = int(y)
        if 0 <= col < data0.shape[1] and 0 <= row < data0.shape[0]:
            z = data0[row, col]
            if 'X_FIRST' in atr0.keys():
                lat, lon = coord.yx2lalo(row, col)
                return f'lon={lon:.4f}, lat={lat:.4f}, x={x:.0f},  y={y:.0f},  value={z:.4f}'
            else:
                return f'x={x:.0f},  y={y:.0f},  value={z:.4f}'
        else:
            return f'x={x:.0f},  y={y:.0f}'
    ax.format_coord = format_coord
    return


def plot_transect(ax, inps):
    print('plot profiles')

    # disp_unit/scale
    if not inps.disp_unit:
        inps.disp_unit = inps.atrList[0]['UNIT']
    inps.disp_unit, inps.disp_scale = pp.scale_data2disp_unit(data=None,
                                                              metadata=inps.atrList[0],
                                                              disp_unit=inps.disp_unit)[1:3]

    # Plot 2.1 - Input Files
    value_min = 0
    value_max = 0
    for i in range(len(inps.file)):
        # Profile Color based on Asc/Desc direction
        if inps.atrList[i]['ORBIT_DIRECTION'][0].upper() == 'A':
            p_color = 'crimson'
        else:
            p_color = 'royalblue'
        # Plot
        distance = inps.transectList[i][:, 0]/1000.0          # km
        value = inps.transectList[i][:, 1]*inps.disp_scale - inps.disp_offset*i
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
    ax.set_ylabel(f'Mean LOS Velocity ({inps.disp_unit})', fontsize=inps.font_size)
    # X axis
    ax.set_xlabel('Distance (km)', fontsize=inps.font_size)
    ax.tick_params(which='both', direction='out', labelsize=inps.font_size)

    # Plot 2.2 - DEM
    if inps.dem:
        ax2 = ax.twinx()
        distance = inps.demTransectList[0][:, 0]/1000.0    # km
        value = inps.demTransectList[0][:, 1]/1000.0    # km
        ax2.fill_between(distance, 0, value, facecolor='gray')

        # Y axis - display DEM in the bottom
        value_min = np.nanmin(value)
        value_max = np.nanmax(value)
        if not inps.dem_disp_min:
            inps.dem_disp_min = np.floor(value_min*2.0)/2.0
        if not inps.dem_disp_max:
            inps.dem_disp_max = np.ceil((value_max +
                                         (value_max-value_min)*(len(inps.file)+0.0))*2.0)/2.0
        ax2.set_ylim(inps.dem_disp_min, inps.dem_disp_max)
        # Show lower part of yaxis
        #dem_tick = ax2.yaxis.get_majorticklocs()
        #dem_tick = dem_tick[:len(dem_tick)/2]
        # ax2.set_yticks(dem_tick)
        ax2.set_ylabel('Elevation (km)', fontsize=inps.font_size)
        ax2.tick_params(which='both', direction='out', labelsize=inps.font_size)

    # X axis - Shared
    distanceMax = np.nanmax(inps.transectList[0][:, 0]/1000.0)   # in km
    plt.xlim(0, distanceMax)
    plt.tight_layout()
    return


def save_transect(fig_list, inps):
    # Output file name
    if not inps.outfile:
        figBase = 'transect_x{}y{}_x{}y{}'.format(inps.start_yx[1],
                                                  inps.start_yx[0],
                                                  inps.end_yx[1],
                                                  inps.end_yx[0])
    else:
        figBase, inps.outfile_ext = os.path.splitext(inps.outfile)
        if not inps.outfile_ext:
            inps.outfile_ext = '.png'

    # save figure
    print('writing >>> '+figBase+inps.outfile_ext)
    fig_list[0].savefig(inps.file[0]+inps.outfile_ext,
                        bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)
    fig_list[1].savefig(figBase+inps.outfile_ext,
                        bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)

    # save data to .mat file
    print('writing >>> '+figBase+'.mat')
    transect_mat = {}
    for i in range(len(inps.file)):
        project = inps.atrList[i]['PROJECT_NAME']
        transect_mat[project] = inps.transectList[i]
    if inps.dem:
        transect_mat['elevation'] = inps.demTransectList[0]
    sio.savemat(figBase+'.mat', {'transection': transect_mat})
    return


############################ Main ###################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    print('\n**************** Transect *********************')
    print(f'input files: ({len(inps.file)})\n{inps.file}')

    get_start_end_points(inps)

    print('extract transect from input files ...')
    inps.transectList, inps.atrList = transect_list(inps.file, inps)
    if inps.dem:
        inps.demTransectList, demAtrList = transect_list([inps.dem], inps)
    print('profile length: '+str(inps.transectList[0][-1, 0]/1000.0)+' km')

    fig0, ax0 = plt.subplots()
    plot_transect_location(ax0, inps)

    fig, ax = plt.subplots(figsize=inps.fig_size)
    plot_transect(ax, inps)

    # save figure and data
    if inps.save_fig:
        save_transect(fig_list=[fig0, fig], inps=inps)

    # Display
    if inps.disp_fig:
        print('showing ...')
        plt.show()
    return


#####################################################################
if __name__ == '__main__':
    main()
