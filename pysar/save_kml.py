#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# Yunjun, Jul 2015: add 'timeseries'/'wrapped' option
# Yunjun, Oct 2015: merge all HDF5 option into one
#                   add support for ROI_PAC product
# Yunjun, Nov 2015: support different fig unit
#                   update colorbar
# Yunjun, Jul 2017: re-write using argparse and more pysar module


import os
import sys
import argparse

try:
    from pykml.factory import KML_ElementMaker as KML
except:
    sys.exit('pykml should be installed!')

from lxml import etree
import h5py
import numpy as np
import matplotlib as mpl;  mpl.use('Agg')
import matplotlib.pyplot as plt

import pysar
import pysar._datetime as ptime
import pysar._readfile as readfile
import pysar._pysar_utilities as ut
import pysar.view as pview
from pysar._readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file


############################################################
def write_kmz_file(data, atr, out_name_base, inps=None):
    ''' Generate Google Earth KMZ file for input data matrix.
    Inputs:
        data - 2D np.array in int/float, data matrix to write
        out_name_base - string, output file name base
        atr  - dict, containing the following attributes:
               WIDTH/FILE_LENGTH : required, file size
               X/Y_FIRST/STEP    : required, for lat/lon spatial converage
               ref_x/y           : optional, column/row number of reference pixel
               PROJECT_NAME      : optional, for KMZ folder name
        inps - Namespace, optional, input options for display
    Output:
        kmz_file - string, output KMZ filename
    Example:
        import pysar._readfile as readfile
        import pysar.view as pview
        import pysar.save_kml as save_kml
        fname = 'geo_velocity_masked.h5'
        data, atr = readfile.read(fname)
        out_name_base = pview.auto_figure_title(fname, None)
        save_kml.write_kmz_file(data, atr, out_name_base)
    '''
    if not inps:
        inps = cmdLineParse()

    if not inps.ylim:
        inps.ylim = [np.nanmin(data), np.nanmax(data)]

    west, east, south, north = ut.four_corners(atr)

    ## 2.1 Make PNG file - Data
    print 'plotting data ...'

    # Figure size
    if not inps.fig_size:
        fig_scale = min(pysar.figsize_single_min/min(data.shape),\
                        pysar.figsize_single_max/max(data.shape))
        inps.fig_size = [np.rint(i*fig_scale*2)/2 for i in data.shape]
    print 'create figure in size: '+str(inps.fig_size)
    fig = plt.figure(figsize=inps.fig_size, frameon=False)
    ax = fig.add_axes([0., 0., 1., 1.])
    ax.set_axis_off()

    print 'colormap: '+inps.colormap
    inps.colormap = plt.get_cmap(inps.colormap)

    # Plot - data matrix
    ax.imshow(data, aspect='auto', cmap=inps.colormap, vmin=inps.ylim[0], vmax=inps.ylim[1])

    # Plot - reference pixel
    if inps.disp_seed == 'yes':
        try:
            xref = int(atr['ref_x'])
            yref = int(atr['ref_y'])
            ax.plot(xref, yref, 'ks', ms=inps.seed_size)
            print 'show reference point'
        except:
            inps.disp_seed = False
            print 'Cannot find reference point info!'

    width = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])
    ax.set_xlim([0,width])
    ax.set_ylim([length,0])

    data_png_file = out_name_base + '.png'
    print 'writing '+data_png_file
    plt.savefig(data_png_file, pad_inches=0.0, transparent=True, dpi=inps.fig_dpi)

    ## 2.2 Making PNG file - colorbar
    pc = plt.figure(figsize=(1,8))
    cax = pc.add_subplot(111)
    norm = mpl.colors.Normalize(vmin=inps.ylim[0], vmax=inps.ylim[1])
    cbar = mpl.colorbar.ColorbarBase(cax, cmap=inps.colormap, norm=norm, orientation='vertical')

    cbar.set_label(inps.cbar_label+' ['+inps.disp_unit+']')
    cbar.locator = mpl.ticker.MaxNLocator(nbins=inps.cbar_bin_num)
    cbar.update_ticks()

    pc.subplots_adjust(left=0.2,bottom=0.3,right=0.4,top=0.7)
    pc.patch.set_facecolor('white')
    pc.patch.set_alpha(0.7)

    cbar_png_file = out_name_base + '_cbar.png'
    print 'writing '+cbar_png_file
    pc.savefig(cbar_png_file, bbox_inches='tight', facecolor=pc.get_facecolor(), dpi=inps.fig_dpi)

    ## 2.3 Generate KML file
    print 'generating kml file ...'
    try:     doc = KML.kml(KML.Folder(KML.name(atr['PROJECT_NAME'])))
    except:  doc = KML.kml(KML.Folder(KML.name('PySAR product')))

    # Add data png file
    slc = KML.GroundOverlay(KML.name(data_png_file), KML.Icon(KML.href(data_png_file)),\
                            KML.altitudeMode('clampToGround'),\
                            KML.LatLonBox(KML.north(str(north)), KML.east(str(east)),\
                                          KML.south(str(south)), KML.west(str(west))))
    doc.Folder.append(slc)

    # Add colorbar png file
    cb_rg = min(north-south, east-west)
    cb_N = (north+south)/2.0 + 0.5*0.5*cb_rg
    cb_W = east  + 0.1*cb_rg

    if inps.cbar_height:
        slc1 = KML.GroundOverlay(KML.name('colorbar'), KML.Icon(KML.href(cbar_png_file)),\
                                 KML.altitude(str(inps.cbar_height)),KML.altitudeMode('absolute'),\
                                 KML.LatLonBox(KML.north(str(cb_N)),KML.south(str(cb_N-0.5*cb_rg)),\
                                               KML.west( str(cb_W)),KML.east( str(cb_W+0.14*cb_rg))))
    else:
        slc1 = KML.GroundOverlay(KML.name('colorbar'), KML.Icon(KML.href(cbar_png_file)),\
                                 KML.altitudeMode('clampToGround'),\
                                 KML.LatLonBox(KML.north(str(cb_N)),KML.south(str(cb_N-0.5*cb_rg)),\
                                               KML.west( str(cb_W)),KML.east( str(cb_W+0.14*cb_rg))))
    doc.Folder.append(slc1)

    # Write KML file
    kmlstr = etree.tostring(doc, pretty_print=True) 
    kml_file = out_name_base + '.kml'
    print 'writing '+kml_file
    f = open(kml_file, 'w')
    f.write(kmlstr)
    f.close()

    ## 2.4 Generate KMZ file
    kmz_file = out_name_base + '.kmz'
    print 'writing '+kmz_file
    cmdKMZ = 'zip '+kmz_file+' '+kml_file+' '+data_png_file+' '+cbar_png_file
    os.system(cmdKMZ)

    cmdClean = 'rm '+kml_file;         print cmdClean;    os.system(cmdClean)
    cmdClean = 'rm '+data_png_file;    print cmdClean;    os.system(cmdClean)
    cmdClean = 'rm '+cbar_png_file;    print cmdClean;    os.system(cmdClean)

    return kmz_file


############################################################
EXAMPLE='''example:
  save_kml.py geo_velocity_masked.h5 
  save_kml.py geo_timeseries_masked.h5  20101120
  save_kml.py geo_unwrapIfgram.h5       101120-110220

  save_kml.py geo_velocity_masked.h5 -u cm --ylim -2 2
  save_kml.py geo_velocity_masked.h5 -u cm --ylim -2.5 0.5 -c jet_r
  save_kml.py geo_velocity_masked.h5 -u cm --ylim -2 2 --ref-size 3 --fig-size 5 8
  save_kml.py demGeo.h5 --cbar-label Elevation
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Generate Google Earth KMZ file.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='file to be converted, in geo coordinate.')
    parser.add_argument('epoch', nargs='?', help='date of timeseries, or date12 of interferograms to be converted')
    parser.add_argument('-o','--output', dest='outfile', help='output file base name. Extension is fixed with .kmz')

    parser.add_argument('--ylim', dest='ylim', nargs=2, metavar=('MIN','MAX'), type=float,\
                        help='Y/value limits for plotting.')
    parser.add_argument('-u', dest='disp_unit', metavar='UNIT',\
                        help='unit for display.')
    parser.add_argument('-c','--cm','--colormap', dest='colormap', default='jet',\
                        help='Colormap for plotting. Default: jet')
    parser.add_argument('--wrap', action='store_true', help='re-wrap data to display data in fringes.')

    # Figure
    fig = parser.add_argument_group('Figure')
    fig.add_argument('--cbar-bin-num', dest='cbar_bin_num', metavar='NUM', type=int, default=9,\
                     help='Colorbar bin number. Default: 9')
    fig.add_argument('--cbar-label', dest='cbar_label', metavar='LABEL', default='Mean LOS velocity',\
                     help='Colorbar label. Default: Mean LOS velocity')
    fig.add_argument('--cbar-height', dest='cbar_height', metavar='NUM', type=float,\
                     help='Colorbar height/elevation/altitude in meters. clampToGround if not specified.')
    fig.add_argument('--dpi', dest='fig_dpi', metavar='NUM', type=int, default=300,\
                     help='Figure DPI (dots per inch). Default: 300')
    fig.add_argument('--figsize', dest='fig_size', metavar=('WID','LEN'), type=float, nargs=2,\
                     help='Figure size in inches - width and length')

    # Reference Pixel
    ref = parser.add_argument_group('Reference Pixel')
    ref.add_argument('--noreference', dest='disp_seed', action='store_false', help='do not show reference point')
    ref.add_argument('--ref-color', dest='seed_color', metavar='COLOR', default='k',\
                     help='marker color of reference point')
    ref.add_argument('--ref-size', dest='seed_size', metavar='NUM', type=int, default=5,\
                     help='marker size of reference point, default: 10')
    ref.add_argument('--ref-symbol', dest='seed_symbol', metavar='SYMBOL', default='s',\
                     help='marker symbol of reference point')

    inps = parser.parse_args()
    return inps


############################################################
def main(argv):
    inps = cmdLineParse()

    ##### 1. Read data
    atr = readfile.read_attribute(inps.file)
    k = atr['FILE_TYPE']
    print 'Input file is '+k

    # Check: file in geo coord
    if 'X_FIRST' not in atr.keys():
        sys.exit('ERROR: Input file is not geocoded.')

    # Check: epoch is required for multi_dataset/group files
    if not inps.epoch and k in multi_group_hdf5_file+multi_dataset_hdf5_file:
        print "No date/date12 input.\nIt's required for "+k+" file"
        sys.exit(1)

    # Read data
    data, atr = readfile.read(inps.file, (), inps.epoch)

    # Output filename
    if not inps.outfile:
        inps.outfile = pview.auto_figure_title(inps.file, inps.epoch, vars(inps))

    # Data Operation - Display Unit & Rewrapping
    data, inps.disp_unit, inps.wrap = pview.scale_data4disp_unit_and_rewrap(data, atr, inps.disp_unit, inps.wrap)
    if inps.wrap:
        inps.ylim = [-np.pi, np.pi]

    ##### 2. Generate Google Earth KMZ
    kmz_file = write_kmz_file(data, atr, inps.outfile, inps)

    print 'Done.'
    return

#######################################################
if __name__ == '__main__':
    main(sys.argv[1:])

