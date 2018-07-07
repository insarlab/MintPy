#!/usr/bin/env python
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2018, Heresh Fattahi, Zhang Yunjun     #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################

import os
import sys
import argparse

try:
    from pykml.factory import KML_ElementMaker as KML
except ImportError:
    raise ImportError('Can not import pykml!')

from lxml import etree
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from pysar.objects import timeseriesKeyNames
from pysar.utils import readfile, ptime, utils as ut, plot as pp


############################################################
EXAMPLE = """example:
  save_kml.py geo_velocity_masked.h5 
  save_kml.py geo_timeseries_masked.h5  20101120
  save_kml.py geo_unwrapIfgram.h5       101120-110220

  save_kml.py geo_velocity_masked.h5 -u cm --ylim -2 2
  save_kml.py geo_velocity_masked.h5 -u cm --ylim -2.5 0.5 -c jet_r
  save_kml.py geo_velocity_masked.h5 -u cm --ylim -2 2 --ref-size 3 --fig-size 5 8
  save_kml.py demGeo.h5 --cbar-label Elevation
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Generate Google Earth KMZ file.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='file to be converted, in geo coordinate.')
    parser.add_argument('dset', nargs='?',
                        help='date of timeseries, or date12 of interferograms to be converted')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output file base name. Extension is fixed with .kmz')

    parser.add_argument('--ylim', dest='ylim', nargs=2, metavar=('MIN', 'MAX'), type=float,
                        help='Y/value limits for plotting.')
    parser.add_argument('-u', dest='disp_unit', metavar='UNIT',
                        help='unit for display.')
    parser.add_argument('-c', '--cm', '--colormap', dest='colormap', default='jet',
                        help='Colormap for plotting. Default: jet')
    parser.add_argument('--wrap', action='store_true',
                        help='re-wrap data to display data in fringes.')

    # Figure
    fig = parser.add_argument_group('Figure')
    fig.add_argument('--cbar-bin-num', dest='cbar_bin_num', metavar='NUM', type=int, default=9,
                     help='Colorbar bin number. Default: 9')
    fig.add_argument('--cbar-label', dest='cbar_label', metavar='LABEL', default='Mean LOS velocity',
                     help='Colorbar label. Default: Mean LOS velocity')
    fig.add_argument('--cbar-height', dest='cbar_height',
                     help='Colorbar height/elevation/altitude in meters;\n' +
                          'if not specified and DEM exists in current directory, use mean DEM height + 1000m;\n' +
                          'if not specified nor DEM exists, clampToGround.')
    fig.add_argument('--dpi', dest='fig_dpi', metavar='NUM', type=int, default=300,
                     help='Figure DPI (dots per inch). Default: 300')
    fig.add_argument('--figsize', dest='fig_size', metavar=('WID', 'LEN'), type=float, nargs=2,
                     help='Figure size in inches - width and length')

    # Reference Pixel
    ref = parser.add_argument_group('Reference Pixel')
    ref.add_argument('--noreference', dest='disp_seed',
                     action='store_false', help='do not show reference point')
    ref.add_argument('--ref-color', dest='seed_color', metavar='COLOR', default='k',
                     help='marker color of reference point')
    ref.add_argument('--ref-size', dest='seed_size', metavar='NUM', type=int, default=5,
                     help='marker size of reference point, default: 10')
    ref.add_argument('--ref-symbol', dest='seed_symbol', metavar='SYMBOL', default='s',
                     help='marker symbol of reference point')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    atr = readfile.read_attribute(inps.file)
    # Check: file in geo coord
    if 'X_FIRST' not in atr.keys():
        raise Exception('ERROR: Input file is not geocoded.')

    # Check: dset is required for multi_dataset/group files
    if not inps.dset and atr['FILE_TYPE'] in ['ifgramStack']+timeseriesKeyNames:
        raise Exception("No date/date12 input.\nIt's required for "+k+" file")

    return inps


############################################################
def write_kmz_file(data, metadata, out_file, inps=None):
    """ Generate Google Earth KMZ file for input data matrix.
    Inputs:
        data - 2D np.array in int/float, data matrix to write
        out_file - string, output file name
        metadata  - dict, containing the following attributes:
               WIDTH/LENGTH      : required, file size
               X/Y_FIRST/STEP    : required, for lat/lon spatial converage
               ref_x/y           : optional, column/row number of reference pixel
               PROJECT_NAME      : optional, for KMZ folder name
        inps - Namespace, optional, input options for display
    Output:
        kmz_file - string, output KMZ filename
    Example:
        from pysar.utils import readfile, plot as pp
        from pysar import save_kml
        fname = 'geo_velocity_masked.h5'
        data, atr = readfile.read(fname)
        out_file = pp.auto_figure_title(fname, None)+'.kmz'
        save_kml.write_kmz_file(data, atr, out_file)
    """
    if not inps:
        inps = cmd_line_parse()

    if not inps.ylim:
        inps.ylim = [np.nanmin(data), np.nanmax(data)]

    west, east, south, north = ut.four_corners(metadata)

    # 2.1 Make PNG file - Data
    print('plotting data ...')

    # Figure size
    if not inps.fig_size:
        plot_shape = [east-west, north-south]
        fig_scale = min(pp.min_figsize_single / min(plot_shape),
                        pp.max_figsize_single / max(plot_shape),
                        pp.max_figsize_height / plot_shape[1])
        inps.fig_size = [np.floor(i*fig_scale*2)/2 for i in plot_shape]
    print('create figure in size: '+str(inps.fig_size))
    fig = plt.figure(figsize=inps.fig_size, frameon=False)
    ax = fig.add_axes([0., 0., 1., 1.])
    ax.set_axis_off()

    print('colormap: '+inps.colormap)
    inps.colormap = plt.get_cmap(inps.colormap)

    # Plot - data matrix
    ax.imshow(data, cmap=inps.colormap,
              vmin=inps.ylim[0], vmax=inps.ylim[1],
              aspect='auto', interpolation='nearest')

    # Plot - reference pixel
    if inps.disp_seed == 'yes':
        try:
            xref = int(metadata['REF_X'])
            yref = int(metadata['REF_Y'])
            ax.plot(xref, yref, 'ks', ms=inps.seed_size)
            print('show reference point')
        except:
            inps.disp_seed = False
            print('Cannot find reference point info!')

    width = int(metadata['WIDTH'])
    length = int(metadata['LENGTH'])
    ax.set_xlim([0, width])
    ax.set_ylim([length, 0])

    out_name_base = os.path.splitext(out_file)[0]
    data_png_file = out_name_base + '.png'
    print('writing {} with dpi={}'.format(data_png_file, inps.fig_dpi))
    plt.savefig(data_png_file, pad_inches=0.0,
                transparent=True, dpi=inps.fig_dpi)

    # 2.2 Making PNG file - colorbar
    pc = plt.figure(figsize=(1, 8))
    cax = pc.add_subplot(111)
    norm = mpl.colors.Normalize(vmin=inps.ylim[0], vmax=inps.ylim[1])
    cbar = mpl.colorbar.ColorbarBase(cax, cmap=inps.colormap,
                                     norm=norm, orientation='vertical')

    cbar.set_label('{} [{}]'.format(inps.cbar_label, inps.disp_unit))
    cbar.locator = mpl.ticker.MaxNLocator(nbins=inps.cbar_bin_num)
    cbar.update_ticks()

    pc.subplots_adjust(left=0.2, bottom=0.3, right=0.4, top=0.7)
    pc.patch.set_facecolor('white')
    pc.patch.set_alpha(0.7)

    cbar_png_file = '{}_cbar.png'.format(out_name_base)
    print('writing '+cbar_png_file)
    pc.savefig(cbar_png_file, bbox_inches='tight',
               facecolor=pc.get_facecolor(), dpi=inps.fig_dpi)

    # 2.3 Generate KML file
    print('generating kml file ...')
    try:
        doc = KML.kml(KML.Folder(KML.name(metadata['PROJECT_NAME'])))
    except:
        doc = KML.kml(KML.Folder(KML.name('PySAR product')))

    # Add data png file
    slc = KML.GroundOverlay(KML.name(data_png_file),
                            KML.Icon(KML.href(data_png_file)),
                            KML.altitudeMode('clampToGround'),
                            KML.LatLonBox(KML.north(str(north)),
                                          KML.east(str(east)),
                                          KML.south(str(south)),
                                          KML.west(str(west))))
    doc.Folder.append(slc)

    # Add colorbar png file
    cb_rg = min(north - south, east - west)
    cb_N = (north + south) / 2.0 + 0.5 * 0.5 * cb_rg
    cb_W = east + 0.1*cb_rg

    # Use mean height from existed DEM file
    if not inps.cbar_height:
        try:
            fileList = ['geo_geometry*.h5', 'INPUTS/geometry*.h5',
                        'dem*.h5', '*.dem', 'radar*.hgt']
            dem_file = ut.get_file_list(fileList)[0]
            print('use mean height from file: {} + 1000 m as colorbar height.'.format(dem_file))
            dem_data = readfile.read(dem_file, datasetName='height')[0]
            inps.cbar_height = np.rint(np.nanmean(dem_data)) + 1000.0
        except:
            pass
    elif str(inps.cbar_height).lower().endswith('ground'):
        inps.cbar_height = None

    if inps.cbar_height:
        print('set colorbar in height: %.2f m' % inps.cbar_height)
        slc1 = KML.GroundOverlay(KML.name('colorbar'),
                                 KML.Icon(KML.href(cbar_png_file)),
                                 KML.altitude(str(inps.cbar_height)),
                                 KML.altitudeMode('absolute'),
                                 KML.LatLonBox(KML.north(str(cb_N)),
                                               KML.south(str(cb_N-0.5*cb_rg)),
                                               KML.west(str(cb_W)),
                                               KML.east(str(cb_W+0.14*cb_rg))))
    else:
        print('set colorbar clampToGround')
        slc1 = KML.GroundOverlay(KML.name('colorbar'),
                                 KML.Icon(KML.href(cbar_png_file)),
                                 KML.altitudeMode('clampToGround'),
                                 KML.LatLonBox(KML.north(str(cb_N)),
                                               KML.south(str(cb_N-0.5*cb_rg)),
                                               KML.west(str(cb_W)),
                                               KML.east(str(cb_W+0.14*cb_rg))))
    doc.Folder.append(slc1)
    kmlstr = etree.tostring(doc, pretty_print=True).decode('utf8')

    # Write KML file
    kml_file = '{}.kml'.format(out_name_base)
    print('writing '+kml_file)
    with open(kml_file, 'w') as f:
        f.write(kmlstr)

    # 2.4 Generate KMZ file
    kmz_file = '{}.kmz'.format(out_name_base)
    cmdKMZ = 'zip {} {} {} {}'.format(kmz_file, kml_file, data_png_file, cbar_png_file)
    print('writing {}\n{}'.format(kmz_file, cmdKMZ))
    os.system(cmdKMZ)

    cmdClean = 'rm {} {} {}'.format(kml_file, data_png_file, cbar_png_file)
    print(cmdClean)
    os.system(cmdClean)

    return kmz_file


############################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    plt.switch_backend('Agg')  # Backend setting
    #print("The Python version is %s.%s.%s" % sys.version_info[:3])

    # Read data
    data, atr = readfile.read(inps.file, datasetName=inps.dset)

    # Data Operation - Display Unit & Rewrapping
    (data,
     inps.disp_unit,
     inps.disp_scale,
     inps.wrap) = pp.scale_data4disp_unit_and_rewrap(data=data,
                                                     metadata=atr,
                                                     disp_unit=inps.disp_unit,
                                                     wrap=inps.wrap)
    if inps.wrap:
        inps.ylim = [-np.pi, np.pi]

    # Output filename
    if not inps.outfile:
        inps.outfile = '{}.kmz'.format(pp.auto_figure_title(inps.file,
                                                            datasetNames=inps.dset,
                                                            inps_dict=vars(inps)))

    # 2. Generate Google Earth KMZ
    kmz_file = write_kmz_file(data,
                              metadata=atr,
                              out_file=inps.outfile,
                              inps=inps)

    print('Done.')
    return


#######################################################
if __name__ == '__main__':
    main()
