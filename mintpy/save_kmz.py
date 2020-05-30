#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


try:
    from pykml.factory import KML_ElementMaker as KML
except ImportError:
    raise ImportError('Can not import pykml!')

import os
import argparse
import numpy as np
from lxml import etree
from zipfile import ZipFile
import matplotlib as mpl
from matplotlib import pyplot as plt
from mintpy.objects import timeseriesKeyNames
from mintpy.utils import readfile, utils as ut, plot as pp


############################################################
EXAMPLE = """example:
  save_kmz.py geo/geo_velocity.h5 
  save_kmz.py geo/geo_velocity.h5 -u cm -v -2 2
  save_kmz.py geo/geo_velocity.h5 -u cm --wrap --wrap-range -3 7

  save_kmz.py geo/geo_timeseries_ERA5_ramp_demErr.h5 20101120
  save_kmz.py geo/geo_timeseries_ERA5_demErr.h5 20200505_20200517

  save_kmz.py geo/geo_ifgramStack.h5 20101120_20110220
  save_kmz.py geo/geo_geometryRadar.h5 --cbar-label Elevation
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Generate Google Earth KMZ file with raster image.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='file to be converted, in geo coordinate.')
    parser.add_argument('dset', nargs='?',
                        help='date of timeseries, or date12 of interferograms to be converted')
    parser.add_argument('-m','--mask', dest='mask_file', metavar='FILE',
                        help='mask file for display')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output file base name. Extension is fixed with .kmz')

    # Data
    parser.add_argument('-v','--vlim', dest='vlim', nargs=2, metavar=('MIN', 'MAX'), type=float,
                        help='Y/value limits for plotting.')
    parser.add_argument('-u', dest='disp_unit', metavar='UNIT',
                        help='unit for display.')
    parser.add_argument('-c', '--cm', '--colormap', dest='colormap', default='jet',
                        help='Colormap for plotting. Default: jet')
    parser.add_argument('--wrap', action='store_true',
                        help='re-wrap data to display data in fringes.')
    parser.add_argument('--wrap-range', dest='wrap_range', type=float, nargs=2,
                        default=[-1.*np.pi, np.pi], metavar=('MIN', 'MAX'),
                        help='range of one cycle after wrapping, default: [-pi, pi]')

    # Figure
    fig = parser.add_argument_group('Figure')
    fig.add_argument('--dpi', dest='fig_dpi', metavar='NUM', type=int, default=600,
                     help='Figure DPI (dots per inch). Default: 600')
    fig.add_argument('--figsize', dest='fig_size', metavar=('WID', 'LEN'), type=float, nargs=2,
                     help='Figure size in inches - width and length')
    fig.add_argument('--cbar-loc', dest='cbar_loc', choices=['lower left','lower right','upper left', 'upper right'],
                     default='lower left', help='Location of colorbar in the screen. Default: lower left.')
    fig.add_argument('--cbar-label', dest='cbar_label', metavar='LABEL', default='Mean LOS velocity',
                     help='Colorbar label. Default: Mean LOS velocity')
    fig.add_argument('--cbar-bin-num', dest='cbar_bin_num', metavar='NUM', type=int, default=7,
                     help='Colorbar bin number. Default: 7')

    # Reference Pixel
    ref = parser.add_argument_group('Reference Pixel')
    ref.add_argument('--noreference', dest='disp_ref_pixel',
                     action='store_false', help='do not show reference point')
    ref.add_argument('--ref-color', dest='ref_marker_color', metavar='COLOR', default='k',
                     help='marker color of reference point')
    ref.add_argument('--ref-size', dest='ref_marker_size', metavar='NUM', type=int, default=5,
                     help='marker size of reference point, default: 10')
    ref.add_argument('--ref-marker', dest='ref_marker', metavar='SYMBOL', default='s',
                     help='marker symbol of reference point')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    atr = readfile.read_attribute(inps.file)
    # Check 1: file in geo coord
    if 'X_FIRST' not in atr.keys():
        raise Exception('ERROR: Input file is not geocoded.')

    # Check 2: dset is required for multi_dataset/group files
    if not inps.dset and atr['FILE_TYPE'] in ['ifgramStack']+timeseriesKeyNames:
        raise Exception("No date/date12 input.\nIt's required for {} file".format(atr['FILE_TYPE']))

    return inps


def plot_colorbar(out_file, vmin, vmax, unit='cm/year', cmap='jet', figsize=(0.18, 3.6),
                  nbins=7, label='Mean LOS velocity'):
    fig, cax = plt.subplots(figsize=figsize)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cbar = mpl.colorbar.ColorbarBase(cax, cmap=plt.get_cmap(cmap), norm=norm, orientation='vertical')
    cbar.set_label('{} [{}]'.format(label, unit), fontsize=12)
    cbar.locator = mpl.ticker.MaxNLocator(nbins=nbins)
    cbar.update_ticks()
    cbar.ax.tick_params(which='both', labelsize=12)
    fig.patch.set_facecolor('white')
    fig.patch.set_alpha(0.7)
    print('writing', out_file)
    fig.savefig(out_file, bbox_inches='tight', facecolor=fig.get_facecolor(), dpi=300)
    return out_file


def generate_cbar_element(cbar_file, vmin, vmax, unit='cm/year', cmap='jet', loc='lower left',
                          nbins=7, label='Mean LOS velocity'):
    # plot colobar and save as an image
    cbar_file = plot_colorbar(out_file=cbar_file,
                              vmin=vmin,
                              vmax=vmax,
                              unit=unit,
                              cmap=cmap,
                              nbins=nbins,
                              label=label)

    # colobar location
    if loc.split()[0] == 'lower':
        oy, sy = '0', '0'
    elif loc.split()[0] == 'upper':
        oy, sy = '1', '1'
    if loc.split()[1] == 'left':
        ox, sx = '0', '0'
    elif loc.split()[1] == 'right':
        ox, sx = '1', '1'

    # generate KML screen overlay object
    cbar_overlay = KML.ScreenOverlay(
        KML.name('colorbar'),
        KML.Icon(
            KML.href("{}".format(os.path.basename(cbar_file))),
            KML.viewBoundScale(0.75)
        ),
        KML.overlayXY(x=ox, y=oy, xunits="fraction", yunits="fraction"),
        KML.screenXY(x=sx, y=sy, xunits="fraction", yunits="fraction"),
        KML.size(x="0", y="250", xunits="pixel", yunits="pixel"),
        KML.rotation(0),
        KML.visibility(1),
        KML.open(0)
    )
    #print('add colorbar.')
    return cbar_overlay


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
        from mintpy.utils import readfile, plot as pp
        from mintpy import save_kmz
        fname = 'geo_velocity_masked.h5'
        data, atr = readfile.read(fname)
        out_file = pp.auto_figure_title(fname, None)+'.kmz'
        save_kmz.write_kmz_file(data, atr, out_file)
    """
    if not inps:
        inps = cmd_line_parse()

    if not inps.vlim:
        inps.vlim = [np.nanmin(data), np.nanmax(data)]

    west, east, south, north = ut.four_corners(metadata)

    # 2.1 Make PNG file - Data
    print('plotting data ...')

    # Figure size
    if not inps.fig_size:
        plot_shape = [east-west, north-south]
        fig_scale = min(pp.min_figsize_single / min(plot_shape),
                        pp.max_figsize_single / max(plot_shape),
                        pp.max_figsize_height / plot_shape[1])
        inps.fig_size = [2.*i*fig_scale for i in plot_shape]
    print('create figure in size: '+str(inps.fig_size))
    fig = plt.figure(figsize=inps.fig_size, frameon=False)
    ax = fig.add_axes([0., 0., 1., 1.])
    ax.set_axis_off()

    inps.colormap = pp.check_colormap_input(metadata, inps.colormap)
    # Plot - data matrix
    ax.imshow(data, cmap=inps.colormap,
              vmin=inps.vlim[0], vmax=inps.vlim[1],
              aspect='auto', interpolation='nearest')

    # Plot - reference pixel
    if inps.disp_ref_pixel:
        try:
            xref = int(metadata['REF_X'])
            yref = int(metadata['REF_Y'])
            ax.plot(xref, yref, inps.ref_marker,
                    color=inps.ref_marker_color,
                    ms=inps.ref_marker_size)
            print('show reference point')
        except:
            inps.disp_ref_pixel = False
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

    # 2.3 Generate KML file
    proj_name = metadata.get('PROJECT_NAME', 'InSAR_MintPy')
    doc = KML.kml(KML.Folder(KML.name(proj_name)))

    # Add data png file
    img_name = os.path.splitext(os.path.basename(data_png_file))[0]
    img = KML.GroundOverlay(
        KML.name(img_name),
        KML.Icon(
            KML.href(os.path.basename(data_png_file))
        ),
        KML.altitudeMode('clampToGround'),
        KML.LatLonBox(KML.north(str(north)),
                      KML.east(str(east)),
                      KML.south(str(south)),
                      KML.west(str(west)))
    )
    doc.Folder.append(img)

    # Add colorbar png file
    cbar_file = '{}_cbar.png'.format(out_name_base)
    cbar_overlay = generate_cbar_element(cbar_file,
                                         vmin=inps.vlim[0],
                                         vmax=inps.vlim[1],
                                         unit=inps.disp_unit,
                                         cmap=inps.colormap,
                                         loc=inps.cbar_loc,
                                         nbins=inps.cbar_bin_num,
                                         label=inps.cbar_label)
    doc.Folder.append(cbar_overlay)
    kmlstr = etree.tostring(doc, pretty_print=True).decode('utf8')

    # Write KML file
    kml_file = '{}.kml'.format(out_name_base)
    print('writing '+kml_file)
    with open(kml_file, 'w') as f:
        f.write(kmlstr)

    # 2.4 Generate KMZ file, by
    # 1) go to the directory of kmz file
    run_dir = os.path.abspath(os.getcwd())
    os.chdir(os.path.abspath(os.path.dirname(out_file)))
    # 2) zip all data files
    kmz_file = '{}.kmz'.format(out_name_base)
    with ZipFile(kmz_file, 'w') as fz:
        for fname in [kml_file, data_png_file, cbar_file]:
            fz.write(os.path.relpath(fname))
            os.remove(fname)
    # 3) go back to the running directory
    os.chdir(run_dir)
    print('merged all files to {}'.format(kmz_file))
    return kmz_file


############################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    plt.switch_backend('Agg')  # Backend setting

    # Read data
    k = readfile.read_attribute(inps.file)['FILE_TYPE']
    if k == 'timeseries' and inps.dset and '_' in inps.dset:
        inps.ref_date, inps.dset = inps.dset.split('_')
    else:
        inps.ref_date = None

    data, atr = readfile.read(inps.file, datasetName=inps.dset)
    if k == 'timeseries' and inps.ref_date:
        data -= readfile.read(inps.file, datasetName=inps.ref_date)[0]

    # mask
    mask = pp.read_mask(inps.file, mask_file=inps.mask_file, datasetName=inps.dset, print_msg=True)[0]
    if mask is not None:
        data = np.ma.masked_where(mask == 0., data)

    # Data Operation - Display Unit & Rewrapping
    (data,
     inps.disp_unit,
     inps.disp_scale,
     inps.wrap) = pp.scale_data4disp_unit_and_rewrap(data,
                                                     metadata=atr,
                                                     disp_unit=inps.disp_unit,
                                                     wrap=inps.wrap,
                                                     wrap_range=inps.wrap_range)
    if inps.wrap:
        inps.vlim = inps.wrap_range

    # Output filename
    inps.fig_title = pp.auto_figure_title(inps.file,
                                          datasetNames=inps.dset,
                                          inps_dict=vars(inps))
    if not inps.outfile:
        inps.outfile = '{}.kmz'.format(inps.fig_title)
    inps.outfile = os.path.abspath(inps.outfile)

    # 2. Generate Google Earth KMZ
    write_kmz_file(data,
                   metadata=atr,
                   out_file=inps.outfile,
                   inps=inps)
    return inps.outfile


#######################################################
if __name__ == '__main__':
    main()
