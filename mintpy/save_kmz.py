#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################
# recommended import
#   from mintpy import save_kmz


import os
import sys
import shutil
import numpy as np
from lxml import etree
from zipfile import ZipFile
from matplotlib import pyplot as plt, colors, colorbar, ticker

try:
    from pykml.factory import KML_ElementMaker as KML
except ImportError:
    raise ImportError('Can not import pykml!')

import mintpy
from mintpy.objects import timeseriesKeyNames
from mintpy.utils import (
    arg_utils,
    attribute as attr,
    ptime,
    readfile,
    utils as ut,
    plot as pp,
)
from mintpy import subset


############################################################
EXAMPLE = """example:
  save_kmz.py geo/geo_velocity.h5 
  save_kmz.py geo/geo_velocity.h5 -u cm --wrap --wrap-range -3 7

  save_kmz.py geo/geo_timeseries_ERA5_ramp_demErr.h5 20101120
  save_kmz.py geo/geo_timeseries_ERA5_demErr.h5 20200505_20200517

  save_kmz.py geo/geo_ifgramStack.h5 20101120_20110220
  save_kmz.py geo/geo_geometryRadar.h5 height --cbar-label Elevation

  # to generate placemarks for the file in radar coordinates, the corresponding
  # geometry file with latitude & longitude in radar coordinates are required,
  # such as provided by ISCE + MintPy workflow
  save_kmz.py velocity.h5 --sub-x 300 800 --sub-y 1000 1500 --step 1
"""


def create_parser(subparsers=None):
    synopsis = 'Generate Google Earth KMZ file (overlay / placemarks for files in geo / radar coordinates).'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = arg_utils.create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', help='file to be converted, in geo or radar coordinate.\n'
                        'Note: for files in radar-coordinate, the corresponding lookup table\n'
                        'in radar-coordinate (as provided by ISCE) is required.')
    parser.add_argument('dset', nargs='?',
                        help='date of timeseries, or date12 of interferograms to be converted')
    parser.add_argument('-m','--mask', dest='mask_file', metavar='FILE',
                        help='mask file for display')
    parser.add_argument('--zero-mask', dest='zero_mask', action='store_true',
                        help='Mask pixels with zero value.')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output file base name. Extension is fixed with .kmz')
    parser.add_argument('--kk','--keep-kml','--keep-kml-file', dest='keep_kml_file', action='store_true',
                        help='Do not remove KML and data/resource files after compressing into KMZ file.')

    # unique for point - file in radar coordinates
    parser.add_argument('-g','--geom', dest='geom_file', metavar='FILE',
                        help='geometry file with lat/lon. [required for file in radar coordinates]')
    parser.add_argument('--step', dest='step', type=int, default=5,
                        help='output one point per {step} pixels, to reduce file size (default: %(default)s).\n'
                             'For file in radar-coordinate ONLY.')

    # Data
    parser.add_argument('-v','--vlim', dest='vlim', nargs=2, metavar=('MIN', 'MAX'), type=float,
                        help='Y/value limits for plotting.')
    parser.add_argument('-u', dest='disp_unit', metavar='UNIT', help='unit for display.')
    parser.add_argument('-c', '--cm', '--colormap', dest='cmap_name', default='jet',
                        help='Colormap for plotting (default: %(default)s), such as jet, RdBu, etc.\n'
                             'More details at https://mintpy.readthedocs.io/en/latest/api/colormaps/')
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
    fig.add_argument('--cbar-loc', dest='cbar_loc', default='lower left',
                     choices=['lower left','lower right','upper left', 'upper right'],
                     help='Location of colorbar in the screen. Default: lower left.')
    fig.add_argument('--cbar-label', dest='cbar_label', metavar='LABEL', default='Mean LOS velocity',
                     help='Colorbar label. Default: Mean LOS velocity')
    fig.add_argument('--cbar-bin-num', dest='cbar_bin_num', metavar='NUM', type=int,
                     help='Colorbar bin number (default: %(default)s).')

    # Reference Pixel
    ref = parser.add_argument_group('Reference Pixel')
    ref.add_argument('--noreference', dest='disp_ref_pixel', action='store_false',
                     help='do not show reference point')
    ref.add_argument('--ref-color', dest='ref_marker_color', metavar='COLOR', default='k',
                     help='marker color of reference point')
    ref.add_argument('--ref-size', dest='ref_marker_size', metavar='NUM', type=int, default=5,
                     help='marker size of reference point (default: %(default)s).')
    ref.add_argument('--ref-marker', dest='ref_marker', metavar='SYMBOL', default='s',
                     help='marker symbol of reference point')

    # subset
    parser = arg_utils.add_subset_argument(parser)

    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    inps.work_dir = os.path.abspath(os.path.dirname(inps.file))

    atr = readfile.read_attribute(inps.file)

    # Check 1: file in radar coord
    if 'Y_FIRST' not in atr.keys():
        geom_ds_list = ['latitude', 'longitude']
        if not inps.geom_file:
            inps.geom_file = ut.get_geometry_file(
                geom_ds_list,
                work_dir=inps.work_dir,
                coord='radar',
            )
        if not inps.geom_file or not os.path.isfile(inps.geom_file):
            raise FileNotFoundError(f'No geometry file with {geom_ds_list} in radar coord found!')

    # Check 2: dset is required for multi_dataset/group files
    if not inps.dset and atr['FILE_TYPE'] in ['ifgramStack'] + timeseriesKeyNames:
        raise Exception("No date/date12 input.\nIt's required for {} file".format(atr['FILE_TYPE']))

    # Backend setting
    plt.switch_backend('Agg')

    return inps


############################################################
def plot_colorbar(out_file, vmin, vmax, unit='cm/year', cmap='jet', figsize=(0.18, 3.6),
                  nbins=None, label='Mean LOS velocity'):
    fig, cax = plt.subplots(figsize=figsize)
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cbar = colorbar.ColorbarBase(cax, cmap=plt.get_cmap(cmap), norm=norm, orientation='vertical')
    cbar.set_label('{} [{}]'.format(label, unit), fontsize=12)

    # update ticks
    if nbins:
        cbar.locator = ticker.MaxNLocator(nbins=nbins)
        cbar.update_ticks()

    cbar.ax.tick_params(which='both', labelsize=12)
    fig.patch.set_facecolor('white')
    fig.patch.set_alpha(0.7)
    print('writing', out_file)
    fig.savefig(out_file, bbox_inches='tight', facecolor=fig.get_facecolor(), dpi=300)
    return out_file


def generate_cbar_element(cbar_file, cmap, vmin, vmax, unit='cm/year', loc='lower left',
                          nbins=None, label='Mean LOS velocity'):
    """Generate colorbar as an screen overlay object.

    Parameters: cbar_file - str, colorbar image file path
                cmap      - matplotlib.colors.Colormap instance
                vmin/vmax - float, min/max value to display
                unit      - str, display unit
                loc       - str, location of colorbar on the screen.
                            lower-left, lower-right, upper-left, upper-right
    Returns:    cbar_overlay - KML.ScreenOverlay object
    """
    # plot colobar and save as an image
    cbar_file = plot_colorbar(
        out_file=cbar_file,
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


def get_hex_color(v, colormap, norm):
    """Get color name in hex format.
    Parameters: v        : float, number of interest
                colormap : matplotlib.colors.Colormap instance
                norm     : matplotlib.colors.Normalize instance
    Returns:    c_hex    : color name in hex format
    """
    # get rgba color components for point velocity
    rgba = colormap(norm(v))
    # rgba to hex
    c_hex = colors.to_hex([rgba[3], rgba[2], rgba[1], rgba[0]], keep_alpha=True)[1:]
    return c_hex


def create_placemark_element(lat, lon, row, col, val, icon_file, inps):
    """Create an KMZ Placemark element.
    Parameters: lat/lon   - float, latitude / longitude in degrees
                row/col   - int, row / column number
                val       - float, value
                icon_file - str, path of the icon file
                inps      - Namespace object
    Returns:    placemark - KMZ.Placemark() object
    """

    placemark = KML.Placemark(
        # style
        KML.Style(
            KML.IconStyle(
                    KML.color(get_hex_color(val, inps.colormap, inps.norm)),
                    KML.scale(0.5),
                    KML.Icon(KML.href(icon_file),
                ),
            ),
        ),

        # extended data
        KML.ExtendedData(
            KML.Data(
                KML.value(f"{lat:.6f}˚"),
                name="Latitude",
            ),
            KML.Data(
                KML.value(f"{lon:.6f}˚"),
                name="Longitude",
            ),
            KML.Data(
                KML.value(f"{row:.0f}"),
                name="Row",
            ),
            KML.Data(
                KML.value(f"{col:.0f}"),
                name="Column",
            ),
            KML.Data(
                KML.value(f"{val:.2f} {inps.disp_unit}"),
                name="Value",
            ),
        ),

        # point coord
        KML.Point(
            KML.coordinates(f"{lon},{lat}"),
        ),
    )

    return placemark


def write_kmz_file(out_file_base, kml_doc, data_files=None, res_files=None, keep_kml_file=False):
    """Write KML and KMZ files.
    Parameters: out_file_base - str, output file name without extension
                kml_doc       - KML.Document() object
                data_files    - list of str, rel path of data files
                res_files     - list of str, rel path of resource files
                keep_kml_file - bool, do not remove KML files after zipping.
    Returns:    kmz_file      - str, zipped KMZ file.
    """
    # default values
    data_files = [] if data_files is None else data_files
    res_files  = [] if res_files  is None else res_files

    work_dir = os.path.dirname(out_file_base)
    kml_file = '{}.kml'.format(out_file_base)
    kmz_file = '{}.kmz'.format(out_file_base)

    # 1. Write KML file
    kml = KML.kml()
    kml.append(kml_doc)

    print('writing '+kml_file)
    with open(kml_file, 'w') as f:
        f.write(etree.tostring(kml, pretty_print=True).decode('utf8'))

    # 2. Copy resource files
    if res_files:
        res_dir = os.path.join(os.path.dirname(mintpy.__file__), "data")
        for fname in res_files:
            src_file = os.path.join(res_dir, os.path.basename(fname))
            shutil.copy2(src_file, work_dir)
            print("copy {} to the local directory".format(src_file))

    # 3. Generate KMZ file, by
    # 1) go to the directory of kmz file
    run_dir = os.path.abspath(os.getcwd())
    os.chdir(work_dir)

    # 2) zip all data files
    with ZipFile(kmz_file, 'w') as fz:
        for fname in [kml_file] + data_files + res_files:
            fz.write(os.path.relpath(fname))
            if not keep_kml_file:
                os.remove(fname)
                print('remove {}'.format(fname))

    # 3) go back to the running directory
    os.chdir(run_dir)
    print('merged all files to {}'.format(kmz_file))

    return kmz_file


############################################################
def write_kmz_overlay(data, meta, out_file, inps):
    """Generate Google Earth Overlay KMZ file for data in GEO coordinates.
    Parameters: data     - 2D np.array in int/float, data matrix to write
                meta     - dict, containing the following attributes:
                           WIDTH/LENGTH      : required, file size
                           X/Y_FIRST/STEP    : required, for lat/lon spatial converage
                           REF_X/Y           : optional, column/row number of reference pixel
                out_file - string, output file name
                inps     - Namespace, optional, input options for display
    Returns:    kmz_file - string, output KMZ filename
    """

    south, north, west, east = ut.four_corners(meta)

    # 1. Make PNG file - Data
    print('plotting data ...')

    # Figure size
    if not inps.fig_size:
        inps.fig_size = pp.auto_figure_size(ds_shape=[north-south, east-west], scale=2.0)
    fig = plt.figure(figsize=inps.fig_size, frameon=False)
    ax = fig.add_axes([0., 0., 1., 1.])
    ax.set_axis_off()

    # Plot - data matrix
    ax.imshow(data, vmin=inps.vlim[0], vmax=inps.vlim[1], cmap=inps.colormap,
              aspect='auto', interpolation='nearest')

    # Plot - reference pixel
    rx = meta.get('REF_X', None)
    ry = meta.get('REF_Y', None)
    if inps.disp_ref_pixel and rx is not None and ry is not None:
        ax.plot(int(rx), int(ry), inps.ref_marker,
                color=inps.ref_marker_color,
                ms=inps.ref_marker_size)
        print('show reference point')
    else:
        print('no plot for reference point.')

    width = int(meta['WIDTH'])
    length = int(meta['LENGTH'])
    ax.set_xlim([0, width])
    ax.set_ylim([length, 0])

    out_file_base = os.path.splitext(out_file)[0]
    data_png_file = out_file_base + '.png'
    print('writing {} with dpi={}'.format(data_png_file, inps.fig_dpi))
    plt.savefig(data_png_file, pad_inches=0.0, transparent=True, dpi=inps.fig_dpi)

    # 2. Generate KML file
    kml_doc = KML.Document()

    # Add data png file
    img_name = os.path.splitext(os.path.basename(data_png_file))[0]
    img_overlay = KML.GroundOverlay(
        KML.name(img_name),
        KML.Icon(
            KML.href(os.path.basename(data_png_file))
        ),
        KML.altitudeMode('clampToGround'),
        KML.LatLonBox(
            KML.north(str(north)),
            KML.east(str(east)),
            KML.south(str(south)),
            KML.west(str(west)),
        ),
    )
    kml_doc.append(img_overlay)

    # Add colorbar png file
    cbar_file = '{}_cbar.png'.format(out_file_base)
    cbar_overlay = generate_cbar_element(
        cbar_file,
        cmap=inps.colormap,
        vmin=inps.vlim[0],
        vmax=inps.vlim[1],
        unit=inps.disp_unit,
        loc=inps.cbar_loc,
        nbins=inps.cbar_bin_num,
        label=inps.cbar_label)
    kml_doc.append(cbar_overlay)

    # Write KML file
    kmz_file = write_kmz_file(
        out_file_base,
        kml_doc,
        data_files=[data_png_file, cbar_file],
        keep_kml_file=inps.keep_kml_file)

    return kmz_file


def write_kmz_placemark(data, meta, out_file, geom_file, inps):
    """Generate Google Earth Placemark KMZ file for data in RADAR coordinates.
    Parameters: data      - 2D np.array in int/float, data matrix to write
                meta      - dict, containing the following attributes:
                            WIDTH/LENGTH      : required, file size
                            X/Y_FIRST/STEP    : required, for lat/lon spatial converage
                            REF_X/Y           : optional, column/row number of reference pixel
                geom_file - str, path of the geometry file with latitude/longitude datasets
                out_file  - string, output file name
                inps      - Namespace, optional, input options for display
    Returns:    kmz_file  - string, output KMZ filename
    """

    out_file_base = os.path.splitext(out_file)[0]
    dot_file = 'shaded_dot.png'
    star_file = 'star.png'

    # read latitude / longitude
    lats = readfile.read(geom_file, datasetName='latitude',  box=inps.pix_box)[0]
    lons = readfile.read(geom_file, datasetName='longitude', box=inps.pix_box)[0]

    ## Generate KML file
    kml_doc = KML.Document()

    # 1. colorbar png file
    print('plot and add colorbar as a ScreenOverlay element')
    cbar_file = '{}_cbar.png'.format(out_file_base)
    cbar_overlay = generate_cbar_element(
        cbar_file,
        cmap=inps.colormap,
        vmin=inps.vlim[0],
        vmax=inps.vlim[1],
        unit=inps.disp_unit,
        loc=inps.cbar_loc,
        nbins=inps.cbar_bin_num,
        label=inps.cbar_label)
    kml_doc.append(cbar_overlay)

    # 2. reference point
    xmin = int(meta.get('SUBSET_XMIN', 0))
    ymin = int(meta.get('SUBSET_YMIN', 0))

    if 'REF_Y' in meta.keys():
        print('add reference point as a star icon')
        ry, rx = int(meta['REF_Y']), int(meta['REF_X'])
        rlat = lats[ry, rx]
        rlon = lons[ry, rx]
        ref_point = create_placemark_element(
            lat=rlat,
            lon=rlon,
            row=ry + ymin,
            col=rx + xmin,
            val=0.0,
            icon_file=star_file,
            inps=inps)
        ref_point.name = 'ReferencePoint'
        ref_point.Style.IconStyle.scale = 1.0
        kml_doc.append(ref_point)

        # do not plot reference point as data again
        data[ry, rx] = np.nan

    # 3. data folder for all points
    data_folder = KML.Folder(KML.name("Data"))

    print('generating point element with step size of {} pixels'.format(inps.step))
    length, width = data.shape
    prog_bar = ptime.progressBar(maxValue=length)
    for y in range(0, length, inps.step):
        for x in range(0, width, inps.step):
            value = data[y, x]
            if not np.isnan(value):
                lat = lats[y, x]
                lon = lons[y, x]

                # create KML icon element
                placemark = create_placemark_element(
                    lat=lat,
                    lon=lon,
                    row=y + ymin,
                    col=x + xmin,
                    val=value,
                    icon_file=dot_file,
                    inps=inps)
                data_folder.append(placemark)

        prog_bar.update(y+1, every=1, suffix=f'row={y+1}/{length}')
    prog_bar.close()
    kml_doc.append(data_folder)

    # Write KML file
    kmz_file = write_kmz_file(
        out_file_base,
        kml_doc,
        data_files=[cbar_file],
        res_files=[dot_file, star_file],
        keep_kml_file=inps.keep_kml_file)

    return kmz_file


############################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    # 1. Read metadata and data
    k = readfile.read_attribute(inps.file)['FILE_TYPE']
    if k == 'timeseries' and inps.dset and '_' in inps.dset:
        inps.ref_date, inps.dset = inps.dset.split('_')
    else:
        inps.ref_date = None
    atr = readfile.read_attribute(inps.file, datasetName=inps.dset)

    # pix_box
    inps.pix_box = subset.subset_input_dict2box(vars(inps), atr)[0]
    inps.pix_box = ut.coordinate(atr).check_box_within_data_coverage(inps.pix_box)
    data_box = (0, 0, int(atr['WIDTH']), int(atr['LENGTH']))
    print('data   coverage in y/x: {}'.format(data_box))
    print('subset coverage in y/x: {}'.format(inps.pix_box))
    atr = attr.update_attribute4subset(atr, inps.pix_box)

    # read data
    data = readfile.read(inps.file, datasetName=inps.dset, box=inps.pix_box)[0]
    if k == 'timeseries' and inps.ref_date:
        data -= readfile.read(inps.file, datasetName=inps.ref_date, box=inps.pix_box)[0]

    # mask
    mask = pp.read_mask(inps.file, mask_file=inps.mask_file,
                        datasetName=inps.dset, box=inps.pix_box)[0]
    if mask is not None:
        print('masking out pixels with zero value in file: {}'.format(inps.mask_file))
        data[mask == 0] = np.nan
    if inps.zero_mask:
        print('masking out pixels with zero value')
        data[data == 0] = np.nan
    del mask

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


    # 2. Generate Google Earth KMZ
    # 2.1 Common settings
    # disp min/max and colormap
    cmap_lut = 256
    if not inps.vlim:
        cmap_lut, inps.vlim = pp.auto_adjust_colormap_lut_and_disp_limit(data)[:2]
    inps.cmap_name = pp.auto_colormap_name(atr, inps.cmap_name)
    inps.colormap = pp.ColormapExt(inps.cmap_name, cmap_lut).colormap
    inps.norm = colors.Normalize(vmin=inps.vlim[0], vmax=inps.vlim[1])

    # Output filename
    inps.fig_title = pp.auto_figure_title(inps.file, datasetNames=inps.dset, inps_dict=vars(inps))
    if not inps.outfile:
        inps.outfile = '{}.kmz'.format(inps.fig_title)
    inps.outfile = os.path.abspath(inps.outfile)

    # 2.2 Write KMZ file
    if 'Y_FIRST' in atr.keys():
        # create ground overlay KML for file in geo-coord
        write_kmz_overlay(
            data,
            meta=atr,
            out_file=inps.outfile,
            inps=inps)

    else:
        # create placemark KML for file in radar-coord
        write_kmz_placemark(
            data,
            meta=atr,
            out_file=inps.outfile,
            geom_file=inps.geom_file,
            inps=inps)

    return


#######################################################
if __name__ == '__main__':
    main(sys.argv[1:])
