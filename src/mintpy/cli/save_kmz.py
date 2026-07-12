#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import math
import os
import sys

from mintpy.utils import arg_utils

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

    # data
    parser.add_argument('-v','--vlim', dest='vlim', nargs=2, metavar=('MIN', 'MAX'), type=float,
                        help='Y/value limits in the unit of {-u/--unit} for plotting.')
    parser.add_argument('-u','--unit', dest='disp_unit', metavar='UNIT', default='cm/year',
                        help='unit for display (default: %(default)s).')
    parser.add_argument('-c', '--cm', '--colormap', dest='cmap_name', default='jet',
                        help='Colormap for plotting (default: %(default)s), such as jet, RdBu, etc.\n'
                             'More details at https://mintpy.readthedocs.io/en/latest/api/colormaps/')
    parser.add_argument('--wrap', action='store_true',
                        help='re-wrap data to display data in fringes.')
    parser.add_argument('--wrap-range', dest='wrap_range', type=float, nargs=2,
                        default=[-1.*math.pi, math.pi], metavar=('MIN', 'MAX'),
                        help='range of one cycle after wrapping, default: [-pi, pi]')

    # figure
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

    # reference pixel
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
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    from mintpy.objects import TIMESERIES_KEY_NAMES
    from mintpy.utils import readfile, utils as ut

    # check
    inps.work_dir = os.path.abspath(os.path.dirname(inps.file))
    atr = readfile.read_attribute(inps.file)

    # default + check: geom_file for file in radar coord
    if 'Y_FIRST' not in atr.keys():
        geom_ds_list = ['latitude', 'longitude']
        # default geom_file
        if not inps.geom_file:
            inps.geom_file = ut.get_geometry_file(
                geom_ds_list,
                work_dir=inps.work_dir,
                coord='radar')
        # check existence
        if not inps.geom_file or not os.path.isfile(inps.geom_file):
            msg = f'No geometry file with {geom_ds_list} in radar coord found!'
            raise FileNotFoundError(msg)

    # check: dset option (required for timeseries and ifgramStack files)
    ftype = atr['FILE_TYPE']
    if not inps.dset and ftype in TIMESERIES_KEY_NAMES + ['ifgramStack']:
        raise Exception(f'No date/date12 specified for {ftype} file!')

    return inps


############################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.save_kmz import save_kmz

    # run
    save_kmz(inps)


#######################################################
if __name__ == '__main__':
    main(sys.argv[1:])
