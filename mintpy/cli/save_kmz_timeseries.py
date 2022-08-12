############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import os
import sys
from mintpy.utils.arg_utils import create_argument_parser


############################################################
EXAMPLE = """example:
  cd $PROJECT_NAME/mintpy/geo
  save_kmz_timeseries.py geo_timeseries_ERA5_ramp_demErr.h5
  save_kmz_timeseries.py geo_timeseries_ERA5_ramp_demErr.h5 -v -5 5 --wrap
  save_kmz_timeseries.py timeseries_ERA5_demErr.h5 --vel velocity.h5 --tcoh temporalCoherence.h5 --mask maskTempCoh.h5
"""

def create_parser(subparsers=None):
    synopsis = 'Generare Google Earth KMZ file for time-series file.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    args = parser.add_argument_group('Input files', 'File/Dataset to display')

    args.add_argument('ts_file', metavar='timeseries_file', help='Timeseries file to generate KML for')
    args.add_argument('--vel', dest='vel_file', metavar='FILE',
                      help='Velocity file, used for the color of dot')
    args.add_argument('--tcoh', dest='tcoh_file', metavar='FILE',
                      help='temporal coherence file, used for stat info')
    args.add_argument('--mask', dest='mask_file', metavar='FILE',
                      help='Mask file')
    args.add_argument('-o','--output', dest='outfile', help='Output KMZ file name.')

    opts = parser.add_argument_group('Display options', 'configurations for the display')
    opts.add_argument('--steps', type=int, nargs=3, default=[20, 5, 2],
                      help='list of steps for output pixel (default: %(default)s).\n'
                           'Set to [20, 5, 0] to skip the 3rd high-resolution level to reduce file size.')
    opts.add_argument('--level-of-details','--lods', dest='lods', type=int, nargs=4, default=[0, 1500, 4000, -1],
                      help='list of level of details to determine the visible range while browering. Default: 0, 1500, 4000, -1.\n'+
                           'Ref: https://developers.google.com/kml/documentation/kml_21tutorial')
    opts.add_argument('--vlim','-v', dest='vlim', nargs=2, metavar=('VMIN', 'VMAX'), type=float,
                      help='min/max range in cm/yr for color coding.')
    opts.add_argument('--wrap', dest='wrap', action='store_true',
                      help='re-wrap data to [VMIN, VMAX) for color coding.')
    opts.add_argument('--colormap','-c', dest='cmap_name', default='jet',
                      help='colormap used for display, i.e. jet, RdBu, hsv, jet_r, temperature, viridis,  etc.\n'
                           'More details at https://mintpy.readthedocs.io/en/latest/api/colormaps/')

    defo = parser.add_argument_group('HD for deforming areas', 'High resolution output for deforming areas')
    defo.add_argument('--cutoff', dest='cutoff', type=int, default=3,
                      help='choose points with velocity >= cutoff * MAD. Default: 3.')
    defo.add_argument('--min-percentage','--min-perc', dest='min_percentage', type=float, default=0.2,
                      help='choose boxes with >= min percentage of pixels are deforming. Default: 0.2.')

    parser.add_argument('--kk','--keep-kml','--keep-kml-file', dest='keep_kml_file', action='store_true',
                        help='Do not remove KML and data/resource files after compressing into KMZ file.')

    return parser


def cmd_line_parse(iargs=None):
    from ..utils import readfile
    from ..save_kmz_timeseries import get_aux_filename

    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check if in geo coordinates
    atr = readfile.read_attribute(inps.ts_file)
    if "Y_FIRST" not in atr.keys():
        raise ValueError("input file {} is NOT geocoded".format(inps.ts_file))

    inps = get_aux_filename(inps)
    for fname in [inps.vel_file, inps.tcoh_file, inps.mask_file]:
        if not os.path.isfile(fname):
            raise FileNotFoundError('auxliary file {} not found.'.format(fname))
    return inps


######################################################################################
def main(iargs=None):
    from ..save_kmz_timeseries import save_kml_timeseries

    inps = cmd_line_parse(iargs)
    inps.work_dir = os.path.abspath(os.path.dirname(inps.ts_file))
    inps.cbar_file = os.path.join(inps.work_dir, 'google_earth_cbar.png')
    inps.star_file = os.path.join(inps.work_dir, "star.png")
    inps.dot_file = os.path.join(inps.work_dir, "shaded_dot.png")
    inps.dygraph_file = os.path.join(inps.work_dir, "dygraph-combined.js")
    inps.kml_data_dir = os.path.join(inps.work_dir, 'kml_data')

    return save_kml_timeseries(inps)


######################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
