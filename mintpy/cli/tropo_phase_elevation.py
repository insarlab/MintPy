############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import os
import sys
import argparse
from mintpy.utils.arg_utils import create_argument_parser


############################################################################
REFERENCE = """reference:
  Doin, M. P., C. Lasserre, G. Peltzer, O. Cavalie, and C. Doubre (2009), Corrections of stratified 
  tropospheric delays in SAR interferometry: Validation with global atmospheric models, J App. Geophy.,
  69(1), 35-50, doi:http://dx.doi.org/10.1016/j.jappgeo.2009.03.010.
"""

EXAMPLE = """example:
  tropo_phase_elevation.py  timeseries_demErr.h5      -g inputs/geometryRadar.h5  -m maskTempCoh.h5    
  tropo_phase_elevation.py  geo_timeseries_demErr.h5  -g geo_geometryRadar.h5     -m geo_maskTempCoh.h5
"""

def create_parser(subparsers=None):
    synopsis = 'Correct Topo-correlated Stratified tropospheric delay'
    epilog = REFERENCE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('timeseries_file', help='time-series file to be corrected')
    parser.add_argument('-g', '--geometry', dest='geom_file', required=True,
                        help='DEM file used for correlation calculation.')
    parser.add_argument('-m', '--mask', dest='mask_file', required=True,
                        help='mask file for pixels used for correlation calculation')

    parser.add_argument('-t', '--threshold', type=float, default=0.,
                        help='correlation threshold to apply phase correction.\n'
                             'if not set, all dates will be corrected.')
    parser.add_argument('-l', '--looks', dest='num_multilook', type=int, default=8,
                        help='number of looks applied to data for empirical estimation (default: %(default)s).')

    parser.add_argument('--poly-order', '-p', dest='poly_order', type=int, default=1, choices=[1, 2, 3],
                        help='polynomial order of phase-height correlation (default: %(default)s).')
    parser.add_argument('-o', '--outfile', help='output corrected timeseries file name')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    if inps.threshold and (not 0.0 <= inps.threshold <= 1.0):
        raise argparse.ArgumentTypeError('%r not in range [0.0, 1.0]' % inps.threshold)
    return inps


############################################################################
def main(iargs=None):
    from ..objects import timeseries
    from ..utils import writefile
    from ..tropo_phase_elevation import read_topographic_data, estimate_phase_elevation_ratio, estimate_tropospheric_delay

    inps = cmd_line_parse(iargs)

    # read timeseries data
    obj = timeseries(inps.timeseries_file)
    obj.open()
    ts_data = obj.read()
    inps.date_list = list(obj.dateList)

    # read topographic data (DEM)
    dem = read_topographic_data(inps.geom_file, obj.metadata)

    # estimate phase/elevation ratio parameters
    X = estimate_phase_elevation_ratio(dem, ts_data, inps)

    # correct trop delay in timeseries
    trop_data = estimate_tropospheric_delay(dem, X, obj.metadata)
    mask = ts_data == 0.
    ts_data -= trop_data
    ts_data[mask] = 0.

    # write time-series file
    meta = dict(obj.metadata)
    meta['mintpy.troposphericDelay.polyOrder'] = str(inps.poly_order)
    if not inps.outfile:
        inps.outfile = '{}_tropHgt.h5'.format(os.path.splitext(inps.timeseries_file)[0])
    writefile.write(ts_data, out_file=inps.outfile, metadata=meta, ref_file=inps.timeseries_file)


############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
