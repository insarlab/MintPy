############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import os
import sys
import time

from mintpy.defaults.template import get_template_content
from mintpy.utils.arg_utils import create_argument_parser


###############################################################
TEMPLATE = get_template_content('correct_SET')

REFERENCE = """reference:
  Milbert, D. (2018), "solid: Solid Earth Tide", [Online]. Available: http://geodesyworld.github.io/
    SOFTS/solid.htm. Accessd on: 2020-09-06.
  Yunjun, Z., Fattahi, H., Pi, X., Rosen, P., Simons, M., Agram, P., & Aoki, Y. (2022). Range
    Geolocation Accuracy of C-/L-band SAR and its Implications for Operational Stack Coregistration.
    IEEE Trans. Geosci. Remote Sens., 60, doi:10.1109/TGRS.2022.3168509.
"""

EXAMPLE = """example:
  solid_earth_tides.py timeseries.h5 -g inputs/geometryRadar.h5
  solid_earth_tides.py timeseries.h5 -g inputs/geometryGeo.h5
  solid_earth_tides.py geo/geo_timeseries_ERA5_demErr.h5 -g geo/geo_geometryRadar.h5
"""

def create_parser(subparsers=None):
    synopsis = 'Solid Earth tides (SET) correction via PySolid'
    epilog = REFERENCE + '\n' + TEMPLATE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('dis_file', help='timeseries HDF5 file, i.e. timeseries.h5')
    parser.add_argument('-g','--geomtry', dest='geom_file', type=str, required=True,
                        help='geometry file including incidence/azimuthAngle.')
    parser.add_argument('--date-wise-acq-time', dest='date_wise_acq_time', action='store_true',
                        help='Use the exact date-wise acquisition time instead of the common one for tides calculation.\n' +
                             'For ISCE-2/topsStack products only, and requires ../reference and ../secondarys folder.\n' +
                             'There is <1 min difference btw. S1A/B -> Negligible impact for InSAR.')

    parser.add_argument('--verbose', dest='verbose', action='store_true', help='Verbose message.')
    parser.add_argument('--update', dest='update_mode', action='store_true', help='Enable update mode.')

    # output
    parser.add_argument('--set-file', dest='set_file', help='line-of-sight solid earth tide file name')
    parser.add_argument('-o', dest='cor_dis_file', help='Output file name for the corrected timeseries.')

    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    from mintpy.utils import readfile

    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check input files - processors & coordinates
    atr1 = readfile.read_attribute(inps.dis_file)
    atr2 = readfile.read_attribute(inps.geom_file)
    coord1 = 'geo' if 'Y_FIRST' in atr1.keys() else 'radar'
    coord2 = 'geo' if 'Y_FIRST' in atr2.keys() else 'radar'
    proc = atr1.get('PROCESSOR', 'isce')

    if coord1 == 'radar' and proc in ['gamma', 'roipac']:
        msg = 'Radar-coded file from {} is NOT supported!'.format(proc)
        msg += '\n    Try to geocode the time-series and geometry files and re-run with them instead.'
        raise ValueError(msg)

    if coord1 != coord2:
        n = max(len(os.path.basename(i)) for i in [inps.dis_file, inps.geom_file])
        msg = 'Input time-series and geometry file are NOT in the same coordinate!'
        msg += '\n    file {f:<{n}} coordinate: {c}'.format(f=os.path.basename(inps.dis_file),  n=n, c=coord1)
        msg += '\n    file {f:<{n}} coordinate: {c}'.format(f=os.path.basename(inps.geom_file), n=n, c=coord2)
        raise ValueError(msg)

    # default SET filename
    if not inps.set_file:
        geom_dir = os.path.dirname(inps.geom_file)
        inps.set_file = os.path.join(geom_dir, 'SET.h5')

    # default corrected time-series filename
    if not inps.cor_dis_file:
        dis_dir = os.path.dirname(inps.dis_file)
        fbase, fext = os.path.splitext(os.path.basename(inps.dis_file))
        inps.cor_dis_file = os.path.join(dis_dir, '{}_SET{}'.format(fbase, fext))

    return inps


###############################################################
def main(iargs=None):
    from mintpy.solid_earth_tides import calc_solid_earth_tides_timeseries, correct_timeseries

    inps = cmd_line_parse(iargs)
    start_time = time.time()

    # calc SET
    calc_solid_earth_tides_timeseries(
        ts_file=inps.dis_file,
        geom_file=inps.geom_file,
        set_file=inps.set_file,
        date_wise_acq_time=inps.date_wise_acq_time,
        update_mode=inps.update_mode,
        verbose=inps.verbose)

    # correct SET
    correct_timeseries(
        dis_file=inps.dis_file,
        set_file=inps.set_file,
        cor_dis_file=inps.cor_dis_file)

    m, s = divmod(time.time() - start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.\n'.format(m, s))


###############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
