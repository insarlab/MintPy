############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import os
import sys
import time

from mintpy.defaults.template import get_template_content
from mintpy.utils import arg_utils


############################################################################
TEMPLATE = get_template_content('velocity')

REFERENCE = """references:
  Fattahi, H., and F. Amelung (2015), InSAR bias and uncertainty due to the systematic and stochastic
    tropospheric delay, J. Geophy. Res. Solid Earth, 120(12), 8758-8773, doi:10.1002/2015JB012419.
  Efron, B., and R. Tibshirani (1986), Bootstrap methods for standard errors, confidence intervals,
    and other measures of statistical accuracy, Statistical Science, 54-75, doi:10.1214/ss/1177013815.
"""

EXAMPLE = """example:
  timeseries2velocity.py timeseries_ERA5_demErr.h5
  timeseries2velocity.py timeseries_ERA5_demErr_ramp.h5 -t KyushuAlosDT73.txt
  timeseries2velocity.py timeseries.h5  --start-date 20080201  --end-date 20100508
  timeseries2velocity.py timeseries.h5  --ex exclude_date.txt

  # complex time functions
  timeseries2velocity.py timeseries_ERA5_demErr.h5 --poly 3 --period 1 0.5 --step 20170910
  timeseries2velocity.py timeseries_ERA5_demErr.h5 --poly 1 --exp 20170910 90
  timeseries2velocity.py timeseries_ERA5_demErr.h5 --poly 1 --log 20170910 60.4
  timeseries2velocity.py timeseries_ERA5_demErr.h5 --poly 1 --log 20170910 60.4 200 --log 20171026 200.7

  # uncertainty quantification of the estimated time functions
  timeseries2velocity.py timeseries_ERA5_demErr.h5 --uq residue
  timeseries2velocity.py timeseries_ERA5_demErr.h5 --uq covariance --ts-cov timeseriesCov.h5
  timeseries2velocity.py timeseries_ERA5_demErr.h5 --uq bootstrap
"""

DROP_DATE_TXT = """exclude_date.txt:
20040502
20060708
20090103
"""


def create_parser(subparsers=None):
    synopsis = 'Estimate velocity / time functions from time-series.'
    epilog = REFERENCE + '\n' + TEMPLATE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = arg_utils.create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    # inputs
    parser.add_argument('timeseries_file', help='Time series file for time function estimation.')
    parser.add_argument('--template', '-t', dest='template_file', help='template file with options')

    # outputs
    parser.add_argument('-o', '--output', dest='outfile', help='output file name')
    parser.add_argument('--update', dest='update_mode', action='store_true',
                        help='Enable update mode, and skip estimation if:\n'+
                             '1) output file already exists, readable '+
                             'and newer than input file\n' +
                             '2) all configuration parameters are the same.')

    # reference in time and space
    # useful for input file without reference info, e.g. ERA5.h5
    parser = arg_utils.add_reference_argument(parser, plot=False)

    # dates of interest
    date = parser.add_argument_group('Dates of interest')
    date.add_argument('-s','--start-date', dest='startDate',
                      help='start date for time function estimation')
    date.add_argument('-e','--end-date', dest='endDate',
                      help='end date for time function estimation')
    date.add_argument('--ex','--ex-date', dest='excludeDate', nargs='+', default=[],
                      help='date(s) not included in time function estimation, i.e.:\n' +
                           '--exclude 20040502 20060708 20090103\n' +
                           '--exclude exclude_date.txt\n'+DROP_DATE_TXT)

    # Uncertainty quantification
    uq = parser.add_argument_group('Uncertainty quantification (UQ)', 'Estimating the time function parameters STD')
    uq.add_argument('--uq', '--uncertainty', dest='uncertaintyQuantification', metavar='VAL',
                    default='residue', choices={'residue', 'covariance', 'bootstrap'},
                    help='Uncertainty quantification method (default: %(default)s).')
    uq.add_argument('--ts-cov','--ts-cov-file', dest='timeSeriesCovFile',
                    help='4D time-series (co)variance file for time function STD calculation')
    uq.add_argument('--bc', '--bootstrap-count', dest='bootstrapCount', type=int, default=400,
                    help='number of iterations for bootstrapping (default: %(default)s).')

    # time functions
    parser = arg_utils.add_timefunc_argument(parser)

    # residual file
    resid = parser.add_argument_group('Residual file', 'Save residual displacement time-series to HDF5 file.')
    resid.add_argument('--save-res', '--save_residual', dest='save_res', action='store_true',
                       help='Save the residual displacement time-series to HDF5 file.')
    resid.add_argument('--res-file', '--residual-file', dest='res_file', default='timeseriesResidual.h5',
                       help='Output file name for the residual time-series file (default: %(default)s).')

    # computing
    parser = arg_utils.add_memory_argument(parser)

    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    from mintpy.utils import readfile, utils as ut
    from mintpy.timeseries2velocity import read_template2inps
    
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check if input file is time series
    inps.file_type = readfile.read_attribute(inps.timeseries_file)['FILE_TYPE']
    if inps.file_type not in ['timeseries', 'giantTimeseries', 'HDFEOS']:
        raise Exception('input file is {}, NOT timeseries!'.format(inps.file_type))

    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    # --uq
    if inps.uncertaintyQuantification == 'bootstrap':
        # check 1 - bootstrap count number
        if inps.bootstrapCount <= 1:
            inps.uncertaintyQuantification = 'residue'
            print('WARNING: bootstrapCount should be larger than 1!')
            print('Change the uncertainty quantification method from bootstrap to residue, and continue.')
        # check 2 - advanced time func
        if (inps.polynomial != 1 or inps.periodic or inps.stepDate or inps.exp or inps.log):
            raise ValueError('bootstrapping support polynomial with the order of 1 ONLY!')

    elif inps.uncertaintyQuantification == 'covariance':
        if not inps.timeSeriesCovFile or not os.path.isfile(inps.timeSeriesCovFile):
            inps.uncertaintyQuantification = 'residue'
            print('WARNING: NO time series covariance file found!')
            print('Change the uncertainty quantification method from covariance to residue, and continue.')

    # --ref-lalo option
    if inps.ref_lalo:
        atr = readfile.read_attribute(inps.timeseries_file)
        coord = ut.coordinate(atr)
        ref_y, ref_x = coord.geo2radar(inps.ref_lalo[0], inps.ref_lalo[1])[:2]
        if ref_y is not None and ref_x is not None:
            inps.ref_yx = [ref_y, ref_x]
            print('input reference point in (lat, lon): ({}, {})'.format(inps.ref_lalo[0], inps.ref_lalo[1]))
            print('corresponding   point in (y, x): ({}, {})'.format(inps.ref_yx[0], inps.ref_yx[1]))

    # --output
    if not inps.outfile:
        # get suffix
        ts_file_base = os.path.splitext(os.path.basename(inps.timeseries_file))[0]
        if ts_file_base in ['timeseriesRg', 'timeseriesAz']:
            suffix = ts_file_base.split('timeseries')[-1]
        else:
            suffix = ''
        # compose default output filename
        inps.outfile = f'velocity{suffix}.h5'

    return inps


############################################################################
def main(iargs=None):
    from mintpy.timeseries2velocity import run_or_skip, run_timeseries2time_func
    inps = cmd_line_parse(iargs)
    start_time = time.time()

    # --update option
    if inps.update_mode and run_or_skip(inps) == 'skip':
        return inps.outfile

    # run
    run_timeseries2time_func(inps)

    # time info
    m, s = divmod(time.time()-start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.'.format(m, s))
    return


############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
