#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Heresh Fattahi, Aug 2022      #
############################################################


import argparse
import os
import sys

from mintpy.defaults.template import get_template_content
from mintpy.utils import arg_utils

############################################################################
TEMPLATE = get_template_content('correct_topography')

REFERENCE = """reference:
  Fattahi, H., and F. Amelung (2013), DEM Error Correction in InSAR Time Series,
    IEEE Trans. Geosci. Remote Sens., 51(7), 4249-4259, doi:10.1109/TGRS.2012.2227761.
"""

EXAMPLE = """example:
  # correct DEM error with pixel-wise geometry parameters [slow]
  dem_error.py  timeseries_ERA5_ramp.h5 -g inputs/geometryRadar.h5 -t smallbaselineApp.cfg

  # correct DEM error with mean geometry parameters [fast]
  dem_error.py  timeseries_ERA5_ramp.h5 -t smallbaselineApp.cfg

  # get updated/corrected DEM
  save_roipac.py inputs/geometryGeo.h5 -o dem.h5   #for dataset in geo coordinates
  mask.py demErr.h5 -m maskTempCoh.h5 -o demErr_msk.h5
  add.py demErr_msk.h5 dem.h5 -o demNew.h5
"""

def create_parser(subparsers=None):
    synopsis = 'DEM Error (Topographic Residual) Correction'
    epilog = REFERENCE + '\n' + TEMPLATE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = arg_utils.create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('ts_file', help='Time-series HDF5 file to be corrected.')
    parser.add_argument('-g', '--geometry', dest='geom_file',
                        help='geometry file including datasets:\n'+
                             'incidence angle\n'+
                             'slant range distance\n' +
                             'and/or 3D perpendicular baseline')
    parser.add_argument('-o', '--outfile', dest='ts_cor_file',
                        help='Output file name for corrected time-series (default: add suffix of "_demErr")')
    parser.add_argument('--dem-err-file','--dem-error-file', dest='dem_err_file', default='demErr.h5',
                        help='Output file name for the estimated DEM error (default: %(default)s).')

    defo_model = parser.add_argument_group('temporal deformation model')
    defo_model.add_argument('-t', '--template', dest='template_file',
                            help='template file with the options')
    defo_model.add_argument('--ex', '--exclude', dest='excludeDate', nargs='*', default=[],
                            help='Exclude date(s) for DEM error estimation.\n' +
                                 'All dates will be corrected for DEM residual phase still.')
    defo_model.add_argument('-p', '--poly-order', dest='polyOrder', type=int, default=2,
                            help='polynomial order number of temporal deformation model (default: %(default)s).')
    defo_model.add_argument('-s', '--step-date', dest='stepFuncDate', nargs='*', default=[],
                            help='Date of step jump for temporal deformation model (default: %(default)s).'+
                                 ' i.e. date of earthquake/volcanic eruption')
    defo_model.add_argument('--periodic', '--period', '--peri', dest='periodic', type=float, nargs='+', default=[],
                            help='periodic functinos of temporal deformation model (default: %(default)s).')

    parser.add_argument('--phase-velocity', dest='phaseVelocity', action='store_true',
                        help='Use phase velocity instead of phase for inversion constrain.')
    parser.add_argument('--update', dest='update_mode', action='store_true',
                        help='Enable update mode, and skip inversion if:\n'+
                             '1) output time-series file already exists, readable '+
                             'and newer than input interferograms file\n' +
                             '2) all configuration parameters are the same.')
    # computing
    parser = arg_utils.add_memory_argument(parser)
    parser = arg_utils.add_parallel_argument(parser)

    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    from mintpy.objects import cluster

    # check
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    # check: --cluster and --num-worker option
    inps.numWorker = str(cluster.DaskCluster.format_num_worker(inps.cluster, inps.numWorker))
    if inps.cluster and inps.numWorker == '1':
        print('WARNING: number of workers is 1, turn OFF parallel processing and continue')
        inps.cluster = None

    # check: --ex option (ignore non-existed exclude_date.txt)
    if inps.excludeDate == 'exclude_date.txt' and not os.path.isfile(inps.excludeDate):
        inps.excludeDate = []

    # check: --poly-order option
    if inps.polyOrder < 1:
        raise argparse.ArgumentTypeError("Minimum polynomial order is 1")

    # check: --output / --dem-error-file option (must be HDF5 file)
    for fname in [inps.ts_cor_file, inps.dem_err_file]:
        if fname:
            fext = os.path.splitext(fname)[1]
            if fext not in ['.h5','.he5']:
                msg = f'--output / --dem-err-file option ({fext}) supports HDF5 file only!'
                raise ValueError(msg)

    # default: --output option
    if not inps.ts_cor_file:
        fbase = os.path.splitext(inps.ts_file)[0]
        inps.ts_cor_file = f'{fbase}_demErr.h5'

    return inps


def read_template2inps(template_file, inps):
    """Read input template file into inps.excludeDate"""
    print('read options from template file:', os.path.basename(template_file))

    from mintpy.dem_error import key_prefix
    from mintpy.utils import ptime, readfile, utils1 as ut

    iDict = vars(inps)
    template = readfile.read_template(template_file, skip_chars=['[', ']'])
    template = ut.check_template_auto_value(template)

    # Read template option
    key_list = [i for i in list(iDict.keys()) if key_prefix+i in template.keys()]
    for key in key_list:
        value = template[key_prefix+key]
        if key in ['phaseVelocity']:
            iDict[key] = value
        elif value:
            if key in ['polyOrder']:
                iDict[key] = int(value)
            elif key in ['excludeDate','stepFuncDate']:
                iDict[key] = ptime.yyyymmdd(value.split(','))

    # computing configurations
    dask_key_prefix = 'mintpy.compute.'
    key_list = [i for i in list(iDict.keys()) if dask_key_prefix+i in template.keys()]
    for key in key_list:
        value = template[dask_key_prefix+key]
        if key in ['cluster', 'config']:
            iDict[key] = value
        elif value:
            if key in ['numWorker']:
                iDict[key] = str(value)
            elif key in ['maxMemory']:
                iDict[key] = float(value)

    return inps


############################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.dem_error import correct_dem_error, run_or_skip

    # run or skip
    if inps.update_mode and run_or_skip(inps) == 'skip':
        return

    # run
    correct_dem_error(inps)


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
