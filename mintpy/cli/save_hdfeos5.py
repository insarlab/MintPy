############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import os
import sys

from mintpy.defaults.template import get_template_content
from mintpy.utils.arg_utils import create_argument_parser


################################################################
TEMPALTE = TEMPLATE = get_template_content('hdfeos5')

EXAMPLE = """example:
  save_hdfeos5.py geo/geo_timeseries_ERA5_ramp_demErr.h5
  save_hdfeos5.py timeseries_ERA5_ramp_demErr.h5 --tc temporalCoherence.h5 --asc avgSpatialCoh.h5 -m maskTempCoh.h5 -g inputs/geometryGeo.h5
  save_hdfeos5.py timeseries_ERA5_ramp_demErr.h5 --tc temporalCoherence.h5 --asc avgSpatialCoh.h5 -m maskTempCoh.h5 -g inputs/geometryRadar.h5
"""

NOTE = """
  https://earthdata.nasa.gov/esdis/eso/standards-and-references/hdf-eos5
  https://mintpy.readthedocs.io/en/latest/hdfeos5/
"""


def create_parser(subparsers=None):
    synopsis = 'Convert MintPy timeseries product into HDF-EOS5 format'
    epilog = TEMPALTE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis+NOTE, epilog=epilog, subparsers=subparsers)

    parser.add_argument('ts_file', default='timeseries.h5', help='Time-series file')
    parser.add_argument('-t', '--template', dest='template_file',
                        help='Template file for 1) arguments/options and 2) missing metadata')

    parser.add_argument('--tc','--temp-coh', dest='tcoh_file',
                        help='Coherence/correlation file, i.e. temporalCoherence.h5')
    parser.add_argument('--asc','--avg-spatial-coh', dest='scoh_file',
                        help='Average spatial coherence file, i.e. avgSpatialCoh.h5')
    parser.add_argument('-m', '--mask', dest='mask_file', help='Mask file')
    parser.add_argument('-g', '--geometry', dest='geom_file', help='geometry file')
    parser.add_argument('--suffix', dest='suffix', help='suffix to be appended to file name (e.g. PS).')

    parser.add_argument('--update', action='store_true',
                        help='Enable update mode, a.k.a. put XXXXXXXX as endDate in filename if endDate < 1 year')
    parser.add_argument('--subset', action='store_true',
                        help='Enable subset mode, a.k.a. put suffix _N31700_N32100_E130500_E131100')
    return parser


def cmd_line_parse(iargs=None):
    from mintpy.utils import readfile

    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # default filenames
    ts_dir = os.path.dirname(inps.ts_file)
    meta = readfile.read_attribute(inps.ts_file)
    if os.path.basename(inps.ts_file).startswith('geo_'):
        tcoh_file = os.path.join(ts_dir, 'geo_temporalCoherence.h5')
        scoh_file = os.path.join(ts_dir, 'geo_avgSpatialCoh.h5')
        mask_file = os.path.join(ts_dir, 'geo_maskTempCoh.h5')
        geom_file = os.path.join(ts_dir, 'geo_geometryRadar.h5')
    else:
        tcoh_file = os.path.join(ts_dir, 'temporalCoherence.h5')
        scoh_file = os.path.join(ts_dir, 'avgSpatialCoh.h5')
        mask_file = os.path.join(ts_dir, 'maskTempCoh.h5')

        if 'Y_FIRST' in meta.keys():
            geom_file = os.path.join(ts_dir, 'inputs/geometryGeo.h5')
        else:
            geom_file = os.path.join(ts_dir, 'inputs/geometryRadar.h5')


    if not inps.tcoh_file:  inps.tcoh_file = tcoh_file
    if not inps.scoh_file:  inps.scoh_file = scoh_file
    if not inps.mask_file:  inps.mask_file = mask_file
    if not inps.geom_file:  inps.geom_file = geom_file

    # check file existence
    for fname in [inps.ts_file, inps.tcoh_file, inps.scoh_file, inps.mask_file, inps.geom_file]:
        if not os.path.isfile(fname):
            raise FileNotFoundError(fname)

    # --subset mode
    if inps.subset and 'Y_FIRST' not in meta.keys():
        raise SystemExit('ERROR: --subset mode is NOT supported for time-series in radar-coordinates!')

    return inps


################################################################
def main(iargs=None):
    from ..save_hdfeos5 import read_template2inps, prep_metadata, get_output_filename, write_hdf5_file

    inps = cmd_line_parse(iargs)
    inps, template = read_template2inps(inps.template_file, inps)

    # Prepare Metadata
    meta = prep_metadata(
        ts_file=inps.ts_file,
        geom_file=inps.geom_file,
        template=template,
        print_msg=True)

    # Get output filename
    out_file = get_output_filename(
        metadata=meta,
        suffix=inps.suffix,
        update_mode=inps.update,
        subset_mode=inps.subset)

    # Open HDF5 File
    write_hdf5_file(
        metadata=meta,
        out_file=out_file,
        ts_file=inps.ts_file,
        tcoh_file=inps.tcoh_file,
        scoh_file=inps.scoh_file,
        mask_file=inps.mask_file,
        geom_file=inps.geom_file)


################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
