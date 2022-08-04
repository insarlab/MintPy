#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Forrest Williams, Mar 2021                       #
############################################################


import os
import sys
from datetime import datetime
import numpy as np

from mintpy.utils import readfile, writefile, utils as ut
from mintpy.utils.arg_utils import create_argument_parser


#########################################################################
EXAMPLE_META_FILE = """
offset1NS.tif  20160206 20161122
offset1EW.tif  20160206 20161122
offset1SNR.tif 20160206 20161122
offset2NS.tif  20160206 20170225
offset2EW.tif  20160206 20170225
offset2SNR.tif 20160206 20170225
...            ...   ...
"""

EXAMPLE = """example:
  prep_cosicorr.py offsets/*offset.tif -m metadata.txt
  prep_cosicorr.py snr/*snr.tif        -m metadata.txt
"""

def create_parser(subparsers=None):
    """Command line parser."""
    synopsis = 'Prepare attributes file for COSI-Corr pixel offset product.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', nargs='+', help='cosicorr file(s)')
    parser.add_argument('-m', '--metadata', type=str, dest='meta_file',
                        help='metadata file with date info. E.g.:'+EXAMPLE_META_FILE)
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    inps.file = ut.get_file_list(inps.file, abspath=True)
    return inps


#########################################################################
def add_cosicorr_metadata(fname, cosicorr_dates, meta):
    '''Read/extract attribute data from cosicorr metadata file and add to metadata dictionary
    Inputs:
        Offset or SNR file name (fname)
        dictionary of file name and date12 pairs (cosicorr_dates)
        Metadata dictionary (meta)
    Output:
        Metadata dictionary (meta)
    '''

    # add general attributes
    meta['PROCESSOR'] = 'cosicorr'
    meta['P_BASELINE_TOP_HDR'] = 0.0 #placeholder
    meta['P_BASELINE_BOTTOM_HDR'] = 0.0 #placeholder
    meta['RANGE_PIXEL_SIZE'] = np.abs(meta['X_STEP'])
    meta['AZIMUTH_PIXEL_SIZE'] = np.abs(meta['Y_STEP'])
    meta['RLOOKS'] = 1
    meta['ALOOKS'] = 1

    # Time attributes
    date1_string, date2_string = cosicorr_dates[os.path.basename(fname)].split('-')
    meta['DATE12'] = f'{date1_string}-{date2_string}'
    date1 = datetime.strptime(date1_string,'%Y%m%d')
    date2 = datetime.strptime(date2_string,'%Y%m%d')
    date_avg = date1 + (date2 - date1) / 2
    date_avg_seconds = (date_avg - date_avg.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds()
    meta['CENTER_LINE_UTC'] = date_avg_seconds

    # add LAT/LON_REF1/2/3/4
    N = float(meta['Y_FIRST'])
    W = float(meta['X_FIRST'])
    S = N + float(meta['Y_STEP']) * int(meta['LENGTH'])
    E = W + float(meta['X_STEP']) * int(meta['WIDTH'])

    meta['LAT_REF1'] = str(S)
    meta['LAT_REF2'] = str(S)
    meta['LAT_REF3'] = str(N)
    meta['LAT_REF4'] = str(N)
    meta['LON_REF1'] = str(W)
    meta['LON_REF2'] = str(E)
    meta['LON_REF3'] = str(W)
    meta['LON_REF4'] = str(E)

    return(meta)


#########################################################################
def main(iargs=None):
    # read in arguments
    inps = cmd_line_parse(iargs)

    # open and read hyp3 metadata
    cosicorr_dates = {}
    with open(inps.meta_file, 'r') as f:
        for line in f:
            name, date1, date2 = line.strip().split(' ')
            cosicorr_dates[name] = f'{date1}-{date2}'

    # for each filename, generate metadata rsc file
    for fname in inps.file:
        meta = readfile.read_gdal_vrt(fname)
        meta = add_cosicorr_metadata(fname, cosicorr_dates, meta)

        rsc_file = fname+'.rsc'
        writefile.write_roipac_rsc(meta, out_file=rsc_file)
    return


#########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
