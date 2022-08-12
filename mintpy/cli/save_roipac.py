############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import sys
from mintpy.utils.arg_utils import create_argument_parser


##############################################################################
EXAMPLE = """example:
  #----- unwrapped phase
  #for velocity: output an interferogram with temporal baseline in DATE12 metadata
  save_roipac.py  velocity.h5
  save_roipac.py  velocity.h5 -m maskTempCoh.h5 maskAoiShinmoe.h5

  #for time-series: specify (date1_)date2
  save_roipac.py  timeseries_ERA5_ramp_demErr.h5  #use the last date
  save_roipac.py  timeseries_ERA5_ramp_demErr.h5  20050601
  save_roipac.py  timeseries_ERA5_ramp_demErr.h5  20040728_20050601

  #for HDF-EOS5: specify displacement-date1_date2
  save_roipac.py  S1_IW12_128_0593_0597_20141213_20180619.he5  displacement-20170904_20170916
  save_roipac.py  S1_IW12_128_0593_0597_20141213_20180619.he5  displacement-20170916

  #for ifgramStack: specify date1_date2
  save_roipac.py  inputs/ifgramStack.h5  unwrapPhase-20091225_20100723
  save_roipac.py  inputs/ifgramStack.h5  unwrapPhase-20091225_20100723  --ref-yx 640 810

  #----- coherence
  save_roipac.py  inputs/ifgramStack.h5  coherence-20091225_20100723
  save_roipac.py  temporalCoherence.h5
  save_roipac.py  S1_IW12_128_0593_0597_20141213_20180619.he5 temporalCoherence -o 20170904_20170916.cor

  #----- DEM
  save_roipac.py  geo_geometryRadar.h5  height -o srtm1.dem
  save_roipac.py  geo_geometryRadar.h5  height -o srtm1.hgt
  save_roipac.py  S1_IW12_128_0593_0597_20141213_20180619.he5 height -o srtm1.dem
"""


def create_parser(subparsers=None):
    synopsis = 'Convert MintPy HDF5 file to ROI_PAC format.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', help='HDF5 file to be converted.')
    parser.add_argument('dset', nargs='?', help='date/date12 of timeseries, or date12 of interferograms to be converted')
    parser.add_argument('-o', '--output', dest='outfile', help='output file name.')
    parser.add_argument('-m','--mask', dest='mask_file', nargs='+', help='mask file')
    parser.add_argument('--ref-yx', dest='ref_yx', type=int, nargs=2, help='custom reference pixel in y/x')
    parser.add_argument('--ref-lalo', dest='ref_lalo', type=float, nargs=2, help='custom reference pixel in lat/lon')
    parser.add_argument('--keep-all-metadata', dest='keepAllMetadata', action='store_true', help='Do not clean the metadata as ROIPAC format')
    return parser


def cmd_line_parse(iargs=None):
    from ..objects import timeseries
    from ..utils import readfile

    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # default dset
    if not inps.dset:
        atr = readfile.read_attribute(inps.file)
        k = atr['FILE_TYPE']
        if k in ['ifgramStack', 'HDFEOS']:
            raise Exception("NO input dataset! It's required for {} file".format(k))

        #for time-series
        if k == 'timeseries':
            inps.dset = timeseries(inps.file).get_date_list()[-1]
            print('NO date specified >>> continue with the last date: {}'.format(inps.dset))
    return inps


##############################################################################
def main(iargs=None):
    from ..utils import writefile
    from ..save_roipac import read_data, clean_metadata4roipac

    inps = cmd_line_parse(iargs)

    data, atr, out_file = read_data(inps)

    if not inps.keepAllMetadata:
        atr = clean_metadata4roipac(atr)

    writefile.write(data, out_file=out_file, metadata=atr)


##########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
