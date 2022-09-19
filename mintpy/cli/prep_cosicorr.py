############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import sys
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
    from mintpy.utils import utils as ut
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    inps.file = ut.get_file_list(inps.file, abspath=True)
    return inps


#########################################################################
def main(iargs=None):
    from mintpy.utils import readfile, writefile
    from mintpy.prep_cosicorr import add_cosicorr_metadata

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


#########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
