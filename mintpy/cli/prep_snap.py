############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import os
import sys
from mintpy.utils.arg_utils import create_argument_parser


##################################################################################################
NOTE = """
  For each interferogram, coherence or unwrapped .dim product this script will prepare.rsc 
  metadata files for for mintpy based on .dim metadata file.

  The SNAP .dim file should contain all the required sensor / baseline metadata needed.
  The baseline metadata gets written during snap back-geocoding (co-registration).
  prep_snap is run separately for unw/ifg/cor files so needs separate .dim/.data products
  with only the relevant band in each product. Use Band Subset > save BEAM-DIMAP file.

  The file name should be yyyymmdd_yyyymmdd_type_tc.dim where type can be filt/unw/coh.

  The DEM should be prepared by adding an elevation file to a coregestered product - 
  then extract the elevation band only. Use Band Subset > save BEAM-DIMAP file

  Currently only works for geocoded (terrain correction step in SNAP) interferograms.
"""

EXAMPLE = """example:
  prep_snap.py  ../interferograms/*/*/Unw_*.img
  prep_snap.py  ../dem_tc.data/dem*.img
"""

def create_parser(subparsers=None):
    synopsis = 'Prepare attributes file for SNAP products.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis+NOTE, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', nargs='+', help='SNAP data file(s) in *.img format.')
    return parser


def cmd_line_parse(iargs=None):
    from mintpy.utils import utils as ut

    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    inps.file = ut.get_file_list(inps.file, abspath=True)
    for fname in inps.file:
        if not fname.endswith('.img'):
            raise ValueError('Input data file does NOT end with .img: {}'.format(fname))

    return inps



##################################################################################################
def main(iargs=None):
    from mintpy.utils import readfile
    from mintpy.prep_snap import write_rsc

    inps = cmd_line_parse(iargs)

    for img_file in inps.file:
        # read metadata from *.dim file
        # the map info from *.img.hdr file is NOT right, thus, not used.
        dim_file = os.path.dirname(img_file)[:-4] + 'dim'
        atr = readfile.read_snap_dim(dim_file)

        # write metadata dict to *.rsc file
        rsc_file = img_file + '.rsc'
        write_rsc(atr, rsc_file)

    return


##################################################################################################
if __name__ == "__main__":
    main(sys.argv[1:])
