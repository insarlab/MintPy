############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import sys
from mintpy.utils.arg_utils import create_argument_parser


#########################################################################
NOTE = """
  For each interferogram, the unwrapped interferogram, coherence, and metadata the file name is required e.g.:
  1) S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_unw_phase.tif
  2) S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_corr.tif
  3) S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2.txt

  A DEM filename is needed and a incidence angle filename is recommended  e.g.:
  1) S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_dem.tif
  2) S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_lv_theta.tif

  This script will read these files, read the geospatial metadata from GDAL,
  find the corresponding HyP3 metadata file (for interferograms and coherence),
  and write to a ROI_PAC .rsc metadata file with the same name as the input file with suffix .rsc,
  e.g. S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_unw_phase.tif.rsc

  Here is an example of how your HyP3 files should look:

  Before loading:
      For each interferogram, 3 files are needed:
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_unw_phase_clip.tif
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_corr_clip.tif
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2.txt
      For the geometry file 2 file are recommended:
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_dem_clip.tif     (required)
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_lv_theta_clip.tif (optional but recommended)

  After running prep_hyp3.py:
      For each interferogram:
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_unw_phase_clip.tif
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_unw_phase_clip.tif.rsc
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_corr_clip.tif
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_corr_clip.tif.rsc
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2.txt
      For the input geometry files:
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_dem_clip.tif
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_dem_clip.tif.rsc
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_lv_theta_clip.tif
          S1AA_20161223T070700_20170116T070658_VVP024_INT80_G_ueF_74C2_lv_theta_clip.tif.rsc

  Notes:
    HyP3 currently only supports generation of Sentinel-1 interferograms, so
    some Sentinel-1 metadata is hard-coded. If HyP3 adds processing of interferograms
    from other satellites, changes will be needed. 
"""

EXAMPLE = """example:
  prep_hyp3.py  interferograms/*/*unw_phase_clip.tif
  prep_hyp3.py  interferograms/*/*corr_clip.tif
  prep_hyp3.py  interferograms/*/*dem_clip.tif
  prep_hyp3.py  interferograms/*/*lv_theta_clip.tif
  prep_hyp3.py  interferograms/*/*clip.tif
"""


def create_parser(subparsers=None):
    synopsis = 'Prepare attributes file for HyP3 InSAR product.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis+NOTE, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', nargs='+', help='HyP3 file(s)')
    return parser


def cmd_line_parse(iargs=None):
    from ..utils import utils as ut
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    inps.file = ut.get_file_list(inps.file, abspath=True)
    return inps


#########################################################################
def main(iargs=None):
    from ..utils import readfile, writefile
    from ..prep_hyp3 import add_hyp3_metadata

    # read in arguments
    inps = cmd_line_parse(iargs)

    # for each filename, generate metadata rsc file
    for fname in inps.file:
        is_ifg = any([x in fname for x in ['unw_phase','corr']])
        meta = readfile.read_gdal_vrt(fname)
        meta = add_hyp3_metadata(fname, meta, is_ifg=is_ifg)

        # write
        rsc_file = fname+'.rsc'
        writefile.write_roipac_rsc(meta, out_file=rsc_file)

    return


###################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
