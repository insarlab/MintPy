############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import os
import sys
import glob
from mintpy.utils.arg_utils import create_argument_parser


#########################################################################
GEOMETRY_PREFIXS = ['hgt', 'lat', 'lon', 'los', 'shadowMask', 'waterMask', 'incLocal']

EXAMPLE = """example:
  ## topsStack
  prep_isce.py -f "./merged/interferograms/*/filt_*.unw" -m ./reference/IW1.xml -b ./baselines/ -g ./merged/geom_reference/

  # topsStack with ionosphere
  prep_isce.py -f "./merged/interferograms/*/filt_*.unw" "./ion/*/ion_cal/filt.ion" -m ./reference/IW1.xml -b ./baselines/ -g ./merged/geom_reference/

  # topsStack for offset
  prep_isce.py -f "./merged/offsets/*/*Off*.bip" -m ./reference/IW1.xml -b ./baselines/ -g ./merged/geom_reference/

  ## stripmapStack
  prep_isce.py -f "./Igrams/*/filt_*.unw" -m ./referenceShelve/data.dat -b ./baselines/ -g ./geom_reference/

  # stripmapApp
  prep_isce.py -m 20120507_slc_crop.xml -g ./geometry 

  ## alosStack
  # where 150408 is the reference date
  prep_isce.py -f "./pairs/*/insar/filt_*.unw" -m "pairs/150408-*/150408.track.xml" -b ./baseline/ -g ./dates_resampled/150408/insar/

  ## UAVSAR
  prep_isce.py -f "./Igrams/*/filt_*.unw" -m ./referenceShelve/data.dat -b ./baselines/ -g ./geometry/

  # UAVSAR for offset
  prep_isce.py -f "./offsets/*/*Off*.bip" -m "SLC/*/data.dat" -b random -g ./geometry/
"""

def create_parser(subparsers=None):
    """Command line parser."""
    synopsis = 'Prepare ISCE metadata files.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    # observations
    parser.add_argument('-f', dest='obs_files', type=str, nargs='+', default='./merged/interferograms/*/filt_*.unw',
                        help='Wildcard path pattern for the primary observation files.\n'
                             'E.g.: topsStack          : {dset_dir}/merged/interferograms/*/filt_*.unw\n'
                             '      topsStack / iono   : {dset_dir}/ion/*/ion_cal/filt.ion\n'
                             '      topsStack / offset : {dset_dir}/merged/offsets/*/*Off*.bip\n'
                             '      stripmapStack      : {dset_dir}/Igrams/*_*/filt_*.unw\n'
                             '      alosStack          : {dset_dir}/pairs/*/insar/filt_*.unw\n'
                             '      UAVSAR / offset    : {dset_dir}/offsets/*/*Off*.bip')

    # metadata
    parser.add_argument('-m', '--meta-file', dest='meta_file', type=str, default=None, required=True,
                        help='Metadata file to extract common metada for the stack.\n'
                             'E.g.: topsStack     : reference/IW3.xml\n'
                             '      stripmapStack : referenceShelve/data.dat\n'
                             '      alosStack     : pairs/{ref_date}-*/{ref_date}.track.xml\n'
                             '      UAVSAR        : SLC/*/data.dat')

    # geometry
    parser.add_argument('-b', '--baseline-dir', dest='baseline_dir', type=str, default=None,
                        help='Directory with baselines. '
                             'Set "random" to generate baseline with random value from [-10,10].')
    parser.add_argument('-g', '--geometry-dir', dest='geom_dir', type=str, default=None, required=True,
                        help='Directory with geometry files ')
    parser.add_argument('--geom-files', dest='geom_files', type=str, nargs='*',
                        default=['{}.rdr'.format(i) for i in GEOMETRY_PREFIXS],
                        help='List of geometry file basenames. Default: %(default)s.\n'
                             'All geometry files need to be in the same directory.')

    parser.add_argument('--force', dest='update_mode', action='store_false',
                        help='Force to overwrite all .rsc metadata files.')
    return parser


def cmd_line_parse(iargs = None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # translate wildcard in meta_file
    if "*" in inps.meta_file:
        fnames = glob.glob(inps.meta_file)
        if len(fnames) > 0:
            inps.meta_file = fnames[0]
        else:
            raise FileNotFoundError(inps.meta_file)

    return inps


#########################################################################
def main(iargs=None):
    from mintpy.utils import isce_utils
    from mintpy.prep_isce import prepare_geometry, prepare_stack, gen_random_baseline_timeseries

    inps = cmd_line_parse(iargs)
    inps.processor = isce_utils.get_processor(inps.meta_file)

    # read common metadata
    metadata = {}
    if inps.meta_file:
        rsc_file = os.path.join(os.path.dirname(inps.meta_file), 'data.rsc')
        metadata = isce_utils.extract_isce_metadata(
            inps.meta_file,
            geom_dir=inps.geom_dir,
            rsc_file=rsc_file,
            update_mode=inps.update_mode)[0]

    # prepare metadata for geometry file
    if inps.geom_dir:
        prepare_geometry(
            inps.geom_dir,
            geom_files=inps.geom_files,
            metadata=metadata,
            processor=inps.processor,
            update_mode=inps.update_mode)

    # read baseline info
    baseline_dict = {}
    if inps.baseline_dir:
        if inps.baseline_dir.startswith('rand') and inps.obs_files:
            baseline_dict = gen_random_baseline_timeseries(inps.obs_files[0])
        else:
            baseline_dict = isce_utils.read_baseline_timeseries(
                inps.baseline_dir,
                processor=inps.processor)

    # prepare metadata for ifgram file
    if inps.obs_files:
        for obs_file in inps.obs_files:
            prepare_stack(
                obs_file,
                metadata=metadata,
                baseline_dict=baseline_dict,
                update_mode=inps.update_mode)

    print('Done.')


#########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
