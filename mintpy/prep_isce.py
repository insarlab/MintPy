#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2018               #
############################################################


import os
import sys
import glob
import numpy as np
from mintpy.utils import (
    attribute as attr,
    isce_utils,
    ptime,
    readfile,
    writefile,
)
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
def add_ifgram_metadata(metadata_in, dates=[], baseline_dict={}):
    """Add metadata unique for each interferogram
    Parameters: metadata_in   : dict, input common metadata for the entire dataset
                dates         : list of str in YYYYMMDD or YYMMDD format
                baseline_dict : dict, output of baseline_timeseries()
    Returns:    metadata      : dict, updated metadata
    """
    # make a copy of input metadata
    metadata = {}
    for k in metadata_in.keys():
        metadata[k] = metadata_in[k]

    # DATE12
    metadata['DATE12'] = '{}-{}'.format(dates[0][2:], dates[1][2:])

    # P_BASELINE*
    if baseline_dict:
        bperp_top = baseline_dict[dates[1]][0] - baseline_dict[dates[0]][0]
        bperp_bottom = baseline_dict[dates[1]][1] - baseline_dict[dates[0]][1]
        metadata['P_BASELINE_TOP_HDR'] = str(bperp_top)
        metadata['P_BASELINE_BOTTOM_HDR'] = str(bperp_bottom)
    return metadata


def prepare_geometry(geom_dir, geom_files=[], metadata=dict(), processor='tops', update_mode=True):
    """Prepare and extract metadata from geometry files.

    Parameters: geom_dir   - str, path to the directorry for the geometry data files
                geom_files - list(str), basenames of geometry data files
                metadata   - dict, common metadata for the stack
                processor  - str, isce-2 stack processor
    """

    print('preparing RSC file for geometry files')
    geom_dir = os.path.abspath(geom_dir)

    # default file basenames
    if not geom_files:
        if processor in ['tops', 'stripmap']:
            geom_files = ['{}.rdr'.format(i) for i in GEOMETRY_PREFIXS]

        elif processor in ['alosStack']:
            alooks = metadata['ALOOKS']
            rlooks = metadata['RLOOKS']
            fexts = ['.hgt', '.lat', '.lon', '.los', '.wbd']
            geom_files = ['*_{}rlks_{}alks{}'.format(rlooks, alooks, fext) for fext in fexts]

        else:
            raise Exception('unknown processor: {}'.format(processor))

    # get absolute file paths
    geom_files = [os.path.join(geom_dir, i) for i in geom_files]

    # check the full resolution version if no multilooked version exists
    if all(not os.path.isfile(i) for i in geom_files):
        geom_files = [i+'.full' for i in geom_files]

    # get existed files
    geom_files = [i for i in geom_files if os.path.isfile(i)]
    # remove duplicates while preserving order
    seen = set()
    seen_add = seen.add
    geom_files = [i for i in geom_files if not (i in seen or seen_add(i))]

    # write rsc file for each file
    for geom_file in geom_files:
        # prepare metadata for current file
        geom_meta = {**metadata}
        if os.path.isfile(geom_file+'.xml'):
            geom_meta.update(readfile.read_attribute(geom_file, metafile_ext='.xml'))
        else:
            geom_meta.update(readfile.read_attribute(geom_file))

        # write .rsc file
        rsc_file = geom_file+'.rsc'
        writefile.write_roipac_rsc(geom_meta, rsc_file,
                                   update_mode=update_mode,
                                   print_msg=True)
    return


def gen_random_baseline_timeseries(obs_file, max_bperp=10):
    """Generate a baseline time series with random values
    with date12 values grabbed from the directory names of the given path pattern from obs_file.
    """
    # list of dates
    date12s = sorted([os.path.basename(os.path.dirname(x)) for x in glob.glob(obs_file)])
    date1s = [x.split('_')[0] for x in date12s]
    date2s = [x.split('_')[1] for x in date12s]
    date_list = sorted(list(set(date1s + date2s)))

    # list of bperp
    bperp_list = [0] + np.random.randint(-max_bperp, max_bperp, len(date_list)-1).tolist()

    # prepare output
    bDict = {}
    for date_str, bperp in zip(date_list, bperp_list):
        bDict[date_str] = [bperp, bperp]

    return bDict


def prepare_stack(obs_file, metadata=dict(), baseline_dict=dict(), update_mode=True):
    """Prepare metadata for a stack of observation data files.

    Parameters: obs_file     : path pattern with wildcards for the primary observation files, e.g. *.unw file.
                metadata     : dict, common metadata for the stack
                baseline_dir : dict, baseline time series
    """
    print(f'preparing RSC file for: {obs_file}')
    isce_files = sorted(glob.glob(obs_file))
    if len(isce_files) == 0:
        raise FileNotFoundError('NO file found with path pattern: {}'.format(obs_file))

    # make a copy
    meta = {**metadata}

    # update A/RLOOKS, RANGE/AZIMUTH_PIXEL_SIZE, NCORRLOOKS
    # for low resolution ionosphere from isce2/topsStack
    keys = ['LENGTH', 'WIDTH']
    if all(x in meta.keys() for x in keys):
        atr = readfile.read_attribute(isce_files[0], metafile_ext='.xml')
        if any(int(meta[x]) != int(atr[x]) for x in keys):
            resize2shape = (int(atr['LENGTH']), int(atr['WIDTH']))
            meta = attr.update_attribute4resize(meta, resize2shape)

    # write .rsc file for each interferogram file
    num_file = len(isce_files)
    print_msg = True if num_file > 5 else False   # do not print progress bar for <=5 files
    prog_bar = ptime.progressBar(maxValue=num_file, print_msg=print_msg)
    for i, isce_file in enumerate(isce_files):
        # get date1/2
        date12 = ptime.get_date12_from_path(isce_file)
        dates = ptime.yyyymmdd(date12.replace('-','_').split('_'))
        prog_bar.update(i+1, suffix=f'{dates[0]}_{dates[1]} {i+1}/{num_file}')

        # merge metadata from: data.rsc, *.unw.xml and DATE12/P_BASELINE_TOP/BOTTOM_HDR
        ifg_meta = {**meta}
        ifg_meta.update(readfile.read_attribute(isce_file, metafile_ext='.xml'))
        ifg_meta = add_ifgram_metadata(ifg_meta, dates, baseline_dict)

        # write .rsc file
        rsc_file = isce_file+'.rsc'
        writefile.write_roipac_rsc(ifg_meta, rsc_file,
                                   update_mode=update_mode,
                                   print_msg=False)
    prog_bar.close()
    return


#########################################################################
def main(iargs=None):
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
    return


#########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
