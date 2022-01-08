#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2018               #
############################################################


import os
import sys
import glob
import argparse
import numpy as np
from mintpy.utils import ptime, readfile, writefile, isce_utils



#########################################################################
GEOMETRY_PREFIXS = ['hgt', 'lat', 'lon', 'los', 'shadowMask', 'waterMask', 'incLocal']

EXAMPLE = """example:
  # interferogram stack
  prep_isce.py -d ./merged/interferograms -m ./reference/IW1.xml -b ./baselines -g ./merged/geom_reference      #for topsStack
  prep_isce.py -d ./Igrams -m ./referenceShelve/data.dat -b ./baselines -g ./geom_reference                     #for stripmapStack
  prep_isce.py -m 20120507_slc_crop.xml -g ./geometry                                                           #for stripmapApp
  prep_isce.py -d "pairs/*-*/insar" -m "pairs/*-*/150408.track.xml" -b baseline -g dates_resampled/150408/insar #for alosStack with 150408 as ref date

  # offset stack
  prep_isce.py -d ./offsets -f *Off*.bip -m ./../reference/IW1.xml -b ./../baselines -g ./offsets/geom_reference  #for topsStack
  prep_isce.py -d ./offsets -f *Off*.bip -m ./SLC/*/data.dat       -b random         -g ./geometry                #for UAVSAR coregStack
"""

def create_parser():
    """Command line parser."""
    parser = argparse.ArgumentParser(description='Prepare ISCE metadata files.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)
    # interferograms
    parser.add_argument('-d', '--ds-dir', '--dset-dir', dest='dsetDir', type=str, default=None, required=True,
                        help='The directory which contains all pairs\n'
                             'e.g.: $PROJECT_DIR/merged/interferograms OR \n'
                             '      $PROJECT_DIR/pairs/*-*/insar OR \n'
                             '      $PROJECT_DIR/merged/offsets')
    parser.add_argument('-f', '--file-pattern', nargs = '+', dest='dsetFiles', type=str, default=['filt_*.unw'],
                        help='List of observation file basenames, e.g.: filt_fine.unw OR filtAz*.off')

    # metadata
    parser.add_argument('-m', '--meta-file', dest='metaFile', type=str, default=None, required=True,
                        help='Metadata file to extract common metada for the stack:\n'
                             'e.g.: for ISCE/topsStack    : reference/IW3.xml;\n'
                             '      for ISCE/stripmapStack: referenceShelve/data.dat;\n'
                             '      for ISCE/alosStack    : pairs/150408-150701/150408.track.xml\n'
                             '          where 150408 is the reference date of stack processing')

    # geometry
    parser.add_argument('-b', '--baseline-dir', dest='baselineDir', type=str, default=None,
                        help='Directory with baselines.'
                             'Set "random" to generate baseline with random value from [-10,10].'
                             'Set "random-100" to generate baseline with random value from [-100,100].')
    parser.add_argument('-g', '--geometry-dir', dest='geometryDir', type=str, default=None, required=True,
                        help='Directory with geometry files ')
    parser.add_argument('--geom-files', dest='geometryFiles', type=str, nargs='*',
                        default=['{}.rdr'.format(i) for i in GEOMETRY_PREFIXS],
                        help='List of geometry file basenames. Default: %(default)s.\n'
                             'All geometry files need to be in the same directory.')

    parser.add_argument('--force', dest='update_mode', action='store_false',
                        help='Force to overwrite all .rsc metadata files.')
    return parser


def cmd_line_parse(iargs = None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # translate wildcard in metaFile
    if "*" in inps.metaFile:
        fnames = glob.glob(inps.metaFile)
        if len(fnames) > 0:
            inps.metaFile = fnames[0]
        else:
            raise FileNotFoundError(inps.metaFile)

    # random baseline input checking
    if inps.baselineDir.lower().startswith('rand'):
        inps.baselineDir = inps.baselineDir.lower().replace('_','-')

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
    """Prepare and extract metadata from geometry files"""

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
        if os.path.isfile(geom_file+'.xml'):
            geom_metadata = readfile.read_attribute(geom_file, metafile_ext='.xml')
        else:
            geom_metadata = readfile.read_attribute(geom_file)
        geom_metadata.update(metadata)

        # write .rsc file
        rsc_file = geom_file+'.rsc'
        writefile.write_roipac_rsc(geom_metadata, rsc_file,
                                   update_mode=update_mode,
                                   print_msg=True)
    return


def gen_random_baseline_timeseries(dset_dir, dset_file, max_bperp=10):
    """Generate a baseline time series with random values.
    """
    # list of dates
    fnames = glob.glob(os.path.join(dset_dir, '*', dset_file))
    date12s = sorted([os.path.basename(os.path.dirname(x)) for x in fnames])
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


def prepare_stack(inputDir, filePattern, metadata=dict(), baseline_dict=dict(), processor='tops', update_mode=True):
    print('preparing RSC file for ', filePattern)
    if processor in ['tops', 'stripmap']:
        isce_files = sorted(glob.glob(os.path.join(os.path.abspath(inputDir), '*', filePattern)))
    elif processor == 'alosStack':
        isce_files = sorted(glob.glob(os.path.join(os.path.abspath(inputDir), filePattern)))
    else:
        raise ValueError('Un-recognized ISCE stack processor: {}'.format(processor))

    if len(isce_files) == 0:
        raise FileNotFoundError('no file found in pattern: {}'.format(filePattern))

    # write .rsc file for each interferogram file
    num_file = len(isce_files)
    prog_bar = ptime.progressBar(maxValue=num_file)
    for i in range(num_file):
        # prepare metadata for current file
        isce_file = isce_files[i]
        if processor in ['tops', 'stripmap']:
            dates = os.path.basename(os.path.dirname(isce_file)).split('_')  # to modify to YYYYMMDDTHHMMSS
        elif processor == 'alosStack':
            dates = os.path.basename(os.path.dirname(os.path.dirname(isce_file))).split('-')  # to modify to YYYYMMDDTHHMMSS
            dates = ptime.yyyymmdd(dates)
        else:
            raise ValueError('Un-recognized ISCE stack processor: {}'.format(processor))

        ifg_metadata = readfile.read_attribute(isce_file, metafile_ext='.xml')
        ifg_metadata.update(metadata)
        ifg_metadata = add_ifgram_metadata(ifg_metadata, dates, baseline_dict)

        # write .rsc file
        rsc_file = isce_file+'.rsc'
        writefile.write_roipac_rsc(ifg_metadata, rsc_file,
                                   update_mode=update_mode,
                                   print_msg=False)
        prog_bar.update(i+1, suffix='{}_{}'.format(dates[0], dates[1]))
    prog_bar.close()
    return


#########################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    inps.processor = isce_utils.get_processor(inps.metaFile)

    # read common metadata
    metadata = {}
    if inps.metaFile:
        rsc_file = os.path.join(os.path.dirname(inps.metaFile), 'data.rsc')
        metadata = isce_utils.extract_isce_metadata(inps.metaFile,
                                                    geom_dir=inps.geometryDir,
                                                    rsc_file=rsc_file,
                                                    update_mode=inps.update_mode)[0]

    # prepare metadata for geometry file
    if inps.geometryDir:
        prepare_geometry(inps.geometryDir,
                         geom_files=inps.geometryFiles,
                         metadata=metadata,
                         processor=inps.processor,
                         update_mode=inps.update_mode)

    # read baseline info
    baseline_dict = {}
    if inps.baselineDir:
        if inps.baselineDir.startswith('rand') and inps.dsetDir and inps.dsetFiles:
            if '-' in inps.baselineDir:
                max_bperp = float(inps.baselineDir.split('-')[1])
            else:
                max_bperp = 10
            baseline_dict = gen_random_baseline_timeseries(dset_dir=inps.dsetDir,
                                                           dset_file=inps.dsetFiles[0],
                                                           max_bperp=max_bperp)

        else:
            baseline_dict = isce_utils.read_baseline_timeseries(inps.baselineDir,
                                                                processor=inps.processor)

    # prepare metadata for ifgram file
    if inps.dsetDir and inps.dsetFiles:
        for namePattern in inps.dsetFiles:
            prepare_stack(inps.dsetDir,
                          namePattern,
                          metadata=metadata,
                          baseline_dict=baseline_dict,
                          processor=inps.processor,
                          update_mode=inps.update_mode)
    print('Done.')
    return


#########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
