#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2018               #
############################################################


import os
import glob
import argparse
import numpy as np
from mintpy.utils import (
    ptime,
    readfile,
    writefile,
    isce_utils,
)


EXAMPLE = """example:
  # interferogram stack
  prep_isce.py -d ./merged/interferograms -m ./reference/IW1.xml -b ./baselines -g ./merged/geom_reference      #for topsStack
  prep_isce.py -d ./Igrams -m ./referenceShelve/data.dat -b ./baselines -g ./geom_reference                     #for stripmapStack
  prep_isce.py -m 20120507_slc_crop.xml -g ./geometry                                                           #for stripmapApp

  # 150408 is the reference date of alosStack processing
  prep_isce.py -d "pairs/*-*/insar" -m "pairs/*-*/150408.track.xml" -b baseline -g dates_resampled/150408/insar #for alosStack

  # offset stack from topsStack
  prep_isce.py -d ./merged/offsets -f filtAz*.off -m ./reference/IW1.xml -b ./baselines -g ./merged/offsets/geom_reference
"""

def create_parser():
    """Command line parser."""
    parser = argparse.ArgumentParser(description='Prepare ISCE metadata files.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)
    parser.add_argument('-d', '--ds-dir', '--dset-dir', dest='dsetDir', type=str, default=None, required=True,
                        help='The directory which contains all pairs\n'
                             'e.g.: $PROJECT_DIR/merged/interferograms OR \n'
                             '      $PROJECT_DIR/pairs/*-*/insar OR \n'
                             '      $PROJECT_DIR/merged/offsets')
    parser.add_argument('-f', '--file-pattern', nargs = '+', dest='dsetFiles', type=str, default=['filt_*.unw'],
                        help='List of observation files, e.g.: filt_fine.unw OR filtAz*.off')
    parser.add_argument('-m', '--meta-file', dest='metaFile', type=str, default=None, required=True,
                        help='Metadata file to extract common metada for the stack:\n'
                             'e.g.: for ISCE/topsStack    : reference/IW3.xml;\n'
                             '      for ISCE/stripmapStack: referenceShelve/data.dat;\n'
                             '      for ISCE/alosStack    : pairs/150408-150701/150408.track.xml\n'
                             '          where 150408 is the reference date of stack processing')
    parser.add_argument('-b', '--baseline-dir', dest='baselineDir', type=str, default=None,
                        help=' directory with baselines ')
    parser.add_argument('-g', '--geometry-dir', dest='geometryDir', type=str, default=None, required=True,
                        help=' directory with geometry files ')
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


def prepare_geometry(geom_dir, metadata=dict(), processor='tops', update_mode=True):
    """Prepare and extract metadata from geometry files"""

    print('prepare .rsc file for geometry files')
    geom_dir = os.path.abspath(geom_dir)

    # grab all existed files
    if processor in ['tops', 'stripmap']:
        fbases = ['hgt', 'lat', 'lon', 'los', 'shadowMask', 'incLocal']
        isce_files = [os.path.join(geom_dir, '{}.rdr'.format(i)) for i in fbases]
        isce_files = [i for i in isce_files if os.path.isfile(i)]
        if len(isce_files) == 0:
            isce_files = [os.path.join(geom_dir, '{}.rdr.full'.format(i)) for i in fbases]
            isce_files = [i for i in isce_files if os.path.isfile(i)]

    elif processor in ['alosStack']:
        alooks = metadata['ALOOKS']
        rlooks = metadata['RLOOKS']
        isce_files = glob.glob(os.path.join(geom_dir, '*_{}rlks_{}alks.hgt'.format(rlooks, alooks)))+\
                     glob.glob(os.path.join(geom_dir, '*_{}rlks_{}alks.lat'.format(rlooks, alooks)))+\
                     glob.glob(os.path.join(geom_dir, '*_{}rlks_{}alks.lon'.format(rlooks, alooks)))+\
                     glob.glob(os.path.join(geom_dir, '*_{}rlks_{}alks.los'.format(rlooks, alooks)))+\
                     glob.glob(os.path.join(geom_dir, '*_{}rlks_{}alks.wbd'.format(rlooks, alooks)))

    else:
        raise Exception('unknown processor: {}'.format(processor))

    # write rsc file for each file
    for isce_file in isce_files:
        # prepare metadata for current file
        if os.path.isfile(isce_file+'.xml'):
            geom_metadata = readfile.read_attribute(isce_file, metafile_ext='.xml')
        else:
            geom_metadata = readfile.read_attribute(isce_file)
        geom_metadata.update(metadata)

        # write .rsc file
        rsc_file = isce_file+'.rsc'
        writefile.write_roipac_rsc(geom_metadata, rsc_file,
                                   update_mode=update_mode,
                                   print_msg=True)
    return


def prepare_stack(inputDir, filePattern, metadata=dict(), baseline_dict=dict(), processor='tops', update_mode=True):
    print('prepare .rsc file for ', filePattern)
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
                         metadata=metadata,
                         processor=inps.processor,
                         update_mode=inps.update_mode)

    # read baseline info
    baseline_dict = {}
    if inps.baselineDir:
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
    """Main driver."""
    main()
