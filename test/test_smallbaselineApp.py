#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Mar 2019                           #
############################################################
# Test example dataset for smallbaselineApp workflow


import os
import sys
from pathlib import Path
import time
import shutil
import argparse
import subprocess
import tarfile

import mintpy
from mintpy import view


CMAP_DICT = {
    'FernandinaSenDT128' : 'jet',
    'WellsEnvD2T399'     : 'jet_r',
    'KujuAlosAT422F650'  : 'jet_r',
}

URL_LIST = [
    'https://zenodo.org/record/3952953/files/FernandinaSenDT128.tar.xz',
    'https://zenodo.org/record/4320005/files/SanFranSenDT42.tar.xz',
    'https://zenodo.org/record/3952950/files/WellsEnvD2T399.tar.xz',
    'https://zenodo.org/record/4318134/files/WCapeSenAT29.tar.xz',
    'https://zenodo.org/record/3952917/files/KujuAlosAT422F650.tar.xz',
]

PROJ_NAME_LIST = [os.path.basename(url).split('.tar.xz')[0] for url in URL_LIST]


#####################################################################################
DSET_INFO = """
  FernandinaSenDT128 for ISCE/topsStack
  SanFranSenDT42     for ARIA
  WellsEnvD2T399     for GAMMA
  WCapeSenAT29       for SNAP
  KujuAlosAT422F650  for ROI_PAC
"""

EXAMPLE = """example:
  $MINTPY_HOME/test/test_smallbaselineApp.py
  $MINTPY_HOME/test/test_smallbaselineApp.py  --dir ~/test
  $MINTPY_HOME/test/test_smallbaselineApp.py  --dset KujuAlosAT422F650
  $MINTPY_HOME/test/test_smallbaselineApp.py  --nofresh
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Test smallbaselineApp workflow with example datasets.'+DSET_INFO,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('--dset', dest='dset_name', nargs='+', metavar='DSET', choices=PROJ_NAME_LIST, default=PROJ_NAME_LIST,
                        help='name(s) of datasets to be tested.\n'+
                             'Available datasets: {}\n'.format(PROJ_NAME_LIST) +
                             '(default: all)')

    parser.add_argument('--test-pyaps', dest='test_pyaps', action='store_true',
                        help='Include testing of PyAPS (default: %(default)s).')

    parser.add_argument('--test-isce', dest='test_isce', action='store_true',
                        help='Include testing of ISCE-2/topsStack (default: %(default)s).')

    parser.add_argument('--dir', dest='test_dir', default='~/data/test',
                        help='test directory (default: %(default)s).')

    parser.add_argument('--nofresh', dest='fresh_start', action='store_false',
                        help='Use exsiting files WITHOUT starting from the scratch (default: %(default)s).')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # expand test_dir
    inps.test_dir = os.path.expanduser(inps.test_dir)
    inps.test_dir = os.path.expandvars(inps.test_dir)

    return inps


#####################################################################################


def test_dataset(dset_name, test_dir, fresh_start=True, test_pyaps=False, test_isce=False):
    print('Go to test directory:', test_dir)
    os.chdir(test_dir)

    dset_idx = PROJ_NAME_LIST.index(dset_name)
    dset_url = URL_LIST[dset_idx]
    mintpy_dir = os.path.dirname(mintpy.__file__)
    template_file = os.path.join(Path(mintpy_dir).parent, 'test/configs/{}.txt'.format(dset_name))

    # download tar file
    tar_file = os.path.basename(dset_url)
    if not os.path.isfile(tar_file):
        print('downloading tar file ...')
        cmd = 'wget {}'.format(dset_url)
        print(cmd)
        os.system(cmd)
    else:
        print('tar file exists, skip re-downloading.')

    # uncompress tar file
    if not fresh_start and os.path.isdir(dset_name):
        print('use existing files in project directory: {}'.format(dset_name))
    else:
        # remove existing directory
        if os.path.isdir(dset_name):
            print('remove existing project directory: {}'.format(dset_name))
            shutil.rmtree(dset_name)

        # uncompress tar file
        print('extracting content from tar file: {}'.format(tar_file))
        tar = tarfile.open(tar_file)
        tar.extractall()
        tar.close()

    # set working directory
    work_dir = os.path.join(test_dir, dset_name, 'mintpy')
    os.chdir(work_dir)
    print('Go to work directory:', work_dir)

    # test pyaps: remove existing tropo delay file
    tropo_file = './inputs/ERA5.h5'
    if test_pyaps and os.path.isfile(tropo_file):
        print('remove existing tropospheric delay file:', tropo_file)
        os.remove(tropo_file)

    # test isce: remove existing metadata file
    meta_file = '../reference/data.rsc'   #for ISCE-2/topsStack
    if test_isce and os.path.isfile(meta_file):
        print('remove existing metadata file:', meta_file)
        os.remove(meta_file)

    # runing smallbaselineApp
    # Note: execute script in command line instead of call main() for a clean run
    # to avoid strange error from prep_aria: not recognized as a supported file format.
    # which only occurs if all datasets are tested in one run
    cmd = 'smallbaselineApp.py {}'.format(template_file)
    print(cmd)
    status = subprocess.Popen(cmd, shell=True).wait()
    if status != 0:
        raise RuntimeError('Test failed for example dataset {}'.format(dset_name))

    # custom plot of velocity map
    vel_files = [os.path.join(work_dir, i) for i in ['geo/geo_velocity.h5', 'velocity.h5']]
    vel_file = [i for i in vel_files if os.path.isfile(i)][0]
    png_file = os.path.join(work_dir, 'pic', '{}.png'.format(os.path.splitext(os.path.basename(vel_file))[0]))
    iargs = [vel_file, 'velocity', '--nodisplay', '--noverbose', '-o', png_file]
    if dset_name in CMAP_DICT.keys():
        iargs += ['-c', CMAP_DICT[dset_name]]
    print('view.py', ' '.join(iargs))
    view.main(iargs)

    # open final velocity map if on mac
    if sys.platform.lower().startswith('darwin'):
        cmd = 'open {}'.format(png_file)
        print(cmd)
        subprocess.Popen(cmd, shell=True).wait()
    return


#####################################################################################
def main(iargs=None):
    start_time = time.time()
    inps = cmd_line_parse(iargs)

    # create test directory
    os.makedirs(inps.test_dir, exist_ok=True)

    # run test
    num_dset = len(inps.dset_name)
    for i in range(num_dset):
        dset_name = inps.dset_name[i]
        print('#'*100)
        print('Start testing smallbaselineApp workflow on exmaple dataset {}/{}: {}'.format(i+1, num_dset, dset_name))
        test_dataset(dset_name,
                     test_dir=inps.test_dir,
                     fresh_start=inps.fresh_start,
                     test_pyaps=inps.test_pyaps,
                     test_isce=inps.test_isce)
        print('#'*100)
        print('   PASS testing of smallbaselineApp workflow on exmaple dataset {}/{}: {}'.format(i+1, num_dset, dset_name))
        print('#'*100+'\n'*3)

    # print message
    if num_dset == len(PROJ_NAME_LIST):
        m, s = divmod(time.time()-start_time, 60)
        msg  = '#'*50
        msg += '\n    PASS ALL testings without running errors.\n'
        msg += '#'*50
        msg += '\nTotal time used: {:02.0f} mins {:02.1f} secs\n'.format(m, s)
        print(msg)

    return


#####################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])

