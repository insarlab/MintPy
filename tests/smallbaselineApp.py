#!/usr/bin/env python3
# Author: Zhang Yunjun, Mar 2019
"""Test example datasets for the small baseline time series analysis workflow."""


import argparse
import os
import shutil
import subprocess
import sys
import tarfile
import time
from pathlib import Path

import mintpy.cli.view
from mintpy.objects.progress import FileProgressObject

CMAP_DICT = {
    'FernandinaSenDT128' : 'jet',
    'WellsEnvD2T399'     : 'jet_r',
    'KujuAlosAT422F650'  : 'jet_r',
}

URL_LIST = [
    'https://zenodo.org/record/5498198/files/FernandinaSenDT128.tar.xz',   # 280 MB
    'https://zenodo.org/record/6662079/files/SanFranSenDT42.tar.xz',       # 280 MB
    'https://zenodo.org/record/3952950/files/WellsEnvD2T399.tar.xz',       # 280 MB
    'https://zenodo.org/record/5502403/files/RidgecrestSenDT71.tar.xz',    # 480 MB
    'https://zenodo.org/record/6661536/files/WCapeSenAT29.tar.xz',         # 240 MB
    'https://zenodo.org/record/3952917/files/KujuAlosAT422F650.tar.xz',    # 230 MB
]

PROJ_NAMES = [os.path.basename(url).split('.tar.xz')[0] for url in URL_LIST]
TEMPLATE_FILES = [Path(__file__).resolve().parent / 'configs' / f'{proj_name}.txt' for proj_name in PROJ_NAMES]


#####################################################################################
DSET_INFO = """
  FernandinaSenDT128 for ISCE/topsStack
  SanFranSenDT42     for ARIA
  WellsEnvD2T399     for GAMMA
  RidgecrestSenDT71  for HyP3
  WCapeSenAT29       for SNAP
  KujuAlosAT422F650  for ROI_PAC
"""

EXAMPLE = """example:
  # regular tests
  $MINTPY_HOME/tests/smallbaselineApp.py

  # fast tests
  $MINTPY_HOME/tests/smallbaselineApp.py --nofresh
  $MINTPY_HOME/tests/smallbaselineApp.py --dset KujuAlosAT422F650

  # change the local test directory
  $MINTPY_HOME/tests/smallbaselineApp.py --dir ~/test

  # the most complete tests
  $MINTPY_HOME/tests/smallbaselineApp.py --test-pyaps --test-isce
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Test smallbaselineApp workflow with example datasets.'+DSET_INFO,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('--dset', dest='dset_name', nargs='+', metavar='DSET', choices=PROJ_NAMES, default=PROJ_NAMES,
                        help='name(s) of datasets to be tested.\n'+
                             f'Available datasets: {PROJ_NAMES}\n' +
                             '(default: all)')

    parser.add_argument('--test-pyaps', dest='test_pyaps', action='store_true',
                        help='Include testing of PyAPS (default: %(default)s).')

    parser.add_argument('--test-isce', dest='test_isce', action='store_true',
                        help='Include testing of ISCE-2/topsStack (default: %(default)s).')

    parser.add_argument('--dir', dest='test_dir', default='~/data/test',
                        help='test directory (default: %(default)s).')

    parser.add_argument('--nofresh', dest='fresh_start', action='store_false',
                        help='Use existing files WITHOUT starting from the scratch (default: %(default)s).')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # expand test_dir
    inps.test_dir = os.path.expanduser(inps.test_dir)
    inps.test_dir = os.path.expandvars(inps.test_dir)
    inps.test_dir = os.path.abspath(inps.test_dir)

    return inps


#####################################################################################
def test_smallbaselineApp(dset_name, test_dir, fresh_start=True, test_pyaps=False, test_isce=False):
    print('Go to test directory:', test_dir)
    os.chdir(test_dir)

    dset_idx = PROJ_NAMES.index(dset_name)
    dset_url = URL_LIST[dset_idx]
    template_file = TEMPLATE_FILES[dset_idx]

    # download tar file
    tar_fbase = os.path.splitext(os.path.basename(dset_url))[0]
    tar_files = [tar_fbase + x for x in ['.xz', '.gz'] if os.path.isfile(tar_fbase + x)]
    if tar_files:
        tar_file = tar_files[0]
        print('tar file exists, skip re-downloading.')
    else:
        tar_file = os.path.basename(dset_url)
        cmd = f'wget {dset_url}'
        print(f'downloading tar file ...\n{cmd}')
        os.system(cmd)

    # uncompress tar file
    if not fresh_start and os.path.isdir(dset_name):
        print(f'use existing files in project directory: {dset_name}')
    else:
        # remove existing directory
        if os.path.isdir(dset_name):
            print(f'remove existing project directory: {dset_name}')
            shutil.rmtree(dset_name)

        # uncompress tar file
        print(f'extracting content from tar file: {tar_file}')
        tar = tarfile.open(fileobj=FileProgressObject(tar_file))
        tar.extractall()
        tar.close()
        print('')

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

    # running smallbaselineApp
    # Note: execute script in command line instead of call main() for a clean run
    # to avoid strange error from prep_aria: not recognized as a supported file format.
    # which only occurs if all datasets are tested in one run
    cmd = f'smallbaselineApp.py {template_file}'
    print(cmd)
    status = subprocess.Popen(cmd, shell=True).wait()
    if status != 0:
        raise RuntimeError(f'Test failed for example dataset {dset_name}')

    # custom plot of velocity map
    vel_files = [os.path.join(work_dir, i) for i in ['geo/geo_velocity.h5', 'velocity.h5']]
    vel_file = [i for i in vel_files if os.path.isfile(i)][0]
    png_file = os.path.join(work_dir, 'pic', f'{os.path.splitext(os.path.basename(vel_file))[0]}.png')
    iargs = [vel_file, 'velocity', '--nodisplay', '--noverbose', '-o', png_file]
    if dset_name in CMAP_DICT.keys():
        iargs += ['-c', CMAP_DICT[dset_name]]
    mintpy.cli.view.main(iargs)

    # open final velocity map if on mac
    if sys.platform.lower().startswith('darwin'):
        cmd = f'open {png_file}'
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
        print(f'Start testing smallbaselineApp workflow on example dataset {i+1}/{num_dset}: {dset_name}')
        test_smallbaselineApp(dset_name,
                              test_dir=inps.test_dir,
                              fresh_start=inps.fresh_start,
                              test_pyaps=inps.test_pyaps,
                              test_isce=inps.test_isce)
        print('#'*100)
        print(f'   PASS testing of smallbaselineApp workflow on example dataset {i+1}/{num_dset}: {dset_name}')
        print('#'*100+'\n'*3)

    # print message
    if num_dset == len(PROJ_NAMES):
        m, s = divmod(time.time()-start_time, 60)
        msg  = '#'*50
        msg += '\n    PASS ALL testings without running errors.\n'
        msg += '#'*50
        msg += f'\nTotal time used: {m:02.0f} mins {s:02.1f} secs\n'
        print(msg)

    return


#####################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
