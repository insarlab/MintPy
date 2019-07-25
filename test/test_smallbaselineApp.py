#!/usr/bin/env python3
# coding: utf-8
# Test example dataset for smallbaselineApp workflow
# Author: Zhang Yunjun, 2019-03-05


import os
import sys
import argparse
import subprocess
import numpy as np
from mintpy import smallbaselineApp


URL_LIST = [
    'https://zenodo.org/record/2748170/files/KujuAlosAT422F650.tar.xz',
    'https://zenodo.org/record/2748560/files/WellsEnvD2T399.tar.xz',
    'https://zenodo.org/record/2748487/files/FernandinaSenDT128.tar.xz',
]

PROJ_NAME_LIST = [os.path.basename(url).split('.')[0] for url in URL_LIST]
TEMPLATE_FILE_LIST = [os.path.join(os.path.dirname(__file__), '{}.txt'.format(proj_name))
                      for proj_name in PROJ_NAME_LIST]


#####################################################################################
EXAMPLE = """example:
  $MINTPY_HOME/test/test_smallbaselineApp.py
  $MINTPY_HOME/test/test_smallbaselineApp.py  --dset KujuAlosAT422F650
  $MINTPY_HOME/test/test_smallbaselineApp.py  --nofresh

  # available test datasets:
  {}
""".format(PROJ_NAME_LIST)

def create_parser():
    parser = argparse.ArgumentParser(description='Test smallbaselineApp workflow with example datasets.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('--dset', dest='dset_name', choices=PROJ_NAME_LIST+['all'], default='all',
                        help='name(s) of datasets to be tested.\n'+
                             'Available datasets: {}\n'.format(PROJ_NAME_LIST)+
                             'Default is "all" for ALL DATASETS.')
    parser.add_argument('--test-pyaps', dest='test_pyaps', action='store_true',
                        help='Include testing of PyAPS.')

    parser.add_argument('--dir', dest='test_dir', default='~/insarlab/test',
                        help='test directory. Default: ~/insarlab/test')

    parser.add_argument('--nofresh', dest='fresh_start', action='store_false',
                        help='Use exsiting files WITHOUT starting from the scratch.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    inps.test_dir = os.path.expanduser(inps.test_dir)

    if inps.dset_name.lower() == 'all':
        inps.dset_name = PROJ_NAME_LIST

    if isinstance(inps.dset_name, str):
        inps.dset_name = [inps.dset_name]
    return inps


#####################################################################################


def test_dataset(dset_name, test_dir, fresh_start=True, test_pyaps=False):
    print('Go to test directory:', test_dir)
    os.chdir(test_dir)

    dset_idx = PROJ_NAME_LIST.index(dset_name)
    dset_url = URL_LIST[dset_idx]
    template_file = TEMPLATE_FILE_LIST[dset_idx]

    # download tar file
    tar_file = os.path.basename(dset_url)
    if not os.path.isfile(tar_file):
        print('downloading tar file ...')
        cmd = 'wget {}'.format(dset_url)
        print(cmd)
        os.system(cmd)
    print('tar file exists, skip re-downloading.')

    # uncompress tar file
    if not fresh_start and os.path.isdir(dset_name):
        print('use existing files in project directory: {}'.format(dset_name))
    else:
        # remove existing directory
        if os.path.isdir(dset_name):
            cmd = 'rm -r {}'.format(dset_name)
            print('removing existing project directory')
            print(cmd)
            os.system(cmd)

        # uncompress tar file
        cmd = 'tar -xJf {}'.format(tar_file)
        print(cmd)
        os.system(cmd)

    # set working directory
    work_dir = os.path.join(test_dir, dset_name, 'mintpy')
    os.chdir(work_dir)
    print('Go to work directory:', work_dir)

    # remove pyaps existing products or not
    if test_pyaps:
        cmd = 'rm ./inputs/ECMWF.h5'
        os.system(cmd)
        print('remove existing tropospheric delay file: ./inputs/ECMWF.h5')

    # runing smallbaselineApp
    smallbaselineApp.main([template_file])
    return


#####################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    num_dset = len(inps.dset_name)
    for i in range(num_dset):
        dset_name = inps.dset_name[i]
        print('-'*50)
        print('Start testing smallbaselineApp workflow on exmaple dataset {}/{}: {}'.format(i+1, num_dset, dset_name))
        test_dataset(dset_name,
                     test_dir=inps.test_dir,
                     fresh_start=inps.fresh_start,
                     test_pyaps=inps.test_pyaps)
        print('PASS testing smallbaselineApp workflow on exmaple dataset {}/{}: {}'.format(i+1, num_dset, dset_name))

    if num_dset == len(PROJ_NAME_LIST):
        print('-'*50)
        print('PASS ALL testings without running errors.')
        print('-'*50)
    return


#####################################################################################
if __name__ == '__main__':
    main()
