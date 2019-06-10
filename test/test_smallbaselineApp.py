#!/usr/bin/env python3
# coding: utf-8
# Test example dataset for smallbaselineApp workflow
# Author: Zhang Yunjun, 2019-03-05


import os
import sys
import argparse
import subprocess
import numpy as np


# example datasets
pDict = {
    'KujuAlosAT422F650' : 'https://zenodo.org/record/2557863/files/KujuAlosAT422F650.tar.xz',
    'WellsEnvD2T399'    : 'https://zenodo.org/record/2613771/files/WellsEnvD2T399.tar.xz',
    'FernandinaSenDT128': 'https://zenodo.org/record/2596744/files/FernandinaSenDT128.tar.xz',
}
pNames = list(pDict.keys())


#####################################################################################
EXAMPLE = """example:
  $MINTPY_HOME/test/test_smallbaselineApp.py
  $MINTPY_HOME/test/test_smallbaselineApp.py  KujuAlosAT422F650
  $MINTPY_HOME/test/test_smallbaselineApp.py  --nofresh
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Test smallbaselineApp workflow with example datasets.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('--dset', dest='dset_name', choices=pNames+['all'], default='all',
                        help='name(s) of datasets to be tested.\n'+
                             'Available datasets: {}\n'.format(pNames)+
                             'Default is "all" for ALL DATASETS.')
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
        inps.dset_name = pNames
    if isinstance(inps.dset_name, str):
        inps.dset_name = [inps.dset_name]
    return inps


#####################################################################################
def run_dataset(dset_name, test_dir, fresh_start=True):
    print('Go to test directory:', test_dir)
    os.chdir(test_dir)

    # download tar file
    tar_file = '{}.tar.xz'.format(dset_name)
    if not os.path.isfile(tar_file):
        print('downloading tar file ...')
        cmd = 'wget {}'.format(pDict[dset_name])
        print(cmd)
        os.system(cmd)
    print('tar file exists, skip downloading.')

    # uncompress tar file
    if not fresh_start and os.path.isdir(dset_name):
        print('use existing files in {}'.format(dset_name))
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

    # runing smallbaselineApp
    work_dir = os.path.join(test_dir, dset_name, 'mintpy')
    os.chdir(work_dir)
    print('Go to work directory:', work_dir)

    cmd = 'rm ./inputs/ECMWF.h5'
    os.system(cmd)
    print('remove existing tropospheric delay file: ./inputs/ECMWF.h5')

    cmd = 'smallbaselineApp.py $MINTPY_HOME/docs/examples/input_files/{}.txt'.format(dset_name)
    cmd = os.path.expandvars(cmd)
    status = subprocess.Popen(cmd, shell=True).wait()
    if status is not 0:
        raise RuntimeError('Test failed for example dataset {}'.format(dset_name))
    return


#####################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    num_dset = len(inps.dset_name)
    for i in range(num_dset):
        dset_name = inps.dset_name[i]
        print('-'*50)
        print('Start testing smallbaselineApp workflow on exmaple dataset {}/{}: {}'.format(i+1, num_dset, dset_name))
        run_dataset(dset_name,
                    test_dir=inps.test_dir,
                    fresh_start=inps.fresh_start)
        print('Pass testing smallbaselineApp workflow on exmaple dataset {}/{}: {}'.format(i+1, num_dset, dset_name))

    if num_dset == len(pNames):
        print('-'*50)
        print('Pass ALL testings without runing errors.')
        print('-'*50)
    return


#####################################################################################
if __name__ == '__main__':
    main()
