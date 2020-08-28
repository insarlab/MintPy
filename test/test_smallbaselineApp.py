#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2019-03-05                         #
############################################################
# Test example dataset for smallbaselineApp workflow


import os
import sys
import time
import shutil
import argparse
import subprocess
import tarfile


CMAP_DICT = {
    'FernandinaSenDT128' : 'jet',
    'WellsEnvD2T399'     : 'jet_r',
    'KujuAlosAT422F650'  : 'jet_r',
}

URL_LIST = [
    'https://zenodo.org/record/3952953/files/FernandinaSenDT128.tar.xz',
    'https://zenodo.org/record/3952950/files/WellsEnvD2T399.tar.xz',
    'https://zenodo.org/record/3952917/files/KujuAlosAT422F650.tar.xz',
]

PROJ_NAME_LIST = [os.path.basename(url).split('.tar.xz')[0] for url in URL_LIST]
TEMPLATE_FILE_LIST = [os.path.join(os.path.dirname(__file__), 'configs/{}.txt'.format(proj_name))
                      for proj_name in PROJ_NAME_LIST]


#####################################################################################
EXAMPLE = """example:
  $MINTPY_HOME/test/test_smallbaselineApp.py
  $MINTPY_HOME/test/test_smallbaselineApp.py  --dir ~/test
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

    parser.add_argument('--dir', dest='test_dir', default='~/data/test',
                        help='test directory. Default: ~/data/test')

    parser.add_argument('--nofresh', dest='fresh_start', action='store_false',
                        help='Use exsiting files WITHOUT starting from the scratch.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # expand test_dir
    inps.test_dir = os.path.expanduser(inps.test_dir)
    inps.test_dir = os.path.expandvars(inps.test_dir)

    # translate --dset all
    if inps.dset_name.lower() == 'all':
        inps.dset_name = PROJ_NAME_LIST

    elif isinstance(inps.dset_name, str):
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

    # remove pyaps existing products or not
    if test_pyaps:
        print('remove existing tropospheric delay file: ./inputs/ERA5.h5')
        os.remove('./inputs/ERA5.h5')

    # runing smallbaselineApp
    cmd = 'smallbaselineApp.py {}'.format(template_file)
    print(cmd)
    status = subprocess.Popen(cmd, shell=True).wait()
    if status != 0:
        raise RuntimeError('Test failed for example dataset {}'.format(dset_name))

    # custom plot of velocity map
    if dset_name in CMAP_DICT.keys():
        cmd = 'view.py geo/geo_velocity.h5 velocity --nodisplay -o pic/geo_velocity.png --noverbose '
        cmd += ' -c {}'.format(CMAP_DICT[dset_name])
        print(cmd)
        subprocess.Popen(cmd, shell=True).wait()

    # open final velocity map if on mac
    if sys.platform.lower().startswith('darwin'):
        cmd = 'open pic/geo_velocity.png'
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
        print('-'*50)
        print('Start testing smallbaselineApp workflow on exmaple dataset {}/{}: {}'.format(i+1, num_dset, dset_name))
        test_dataset(dset_name,
                     test_dir=inps.test_dir,
                     fresh_start=inps.fresh_start,
                     test_pyaps=inps.test_pyaps)
        print('PASS testing smallbaselineApp workflow on exmaple dataset {}/{}: {}'.format(i+1, num_dset, dset_name))

    # print message
    if num_dset == len(PROJ_NAME_LIST):
        m, s = divmod(time.time()-start_time, 60)
        msg = '-'*50
        msg += '\nPASS ALL testings without running errors.'
        msg += '\n'+'-'*50
        msg += '\nTotal time used: {:02.0f} mins {:02.1f} secs\n'.format(m, s)
        print(msg)

    return


#####################################################################################
if __name__ == '__main__':
    main()
