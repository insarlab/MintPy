#!/usr/bin/env python3
# coding: utf-8
# Test example dataset for pysarApp workflow
# Author: Zhang Yunjun, 2019-03-05


import os
import sys
import subprocess
import numpy as np


test_dir = os.path.expanduser('~/insarlab/test')
proj_names = ['KujuAlosAT422F650', 'FernandinaSenDT128']
urls = ['https://zenodo.org/record/2557863/files/KujuAlosAT422F650.tar.xz',
        'https://zenodo.org/record/2571980/files/FernandinaSenDT128.tar.xz']
num_proj = len(proj_names)

# check input argument
if len(sys.argv) > 1:
    if sys.argv[1] in proj_names:
        print('Test on selected project:', sys.argv[1])
        idx = [proj_names.index(sys.argv[1])]
    else:
        msg = 'ERROR: input dataset not found: {}'.format(sys.argv[1])
        msg += '\navailable datasets: {}'.format(proj_names)
        raise SystemExit(msg)
else:
    idx = np.arange(num_proj)

for i in idx:
    proj_name = proj_names[i]
    print('-'*50)
    print('Start testing pysarApp workflow on exmaple dataset {}/{}: {}'.format(i+1, len(proj_names), proj_name))
    work_dir = os.path.join(test_dir, proj_name, 'PYSAR')

    # download tar file
    os.chdir(test_dir)
    print('Go to test directory:', test_dir)
    tar_file = '{}.tar.xz'.format(proj_name)
    if not os.path.isfile(tar_file):
        print('downloading tar file ...')
        cmd = 'wget {}'.format(urls[i])
        print(cmd)
        os.system(cmd)
    print('tar file exists, skip downloading.')

    # uncompress tar file
    if os.path.isdir(proj_name):
        cmd = 'rm -r {}'.format(proj_name)
        print('removing existing project directory')
        print(cmd)
        os.system(cmd)
    cmd = 'tar -xJf {}'.format(tar_file)
    print(cmd)
    os.system(cmd)

    # runing pysarApp
    os.chdir(work_dir)
    print('Go to work directory:', work_dir)

    cmd = 'pysarApp.py {}.txt'.format(proj_name)
    status = subprocess.Popen(cmd, shell=True).wait()
    if status is not 0:
        raise RuntimeError('Test failed for example dataset {}'.format(proj_name))
    print('Test passed for exmaple dataset {}/{}: {}'.format(i+1, len(proj_names), proj_name))
print('-'*50)
print('Passed all testings without runing errors.')
print('-'*50)

