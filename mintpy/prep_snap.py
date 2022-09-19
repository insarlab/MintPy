############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Andre Theron, Zhang Yunjun, Jun 2019             #
# Email: andretheronsa@gmail.com                           #
############################################################


import os
from mintpy.utils import readfile, writefile, utils1 as ut


##################################################################################################
def prep_snap_metadata(snap_file):
    """Prepare SNAP metadata files.

    Parameters: snap_file - str, SNAP data file with .img file extension.
    Returns:    rsc_file  - str, metadata file
    """

    # read metadata from *.dim file
    # the map info from *.img.hdr file is NOT right, thus, not used.
    dim_file = os.path.dirname(snap_file)[:-4] + 'dim'
    atr = readfile.read_snap_dim(dim_file)

    rsc_file = snap_file + '.rsc'

    # grab atr from existing rsc file
    if os.path.isfile(rsc_file):
        atr_orig = readfile.read_roipac_rsc(rsc_file)
    else:
        atr_orig = dict()

    # (over)write to rsc file if input atr has more items
    if not set(atr.items()).issubset(set(atr_orig.items())):
        atr_out = {**atr_orig, **atr}
        print('write metadata to {} '.format(os.path.basename(rsc_file)))
        writefile.write_roipac_rsc(atr_out, out_file=rsc_file)

    return rsc_file


def run_prep_snap(inps):
    # grab all input files
    inps.file = ut.get_file_list(inps.file, abspath=True)

    # check: input file extensions
    for fname in inps.file:
        if not fname.endswith('.img'):
            raise ValueError('Input data file does NOT end with .img: {}'.format(fname))

    # loop over each file
    for img_file in inps.file:
        prep_snap_metadata(img_file)

    return
