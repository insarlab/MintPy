############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Andre Theron, Zhang Yunjun, Jun 2019             #
# Email: andretheronsa@gmail.com                           #
############################################################


import os
from mintpy.utils import readfile, writefile


SPEED_OF_LIGHT = 299792458  # m / s


##################################################################################################
def write_rsc(atr, rsc_file):
    """Write to rsc file"""
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
