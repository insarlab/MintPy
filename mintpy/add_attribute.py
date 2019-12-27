#! /usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2016                               #
############################################################


import os
import sys
from mintpy.utils import readfile, writefile, utils as ut
import mintpy.info as info


################################################################################
USAGE = """
usage: add_attribute.py file metadata_file
       add_attribute.py file key1=value1 [key2=value2 [...]]

Add/Update attributes to file.

Example:
  add_attribute.py timeseries.h5  unavco_attribute.txt
  add_attribute.py timeseries.h5  track=422   frame=650
  add_attribute.py ifgramStack.h5 ref_y=None  ref_x=None  #Use None value to delete attribute
"""


def usage():
    print(USAGE)
    return


def read_input_attribute(argv, print_msg=True):
    atr_new = dict()
    for i in range(1, len(argv)):
        if os.path.isfile(argv[i]):
            atr_tmp = readfile.read_template(argv[i])
            atr_new.update(atr_tmp)
        else:
            atr_tmp = argv[i].split('=')
            atr_new[atr_tmp[0].strip()] = atr_tmp[1].strip()

    if print_msg:
        print("The following attributes will be added/updated, or removed if new value is 'None':")
        info.print_attributes(atr_new)
    return atr_new


def update_file_attribute(fname, atr_new):
    # Read Original Attributes
    atr = readfile.read_attribute(fname)
    print('update {} file attribute: {}'.format(atr['FILE_TYPE'], fname))

    ext = os.path.splitext(fname)[1]
    if ext in ['.h5', '.he5']:
        fname = ut.add_attribute(fname, atr_new)
    else:
        if not ut.update_attribute_or_not(atr_new, atr):
            print('All updated (removed) attributes already exists (do not exists) and have the same value, skip update.')
        else:
            for key, value in iter(atr_new.items()):
                if value == 'None':
                    try:
                        atr.pop(key)
                    except:
                        pass
                else:
                    atr[key] = value

            rsc_file = '{}.rsc'.format(fname)
            print('writing >>> {}'.format(rsc_file))
            writefile.write_roipac_rsc(atr, out_file=rsc_file)
    return fname


def main(argv):
    # Check Inputs
    if not argv or argv[0] in ['-h', '--help']:
        usage()
        sys.exit(1)
    if len(argv) < 2 or not argv[1]:
        raise Exception('\nAt lease 2 inputs are needed.\n')
    infile = argv[0]

    atr_new = read_input_attribute(argv)

    update_file_attribute(fname=infile, atr_new=atr_new)

    return infile


################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
