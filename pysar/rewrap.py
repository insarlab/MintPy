#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################


import sys

import h5py
import numpy as np


def usage():
    print('''usage: rewrap.py  ifgram_file   [output_name]

Re-wrap unwraped interferograms to wrapped interferograms.

example:
  rewrap.py  interferograms.h5
    ''')
    return


def rewrap(unw):
    rewrapped = unw - np.round(unw/(2*np.pi)) * 2*np.pi
    return rewrapped


def main(argv):
    try:
        file = argv[0]
    except:
        usage(); sys.exit(1)

    try:    outfile = argv[1]
    except: outfile = 'rewrapped_'+file

    print('writing >>> '+outfile)
    h5 = h5py.File(file, 'r')
    h5out = h5py.File(outfile,'w')
    gg = h5out.create_group('interferograms')
    ifgramList = list(h5['interferograms'].keys())
    print('number of interferograms: '+str(len(ifgramList)))
    for ifgram in ifgramList:
        print(ifgram)
        unw = h5['interferograms'][ifgram].get(ifgram)[:]
        rewrapped = rewrap(unw)
        group = gg.create_group(ifgram)
        dset = group.create_dataset(ifgram, data=rewrapped, compression='gzip')
        for key, value in h5['interferograms'][ifgram].attrs.items():
            group.attrs[key] = value

    h5.close()
    h5out.close()
    print('Done.')
    return outfile


if __name__ == '__main__':
    main(sys.argv[1:])


