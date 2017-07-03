#! /usr/bin/env python

import os
import sys
import h5py
import pysar._pysar_utilities as ut


def usage():
    print '''
*******************************************************************

  igram_closure.py Seeded_unwrapIfgram.h5

*******************************************************************
    '''

def main(argv):

    try:
        file=argv[0]
    except:
        usage(); sys.exit(1)
 
    h5file = h5py.File(file)
    curls,Triangles,C = ut.get_triangles(h5file)
 
    try:     curlfile = argv[1]
    except:  curlfile = 'curls.h5'
 
    if not os.path.isfile(curlfile):
        ut.generate_curls(curlfile,h5file,Triangles,curls)
    else:
        print curlfile + " already exists!"


if __name__ == '__main__':
    main(sys.argv[1:])


