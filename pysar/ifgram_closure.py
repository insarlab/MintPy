#!/usr/bin/env python3

import os
import sys
import h5py
import pysar.utils.utils as ut


############################################################
USAGE = """
usage: ifgram_closure.py  ifgram_file  [output_file]

Generate closure file for interferograms.

example:
  ifgram_closure.py  INPUTS/ifgramStack.h5
  ifgram_closure.py  INPUTS/ifgramStack.h5  curls.h5
"""

def usage():
    print(USAGE)
    return


############################################################
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
        print(curlfile + " already exists!")

    return curlfile


############################################################
if __name__ == '__main__':
    main(sys.argv[1:])

