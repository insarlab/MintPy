#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2016, Yunjun Zhang                          #
# Author:  Yunjun Zhang                                    #
############################################################
# Modified from load_data.py written by Heresh Fattahi.
#

import os
import sys

import _pysar_utilities as ut


#################################  Usage  ####################################
def usage():
    print('''usage:  temporal_average.py  file  [output_file]

Calculate temporal average/mean of multi-temporal datasets.

arguments:
  file        : string, file with multiple datasets/groups, e.g. timeseries.h5, coherence.h5, unwrapIfgram.h5
  output_file : string, path/name of output temporal average file

example:
  temporal_average.py  coherence.h5  averageSpatialCoherence.h5
    ''')
    return


#############################  Main Function  ################################
def main(argv):

    try:
        File = argv[0]
    except:
        usage(); sys.exit(1)

    try:    outName = argv[1]
    except: outName = os.path.splitext(File)[0]+'_tempAverage.h5'

    outName = ut.temporal_average(File, outName)
    print('Done.')
    return outName


##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
