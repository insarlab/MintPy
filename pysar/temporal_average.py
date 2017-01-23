#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2016, Yunjun Zhang                          #
# Author:  Yunjun Zhang                                    #
############################################################
# Modified from load_data.py written by Heresh Fattahi.
#

import os
import sys

import pysar._pysar_utilities as ut


#################################  Usage  ####################################
def usage():
    print '''
******************************************************************************
  Calculate temporal average/mean of multi-temporal datasets.

  Usage:
      temporal_average.py multi_temporal_file [output_filename]

  Example:
      temporal_average.py Coherence.h5 average_spatial_coherence.h5

******************************************************************************
    '''


#############################  Main Function  ################################
def main(argv):

    try: File = argv[0]
    except: usage(); sys.exit(1)

    try:    outName = argv[1]
    except: outName = os.path.splitext(File)[0]+'_tempAverage.h5'
    #except: outName = 'average_spatial_coherence.h5'

    print '\n*************** Average in Time Domain ******************'
    ut.temporal_average(File, outName)

##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
