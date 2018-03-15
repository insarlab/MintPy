#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2017, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################
# 

import sys
import argparse

import h5py

import pysar._pysar_utilities as ut


#################################  Usage  ####################################
EXAMPLE='''example:
  stacking.py unwrapIfgram.h5
  stacking.py coherence.h5 -m mask.h5
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Stack multiple layers dataset into one.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs='+', help='File(s) to be stacked')
    parser.add_argument('-m','--mask', dest='mask_file', help='Mask file for the calculation')

    inps = parser.parse_args()
    return inps


#############################  Main Function  ################################
def main(argv):
    inps = cmdLineParse()
    print '\n*************** Stacking ******************'
    for File in inps.file:
        ut.get_file_stack(File, inps.mask_file)
    return


##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
