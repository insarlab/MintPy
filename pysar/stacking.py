#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2017, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################


import sys
import argparse
import h5py
from pysar.utils import utils as ut
from pysar.objects import ifgramDatasetNames


#################################  Usage  ####################################
EXAMPLE='''example:
  stacking.py ifgramStack.h5 unwrapPhase -o averagePhaseVelocity.h5
  stacking.py ifgramStack.h5 coherence   -o averageSpatialCoherence.h5
'''

def createParser():
    parser = argparse.ArgumentParser(description='Stack multiple layers dataset into one.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs='+', help='File(s) to be stacked')
    parser.add_argument('-d','--dataset', dest='dataset_name', default=ifgramDatasetNames[0],\
                        help='Dataset to be used for stacking, when input file is ifgramStack')
    parser.add_argument('-m','--mask', dest='mask_file', help='Mask file for the calculation')
    parser.add_argument('-o','--output', dest='outfile', help='output file name')
    return parser

def cmdLineParse(iargs=None):
    parser = createParser()
    inps = parser.parse_args(args=iargs)
    return inps


#############################  Main Function  ################################
def main(iargs=None):
    inps = cmdLineParse(iargs)
    print('\n*************** Stacking ******************')
    for File in inps.file:
        ut.temporal_average(File, datasetName=inps.dataset_name, maskFile=inps.mask_file)
    return


##############################################################################
if __name__ == '__main__':
    main()
