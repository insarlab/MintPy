#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2016, Yunjun Zhang                          #
# Author:  Yunjun Zhang                                    #
############################################################
# 

import sys
import argparse

import h5py
import matplotlib.pyplot as plt

import pysar._readfile as readfile
import pysar._pysar_utilities as ut
import pysar._datetime as ptime
from pysar._readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file


#################################  Usage  ####################################
EXAMPLE='''example:
  spatial_average.py coherence.h5
  spatial_average.py unwrapIfgram.h5 -m Mask.h5
  spatial_average.py sum_timeseries_ECMWF_demCor.h5 -m Mask_tempCoh.h5
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Calculate average in space',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)
    
    parser.add_argument('file', nargs='+', help='File(s) to calculate spatial average')
    parser.add_argument('-m','--mask', dest='mask_file', help='Mask file for the calculation')
    parser.add_argument('--nodisplay', dest='disp_fig', action='store_false', help='save and do not display the figure')

    inps = parser.parse_args()
    return inps


#############################  Main Function  ################################
def main(argv):
    inps = cmdLineParse()
    print '\n*************** Spatial Average ******************'

    if inps.mask_file:
        print 'reading mask file: '+inps.mask_file
        mask, mask_atr = readfile.read(inps.mask_file)
    else:
        mask = None
    
    for File in inps.file:
        mean_list = ut.spatial_average(File, mask, saveList=True)
        atr = readfile.read_attribute(File)
        k = atr['FILE_TYPE']
        if inps.disp_fig and k == 'timeseries':
            # Get date list
            h5file = h5py.File(File)
            dateList = sorted(h5file[k].keys())
            h5file.close()
            dates, datevector = ptime.date_list2vector(dateList)

            # plot
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(dates, mean_list, '-ko', lw=2, ms=16, alpha=0.7, mfc='crimson')
            ax.set_title('Spatial Average',fontsize=12)
            ax = ptime.auto_adjust_xaxis_date(ax, datevector)
            ax.set_xlabel('Time [years]',fontsize=12)
            ax.set_ylabel('Mean',fontsize=12)
            plt.show()


##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
