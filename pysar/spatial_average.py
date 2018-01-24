#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2016, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################
# 

import sys
import argparse

import matplotlib.pyplot as plt

import _readfile as readfile
import _pysar_utilities as ut
import _datetime as ptime


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
    print('\n*************** Spatial Average ******************')
    for File in inps.file:
        mean_list, date_list = ut.spatial_average(File, inps.mask_file, saveList=True)
        atr = readfile.read_attribute(File)
        k = atr['FILE_TYPE']
        if inps.disp_fig and k == 'timeseries':
            dates, datevector = ptime.date_list2vector(date_list)

            # plot
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(dates, mean_list, '-ko', lw=2, ms=16, alpha=0.7, mfc='crimson')
            ax.set_title('Spatial Average',fontsize=12)
            ax = ptime.auto_adjust_xaxis_date(ax, datevector)[0]
            ax.set_xlabel('Time [years]',fontsize=12)
            ax.set_ylabel('Mean',fontsize=12)
            plt.show()


##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
