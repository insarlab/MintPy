#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2016                               #
############################################################


import argparse
import matplotlib.pyplot as plt
from mintpy.utils import readfile, ptime, utils as ut, plot as pp


#################################  Usage  ####################################
EXAMPLE = """example:
  spatial_average.py inputs/ifgramStack.h5  -d coherence -m maskConnComp.h5
  spatial_average.py timeseries_ERA5_demErr.h5 -m maskTempCoh.h5
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Calculate average in space',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)
    parser.add_argument('file', help='File to calculate spatial average')
    parser.add_argument('-d', '--dset', '--dataset', dest='datasetName',
                        help='dataset used to calculate, for ifgramStack file only.')
    parser.add_argument('-m', '--mask', dest='mask_file',
                        help='Mask file for the calculation')
    parser.add_argument('--nodisplay', dest='disp_fig',
                        action='store_false', help='save and do not display the figure')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


#############################  Main Function  ################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    print('\n*************** Spatial Average ******************')
    mean_list, date_list = ut.spatial_average(inps.file,
                                              datasetName=inps.datasetName,
                                              maskFile=inps.mask_file,
                                              saveList=True)
    atr = readfile.read_attribute(inps.file)
    k = atr['FILE_TYPE']
    if inps.disp_fig and k == 'timeseries':
        dates, datevector = ptime.date_list2vector(date_list)
        # plot
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(dates, mean_list, '-o')#, lw=2, ms=16, alpha=0.7) #, mfc='crimson')
        ax.set_title('Spatial Average', fontsize=12)
        ax = pp.auto_adjust_xaxis_date(ax, datevector)[0]
        ax.set_xlabel('Time [years]', fontsize=12)
        ax.set_ylabel('Mean', fontsize=12)
        plt.show()
    return


##############################################################################
if __name__ == '__main__':
    main()
