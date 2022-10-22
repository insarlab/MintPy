#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Antonio Valentino, Aug 2016        #
############################################################


import sys

from mintpy.utils.arg_utils import create_argument_parser

#################################  Usage  ####################################
EXAMPLE = """example:
  spatial_average.py inputs/ifgramStack.h5  -d coherence -m maskConnComp.h5
  spatial_average.py timeseries_ERA5_demErr.h5 -m maskTempCoh.h5
"""

def create_parser(subparsers=None):
    synopsis = 'Calculate the spatial average.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

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


##############################  Sub Function  ################################
def plot_spatial_average_ts(date_list:list, avg_list:list):
    """Plot the spatial average time-series."""
    import matplotlib.pyplot as plt

    from mintpy.utils import plot as pp, ptime

    dates = ptime.date_list2vector(date_list)[0]

    # plot
    _, ax = plt.subplots()
    ax.plot(dates, avg_list, '-o')

    # axis format
    ax = pp.auto_adjust_xaxis_date(ax, dates)[0]
    ax.set_title('Spatial Average', fontsize=12)
    ax.set_xlabel('Time [years]', fontsize=12)
    ax.set_ylabel('Mean', fontsize=12)
    plt.show()

    return


#############################  Main Function  ################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.utils import readfile, utils1 as ut

    # run
    avg_list, date_list = ut.spatial_average(
        inps.file,
        datasetName=inps.datasetName,
        maskFile=inps.mask_file,
        saveList=True,
    )

    ftype = readfile.read_attribute(inps.file)['FILE_TYPE']
    if inps.disp_fig and ftype == 'timeseries':
        plot_spatial_average_ts(date_list, avg_list)

    return


##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
