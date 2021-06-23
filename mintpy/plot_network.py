#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from mintpy.objects import ifgramStack
from mintpy.utils import (readfile,
                          utils as ut,
                          plot as pp)
# suppress UserWarning from matplotlib
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")


###########################  Sub Function  #############################
BL_LIST = """
070106     0.0   0.03  0.0000000  0.00000000000 2155.2 /scratch/SLC/070106/
070709  2631.9   0.07  0.0000000  0.00000000000 2155.2 /scratch/SLC/070709/
070824  2787.3   0.07  0.0000000  0.00000000000 2155.2 /scratch/SLC/070824/
"""

DATE12_LIST = """
20070709_20100901
20070709_20101017
20070824_20071009
"""

EXAMPLE = """example:
  plot_network.py inputs/ifgramStack.h5
  plot_network.py inputs/ifgramStack.h5 -t smallbaselineApp.cfg --nodisplay   #Save figures to files without display
  plot_network.py inputs/ifgramStack.h5 -t smallbaselineApp.cfg --show-kept   #Do not plot dropped ifgrams
  plot_network.py coherenceSpatialAvg.txt

  # offsetSNR
  plot_network.py inputs/ifgramStack.h5 -d offsetSNR -v 0 20 --cmap-vlist 0 0.2 1
"""

TEMPLATE = """
mintpy.network.maskFile  = auto  #[file name, no], auto for waterMask.h5 or no for all pixels
mintpy.network.aoiYX     = auto  #[y0:y1,x0:x1 / no], auto for no, area of interest for coherence calculation
mintpy.network.aoiLALO   = auto  #[lat0:lat1,lon0:lon1 / no], auto for no - use the whole area
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Display Network of Interferograms',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file',
                        help='file with network information, ifgramStack.h5 or coherenceSpatialAvg.txt')
    parser.add_argument('--show-kept', dest='disp_drop', action='store_false',
                        help='display kept interferograms only, without dropped interferograms')
    parser.add_argument('-d', '--dset', type=str, dest='dsetName', default='coherence',
                        help='dataset used to calculate the mean')

    # Display coherence
    coh = parser.add_argument_group('Display Coherence', 'Show coherence of each interferogram pair with color')
    coh.add_argument('-t', '--template', dest='template_file',
                     help='template file with options below:\n'+TEMPLATE)
    coh.add_argument('-c', '--colormap', dest='cmap_name', default='RdBu_truncate',
                     help='colormap name for the network display. Default: RdBu_truncate')
    coh.add_argument('--cmap-vlist', dest='cmap_vlist', type=float, nargs=3, default=[0.2, 0.4, 1.0],
                     help='normalized start/jump/end value for truncated colormap. Default: 0.2 0.4 1.0')
    coh.add_argument('--mask', dest='maskFile', default='waterMask.h5',
                     help='mask file used to calculate the coherence. Default: waterMask.h5 or None.')

    coh.add_argument('-v', '--vlim', nargs=2, type=float, default=(0.2, 1.0),
                     help='display range')

    # Figure  Setting
    fig = parser.add_argument_group('Figure', 'Figure settings for display')
    fig.add_argument('--fs', '--fontsize', type=int,
                     default=12, help='font size in points')
    fig.add_argument('--lw', '--linewidth', dest='linewidth',
                     type=int, default=2, help='line width in points')
    fig.add_argument('--mc', '--markercolor', dest='markercolor',
                     default='orange', help='marker color')
    fig.add_argument('--ms', '--markersize', dest='markersize',
                     type=int, default=12, help='marker size in points')
    fig.add_argument('--every-year', dest='every_year', type=int,
                     default=1, help='number of years per major tick on x-axis')

    fig.add_argument('--dpi', dest='fig_dpi', type=int, default=150,
                     help='DPI - dot per inch - for display/write')
    fig.add_argument('--figsize', dest='fig_size', type=float, nargs=2,
                     help='figure size in inches - width and length')
    fig.add_argument('--notitle', dest='disp_title', action='store_false',
                     help='Do not display figure title.')
    fig.add_argument('--number', dest='number', type=str,
                     help='number mark to be plot at the corner of figure.')
    fig.add_argument('--nosplit-cmap', dest='split_cmap', action='store_false',
                     help='do not split colormap for coherence color')

    fig.add_argument('--save', dest='save_fig',
                     action='store_true', help='save the figure')
    fig.add_argument('--nodisplay', dest='disp_fig',
                     action='store_false', help='save and do not display the figure')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check input file type
    if inps.file.endswith(('.h5','.he5')):
        k = readfile.read_attribute(inps.file)['FILE_TYPE']
        if k != 'ifgramStack':
            raise ValueError('input HDF5 file is NOT ifgramStack.')

    if not inps.disp_fig:
        inps.save_fig = True

    if not inps.disp_fig:
        plt.switch_backend('Agg')

    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)
    else:
        inps.template = {}

    if not os.path.isfile(inps.maskFile):
        inps.maskFile = None

    return inps


def read_template2inps(template_file, inps=None):
    """Read input template options into Namespace inps"""
    if not inps:
        inps = cmd_line_parse()
    inpsDict = vars(inps)
    print('read options from template file: '+os.path.basename(template_file))
    inps.template = readfile.read_template(inps.template_file)
    inps.template = ut.check_template_auto_value(inps.template)

    # Coherence-based network modification
    prefix = 'mintpy.network.'
    key = prefix+'maskFile'
    if key in inps.template.keys():
        if inps.template[key]:
            inps.maskFile = inps.template[key]
    return inps


def read_network_info(inps):
    ext = os.path.splitext(inps.file)[1]
    print('read temporal/spatial baseline info from file:', inps.file)

    ## 1. Read date, pbase, date12 and coherence
    if ext == '.h5':
        inps.date12List = ifgramStack(inps.file).get_date12_list(dropIfgram=False)
        inps.dateList = ifgramStack(inps.file).get_date_list(dropIfgram=False)
        inps.pbaseList = ifgramStack(inps.file).get_perp_baseline_timeseries(dropIfgram=False)
        inps.cohList = ut.spatial_average(inps.file,
                                          datasetName=inps.dsetName,
                                          maskFile=inps.maskFile,
                                          saveList=True,
                                          checkAoi=False)[0]
    elif ext == '.txt':
        inps.date12List = np.loadtxt(inps.file, dtype=bytes).astype(str)[:,0].tolist()

        # date12List --> dateList
        mDates = [i.split('_')[0] for i in inps.date12List]
        sDates = [i.split('_')[1] for i in inps.date12List]
        inps.dateList = sorted(list(set(mDates + sDates)))

        # pbase12List + date12List --> pbaseList
        pbase12List = np.loadtxt(inps.file, dtype=bytes).astype(float)[:,3]
        A = ifgramStack.get_design_matrix4timeseries(inps.date12List)[0]
        inps.pbaseList = np.zeros(len(inps.dateList), dtype=np.float32)
        inps.pbaseList[1:] = np.linalg.lstsq(A, np.array(pbase12List), rcond=None)[0]

        # cohList
        inps.cohList = np.loadtxt(inps.file, dtype=bytes).astype(float)[:,1]
    else:
        raise ValueError('un-recognized input file extention:', ext)
    print('number of acquisitions: {}'.format(len(inps.dateList)))
    print('number of interferograms: {}'.format(len(inps.date12List)))

    print('shift all perp baseline by {} to zero mean for plotting'.format(np.mean(inps.pbaseList)))
    inps.pbaseList -= np.mean(inps.pbaseList)

    # Optional: Read dropped date12 / date
    inps.dateList_drop = []
    inps.date12List_drop = []
    if ext == '.h5':
        inps.date12List_keep = ifgramStack(inps.file).get_date12_list(dropIfgram=True)
        inps.date12List_drop = sorted(list(set(inps.date12List) - set(inps.date12List_keep)))
        print('-'*50)
        print('number of interferograms marked as drop: {}'.format(len(inps.date12List_drop)))
        print('number of interferograms marked as keep: {}'.format(len(inps.date12List_keep)))

        mDates = [i.split('_')[0] for i in inps.date12List_keep]
        sDates = [i.split('_')[1] for i in inps.date12List_keep]
        inps.dateList_keep = sorted(list(set(mDates + sDates)))
        inps.dateList_drop = sorted(list(set(inps.dateList) - set(inps.dateList_keep)))
        print('number of acquisitions marked as drop: {}'.format(len(inps.dateList_drop)))
        if len(inps.dateList_drop) > 0:
            print(inps.dateList_drop)

    return inps


def check_colormap(inps):
    """Check input colormaps."""
    # adjust color jump if no --cmap-vlist is input
    if inps.cohList is not None and '--cmap-vlist' not in sys.argv:
        # use template value
        key = 'mintpy.network.minCoherence'
        if key in inps.template.keys():
            inps.cmap_vlist[1] = max(0.01, float(inps.template[key]))

        # calculate from dropped interferograms, if there is any
        elif len(inps.date12List_drop) > 0:
            idx_drop = [inps.date12List.index(i) for i in inps.date12List_drop]
            coh_drop = [inps.cohList[i] for i in idx_drop]
            inps.cmap_vlist[1] = max(coh_drop)
            print('max coherence of excluded interferograms: {}'.format(max(coh_drop)))

        # update default cmap_vlist[0] value if the manually input cmap_vlist[1] is too small
        if inps.cmap_vlist[1] <= inps.cmap_vlist[0]:
            inps.cmap_vlist[0] = min(0, inps.cmap_vlist[1])

    # in case the manually input list is not in order
    inps.cmap_vlist = sorted(inps.cmap_vlist)
    inps.colormap = pp.ColormapExt(inps.cmap_name, vlist=inps.cmap_vlist).colormap

    return inps


##########################  Main Function  ##############################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    # read / calculate
    inps = read_network_info(inps)

    # Plot settings
    inps = check_colormap(inps)
    if inps.dsetName == 'coherence':
        inps.ds_name = 'Coherence'
        figNames = [i+'.pdf' for i in ['bperpHistory', 'coherenceMatrix', 'coherenceHistory', 'network']]
    elif inps.dsetName == 'offsetSNR':
        inps.ds_name = 'SNR'
        figNames = [i+'.pdf' for i in ['bperpHistory', 'SNRMatrix', 'SNRHistory', 'network']]
    inps.cbar_label = 'Average Spatial {}'.format(inps.ds_name)

    # Fig 1 - Baseline History
    fig, ax = plt.subplots(figsize=inps.fig_size)
    ax = pp.plot_perp_baseline_hist(ax,
                                    inps.dateList,
                                    inps.pbaseList,
                                    vars(inps),
                                    inps.dateList_drop)
    if inps.save_fig:
        fig.savefig(figNames[0], bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)
        print('save figure to {}'.format(figNames[0]))

    if inps.cohList is not None:
        # Fig 2 - Coherence Matrix
        fig, ax = plt.subplots(figsize=inps.fig_size)
        ax = pp.plot_coherence_matrix(ax,
                                      inps.date12List,
                                      inps.cohList,
                                      inps.date12List_drop,
                                      p_dict=vars(inps))[0]
        if inps.save_fig:
            fig.savefig(figNames[1], bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)
            print('save figure to {}'.format(figNames[1]))

        # Fig 3 - Min/Max Coherence History
        fig, ax = plt.subplots(figsize=inps.fig_size)
        ax = pp.plot_coherence_history(ax,
                                       inps.date12List,
                                       inps.cohList,
                                       p_dict=vars(inps))
        if inps.save_fig:
            fig.savefig(figNames[2], bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)
            print('save figure to {}'.format(figNames[2]))

    # Fig 4 - Interferogram Network
    fig, ax = plt.subplots(figsize=inps.fig_size)
    ax = pp.plot_network(ax,
                         inps.date12List,
                         inps.dateList,
                         inps.pbaseList,
                         vars(inps),
                         inps.date12List_drop)
    if inps.save_fig:
        fig.savefig(figNames[3], bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)
        print('save figure to {}'.format(figNames[3]))

    if inps.disp_fig:
        print('showing ...')
        plt.show()


############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
