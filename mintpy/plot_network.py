############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import warnings

import numpy as np
from matplotlib import pyplot as plt

from mintpy.objects import ifgramStack
from mintpy.utils import plot as pp, readfile, utils as ut

# suppress UserWarning from matplotlib
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")


###########################  Sub Function  #############################
def read_network_info(inps):
    ext = os.path.splitext(inps.file)[1]
    print('read temporal/spatial baseline info from file:', inps.file)

    ## 1. Read date, pbase, date12 and coherence
    if ext == '.h5':
        stack_obj = ifgramStack(inps.file)
        stack_obj.open()
        inps.date12List = stack_obj.get_date12_list(dropIfgram=False)
        inps.dateList = stack_obj.get_date_list(dropIfgram=False)
        inps.pbaseList = stack_obj.get_perp_baseline_timeseries(dropIfgram=False)

        if inps.dsetName in readfile.get_dataset_list(inps.file):
            inps.cohList = ut.spatial_average(inps.file,
                                              datasetName=inps.dsetName,
                                              maskFile=inps.maskFile,
                                              saveList=True,
                                              checkAoi=False)[0]
        elif inps.dsetName == 'pbase':
            inps.cohList = np.abs(stack_obj.pbaseIfgram).tolist()

        elif inps.dsetName == 'tbase':
            inps.cohList = stack_obj.tbaseIfgram.tolist()

        else:
            inps.cohList = None
            msg = f'{inps.dsetName} NOT found in file: {inps.file}! '
            msg += 'Disable the color coding and continue'
            print(msg)

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
    print(f'number of acquisitions: {len(inps.dateList)}')
    print(f'number of interferograms: {len(inps.date12List)}')

    print(f'shift all perp baseline by {np.mean(inps.pbaseList)} to zero mean for plotting')
    inps.pbaseList -= np.mean(inps.pbaseList)

    # Optional: Read dropped date12 / date
    inps.dateList_drop = []
    inps.date12List_drop = []
    if ext == '.h5':
        inps.date12List_keep = ifgramStack(inps.file).get_date12_list(dropIfgram=True)
        inps.date12List_drop = sorted(list(set(inps.date12List) - set(inps.date12List_keep)))
        print('-'*50)
        print(f'number of interferograms marked as drop: {len(inps.date12List_drop)}')
        print(f'number of interferograms marked as keep: {len(inps.date12List_keep)}')

        mDates = [i.split('_')[0] for i in inps.date12List_keep]
        sDates = [i.split('_')[1] for i in inps.date12List_keep]
        inps.dateList_keep = sorted(list(set(mDates + sDates)))
        inps.dateList_drop = sorted(list(set(inps.dateList) - set(inps.dateList_keep)))
        print(f'number of acquisitions marked as drop: {len(inps.dateList_drop)}')
        if len(inps.dateList_drop) > 0:
            print(inps.dateList_drop)

    return inps


def check_colormap(inps):
    """Check input colormaps."""
    # adjust color jump if no --cmap-vlist is input
    if inps.cohList is not None and '--cmap-vlist' not in inps.argv:
        # use template value
        key = 'mintpy.network.minCoherence'
        if key in inps.template.keys():
            inps.cmap_vlist[1] = max(0.01, float(inps.template[key]))

        # calculate from dropped interferograms, if there is any
        elif len(inps.date12List_drop) > 0:
            idx_drop = [inps.date12List.index(i) for i in inps.date12List_drop]
            coh_drop = [inps.cohList[i] for i in idx_drop]
            inps.cmap_vlist[1] = max(coh_drop)
            print(f'max coherence of excluded interferograms: {max(coh_drop)}')

        # update default cmap_vlist[0] value if the manually input cmap_vlist[1] is too small
        if inps.cmap_vlist[1] <= inps.cmap_vlist[0]:
            inps.cmap_vlist[0] = min(0, inps.cmap_vlist[1])

    # in case the manually input list is not in order
    inps.cmap_vlist = sorted(inps.cmap_vlist)
    inps.colormap = pp.ColormapExt(inps.cmap_name, vlist=inps.cmap_vlist).colormap

    return inps


########################################################################
def plot_network(inps):
    """Plot all the network info."""

    # matplotlib backend setting
    if not inps.disp_fig:
        plt.switch_backend('Agg')

    # read / calculate
    inps = read_network_info(inps)

    # plot settings
    inps = check_colormap(inps)
    ext = '.pdf'
    if os.path.basename(inps.file).startswith('ion'):
        ext = f'_ion{ext}'
    kwargs = dict(bbox_inches='tight', transparent=True, dpi=inps.fig_dpi)

    if inps.dsetName == 'coherence':
        fig_names = [i+ext for i in ['pbaseHistory', 'coherenceHistory', 'coherenceMatrix', 'network']]
        inps.ds_name = 'Coherence'
        inps.cbar_label = 'Average Spatial Coherence'

    elif inps.dsetName == 'offsetSNR':
        fig_names = [i+ext for i in ['pbaseHistory', 'SNRHistory', 'SNRMatrix', 'network']]
        inps.ds_name = 'SNR'
        inps.cbar_label = 'Average Spatial SNR'

    elif inps.dsetName == 'tbase':
        fig_names = [i+ext for i in ['pbaseHistory', 'tbaseHistory', 'tbaseMatrix', 'network']]
        inps.ds_name = 'Temporal Baseline'
        inps.cbar_label = 'Temporal Baseline [day]'

    elif inps.dsetName == 'pbase':
        fig_names = [i+ext for i in ['pbaseHistory', 'pbaseRangeHistory', 'pbaseMatrix', 'network']]
        inps.ds_name = 'Perp Baseline'
        inps.cbar_label = 'Perp Baseline [m]'

    # Fig 1 - Baseline History
    fig, ax = plt.subplots(figsize=inps.fig_size)
    ax = pp.plot_perp_baseline_hist(
        ax,
        inps.dateList,
        inps.pbaseList,
        vars(inps),
        inps.dateList_drop,
    )
    if inps.save_fig:
        fig.savefig(fig_names[0], **kwargs)
        print(f'save figure to {fig_names[0]}')

    if inps.cohList is not None:
        # Fig 2 - Min/Max Coherence History
        fig, ax = plt.subplots(figsize=inps.fig_size)
        ax = pp.plot_coherence_history(
            ax,
            inps.date12List,
            inps.cohList,
            p_dict=vars(inps),
        )
        if inps.save_fig:
            fig.savefig(fig_names[1], **kwargs)
            print(f'save figure to {fig_names[2]}')

        # Fig 3 - Coherence Matrix
        fig, ax = plt.subplots(figsize=inps.fig_size)
        ax = pp.plot_coherence_matrix(
            ax,
            inps.date12List,
            inps.cohList,
            inps.date12List_drop,
            p_dict=vars(inps),
        )[0]
        if inps.save_fig:
            fig.savefig(fig_names[2], **kwargs)
            print(f'save figure to {fig_names[1]}')

    # Fig 4 - Interferogram Network
    fig, ax = plt.subplots(figsize=inps.fig_size)
    ax = pp.plot_network(
        ax,
        inps.date12List,
        inps.dateList,
        inps.pbaseList,
        vars(inps),
        inps.date12List_drop,
    )
    if inps.save_fig:
        fig.savefig(fig_names[3], **kwargs)
        print(f'save figure to {fig_names[3]}')

    if inps.disp_fig:
        print('showing ...')
        plt.show()
    else:
        plt.close()

    return
