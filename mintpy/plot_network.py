############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import numpy as np

from mintpy.objects import ifgramStack
from mintpy.utils import readfile, utils as ut, plot as pp

# suppress UserWarning from matplotlib
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")


###########################  Sub Function  #############################
def read_template2inps(template_file, inps):
    """Read input template options into Namespace inps"""
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
            print('max coherence of excluded interferograms: {}'.format(max(coh_drop)))

        # update default cmap_vlist[0] value if the manually input cmap_vlist[1] is too small
        if inps.cmap_vlist[1] <= inps.cmap_vlist[0]:
            inps.cmap_vlist[0] = min(0, inps.cmap_vlist[1])

    # in case the manually input list is not in order
    inps.cmap_vlist = sorted(inps.cmap_vlist)
    inps.colormap = pp.ColormapExt(inps.cmap_name, vlist=inps.cmap_vlist).colormap

    return inps
