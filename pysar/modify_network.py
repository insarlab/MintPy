#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# Yunjun, Jul 2015: Update Coherence file
# Yunjun, Oct 2015: Add 'T' option for template file input
#                   and pysar.dropIfgIndex in template content
# Yunjun, Jan 2017: add spatial_average(), update_inps_from_template()
#                   modify_file_date12_list(), cmdLineParse()
#                   merge update_network() into this.
#                   add coherence-based network modification.


import os
import sys
import argparse
import getopt
import time
import datetime

import h5py
import numpy as np
import matplotlib.pyplot as plt

import pysar._datetime as ptime
import pysar._network as pnet
import pysar._pysar_utilities as ut
import pysar._readfile as readfile
from pysar._readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file


###########################  Sub Function  #############################
######################################
def nearest_neighbor(x,y, tbase, pbase):
    """ find nearest neighbour
    Input:
        tbase/pbase : numpy.array, temporal/perpendicular spatial baseline
    """
    dist = sqrt((tbase -x)**2+(pbase -y)**2)
    indx = dist==min(dist)
    return indx


def spatial_average(File, mask=None):
    '''Calculate spatial average value for multiple dataset/group file
    Example:
        meanList = spatial_average('coherence.h5')
        meanList = spatial_average('timeseries_ECMWF_demCor_quadratic.h5', mask)
    '''
    atr = readfile.read_attribute(File)
    k = atr['FILE_TYPE']
    # Read Epoch List
    h5file = h5py.File(inps.file,'r')
    epochList = h5file[k].keys()
    epochList = sorted(epochList)
    
    meanList = np.zeros(len(epochList))
    for epoch in epochList:
        # Read data
        if k in multi_group_hdf5_file:
            data = h5file[k][epoch].get(epoch)[:]
        elif k in multi_dataset_hdf5_file:
            data = h5file[k].get(epoch)[:]
        else:
            print 'Un-supported data type: '+k
            data = 0
        # Mask
        if not mask:
            data[mask==0] = np.nan
        # Calculate average in space
        meanList[epochList.index(epoch)] = np.nanmean(data)
    h5file.close()
    
    return meanList
    


def update_inps_with_template(inps, template_file):
    template_dict = readfile.read_template(inps.template_file)
    keyList = template_dict.keys()
    
    if not inps.max_temp_baseline and 'pysar.network.maxTempBaseline' in keyList:
        inps.max_temp_baseline = float(template_dict['pysar.network.maxTempBaseline'])
    if not inps.max_perp_baseline and 'pysar.network.maxPerpBaseline' in keyList:
        inps.max_perp_baseline = float(template_dict['pysar.network.maxPerpBaseline'])
    
    if not inps.drop_date and 'pysar.network.dropDate' in keyList:
        inps.drop_date = [i for i in template_dict['pysar.network.dropDate'].replace(',',' ').split()]
    if not inps.drop_ifg_index and 'pysar.network.dropIfgramIndex' in keyList:
        inps.drop_ifg_index = [i for i in template_dict['pysar.network.dropDate'].replace(',',' ').split()]
    
    if not inps.reference_file and 'pysar.network.reference' in keyList:
        inps.reference_file = template_dict['pysar.network.reference']
    
    if not inps.coherence_file and 'pysar.network.coherenceBase.coherence' in keyList:
        inps.coherence_file = template_dict['pysar.network.coherenceBase.coherence']
    if not inps.mask_file and 'pysar.network.coherenceBase.mask' in keyList:
        inps.mask_file = template_dict['pysar.network.coherenceBase.mask']
    if not inps.min_coherence and 'pysar.network.coherenceBase.minCoherence' in keyList:
        inps.min_coherence = float(template_dict['pysar.network.coherenceBase.minCoherence'])
    
    return inps


def modify_file_date12_list(File, date12_to_rmv, outFile=None):
    '''Update multiple group hdf5 file using date12 to remove/keep'''
    date12_orig = pnet.get_date12_list(File)
    
    if not outFile:
        outFile = 'Modified_'+os.path.basename(File)
    print 'writing >>> '+outFile
    h5out = h5py.File(outFile, 'w')
    gg = h5out.create_group(k)
    
    k = readfile.read_attribute(File)['FILE_TYPE']
    h5 = h5py.File(File, 'r')
    igramList = h5[k].keys()
    for igram in igramList:
        date12 = h5[k][igram].attrs['DATE12']
        if not date12 in date12_to_rmv:
            print igram
            data = h5[k][igram].get(igram)[:]
            
            group = gg.create_group(igram)
            dset = group.create_dataset(igram, data=data, compression='gzip')
            for key, value in h5file[k][igram].attrs.iteritems():
                group.attrs[key] = value
    h5.close()
    h5out.close()
    
    return outFile



##############  Usage  ###############

EXAMPLE='''example:
  modify_network.py LoadedData.h5
  modify_network.py -f Seeded_LoadedData.h5 -T $TE/KirishimaT246EnvD2.template
  modify_network.py -f LoadedData.h5 -b 400 -n no
  modify_network.py -f LoadedData.h5 -b 500 -t 600 -d 080307
  modify_network.py -f LoadedData.h5 -d '080307 091023'           
  modify_network.py -f LoadedData.h5 -N '1 4 20 76 89 100'
  modify_network.py -f LoadedData.h5 -C Coherence.h5 -N '1 4 20 76 89 100'

************************************************************************************
'''

TEMPLATE='''template:
  pysar.network.dropIfgramIndex = 7:9 15 25 26
  pysar.network.dropDate        = 20080520 20090816
  
  pysar.network.maxTempBaseline = 720
  pysar.network.maxPerpBaseline = 2000
  
  pysar.network.reference   = Modified_unwrapIfgram.h5
  pysar.network.reference   = Paris.list
  
  pysar.network.coherenceBase.coherence    =  coherence.h5
  pysar.network.coherenceBase.mask         =  maskAmp.h5
  pysar.network.coherenceBase.minCoherence =  0.7
'''


def cmdLineParse():
    parser = argparse.ArgumentParser(description='Modify the network of interferograms',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)
    parser.add_argument('file', nargs='+',\
                        help='Files to modify/drop network.\n'\
                             'i.e. unwrapIfgram.h5, wrapIfgram.h5, coherence.h5, ...')
    parser.add_argument('-r','--reference', dest='reference_file',\
                        help='Reference hdf5 / list file with network information.\n'\
                             'i.e. Modified_unwrapIfgram.h5, Pairs.list')
    parser.add_argument('--template', dest='template_file', help='Template file with ')
    
    parser.add_argument('-t', dest='max_temp_baseline', type=float, help='temporal baseline threshold/maximum in days')
    parser.add_argument('-b', dest='max_perp_baseline', type=float, help='perpendicular baseline threshold/maximum in meters')
    parser.add_argument('--drop-ifg-index', dest='drop_ifg_index', type=int, nargs='*',\
                        help='index of interferograms to remove/drop.\n1 as the first')
    parser.add_argument('--drop-date', dest='drop_date', nargs='*',\
                        help='date(s) to remove/drop, all interferograms included date(s) will be removed')
    
    # Coherence-based network    
    coherenceGroup = parser.add_argument_group('Coherence-based Network',\
                                               'Drop/modify network based on spatial coherence')
    coherenceGroup.add_argument('--coherence-base', dest='coherence_file',\
                                help='Enable coherence-based network modification\n'\
                                     'Input coherence file should have the same network as input file(s)')
    coherenceGroup.add_argument('--mask', dest='mask_file',\
                                help='Mask file used to calculate the spatial coherence\n'\
                                     'Will use the whole area if not assigned')
    coherenceGroup.add_argument('--min-coherence', dest='min_coherence', type=float, default=0.7,\
                                help='Minimum coherence value')
    
    # Manually select network
    manualGroup = parser.add_argument_group('Manual Network', 'Manually select/drop/modify network')
    manualGroup.add_argument('--disp-network', dest='disp_network', action='store_true',\
                             help='display network to manually choose line/interferogram to remove')

    inps = parser.parse_args()
    return inps

#########################  Main Function  ##############################
def main(argv):
    ##### Read Inputs
    inps = cmdLineParse()
    if (not inps.reference_file and not inps.template_file and\
        not inps.max_temp_baseline and not inps.max_perp_baseline and\
        not inps.drop_ifg_index and not inps.drop_date and \
        not inps.coherence_file):
        # Display network for manually modification when there is no other modification input.
        inps.disp_network = True
    
    # Update inps if template is input
    if not inps.template_file:
        inps = update_inps_with_template(inps, inps.template_file)
    
    # Convert index : input to continous index list
    if not inps.drop_ifg_index:
        drop_ifg_index = list(inps.drop_ifg_index)
        inps.drop_ifg_index = []
        for index in drop_ifg_index:
            index_temp = [int(i) for i in index.split(':')]
            index_temp.sort()
            if len(index_temp)==2:
                for j in range(index_temp[0], index_temp[1]+1):
                    inps.drop_ifg_index.append(str(j))
            elif len(index_temp)==1:
                inps.drop_ifg_index.append(index)
            else:
                print 'Unrecoganized input: '+index
        inps.drop_ifg_index = sorted(inps.drop_ifg_index)

    ##### Basic Info of Input File
    date8_list = ptime.date_list(inps.file[0])
    tbase_list = ptime.date_list2tbase(date8_list)[0]
    date6_list = ptime.yymmdd(date8_list)
    ifg_bperp_list = pnet.igram_perp_baseline_list(inps.file[0])
    
    ##### Convert 
    date12_orig = pnet.get_date12_list(inps.file[0])
    date12_to_rmv = []
    
    # 1. Update date12_to_keep from reference file
    if not inps.reference_file:
        date12_to_keep = pnet.get_date12_list(inps.reference_file)
        for date12 in date12_orig:
            if date12 not in date12_to_keep:
                date12_to_rmv.append(date12)
    
    # 2.1 Update date12_to_rmv from coherence file
    if not inps.coherence_file:
        # Calculate spatial average coherence
        if not inps.mask_file:
            mask = readfile.read(inps.mask_file)[0]
        else:
            mask = None
        mean_coherence_list = spatial_average(inps.coherence_file, mask)
        
        coh_date12_list = pnet.get_date12_list(inps.coherence_file)
        for i in range(len(coh_date12_list)):
            if mean_coherence_list[i] < inps.min_coherence:
                date12_to_rmv.append(coh_date12_list[i])

    # 2.2 Update date12_to_rmv from perp baseline threshold
    if not inps.max_perp_baseline:
        print 'Drop pairs with perpendicular spatial baseline > '+str(inps.max_perp_baseline)+' meters'
        for i in range(len(ifg_bperp_list)):
            if ifg_bperp_list[i] > inps.max_perp_baseline:
                date12_to_rmv.append(date12_orig[i])
    
    # 2.3 Update date12_to_rmv from temp baseline threshold
    if not inps.max_temp_baseline:
        print 'Drop pairs with temporal baseline > '+str(inps.max_temp_baseline)+' days'
        for i in range(len(date12_orig)):
            date1, date2 = date12_orig.split('-')
            idx1 = date6_list.index(date1)
            idx2 = date6_list.index(date2)
            t_diff = tbase_list[idx2] - tbase_list[idx1]
            if t_diff > inps.max_temp_baseline:
                date12_to_rmv.append(date12_orig[i])
    
    # 2.4 Update date12_to_rmv from drop_ifg_index
    if not inps.drop_ifg_index:
        for index in inps.drop_ifg_index:
            date12_to_rmv.append(date12_orig[index-1])
    
    # 2.5 Update date12_to_rmv from drop_date
    if not inps.drop_date:
        inps.drop_date = ptime.yymmdd(inps.drop_date)
        print 'Drop pairs including the following dates: \n'+str(inps.drop_date)
        for i in range(len(date12_orig)):
            date1, date2 = date12_orig.split('-')
            if (date1 in inps.drop_date) or (date2 in inps.drop_date):
                date12_to_rmv.append(date12_orig[i])
    
    # 3. Manually drop pairs
    if inps.disp_network:
        pairs_idx = pnet.read_igram_pairs(inps.file[0])
        bperp_list = ut.Baseline_timeseries(inps.file[0])
        tbase_array = np.array(tbase_list)
        bperp_array = np.array(bperp_list)
        
        # Display the network
        fig1 = plt.figure(1)
        ax1 = fig1.add_subplot(111)
        ax1 = pnet.plot_network(ax1, pairs_idx, date8_list, bperp_list)

        date12_idx = []
        def onclick(event):
            if event.button==1:
                print 'click'
                xClick = event.xdata
                yClick = event.ydata
                idx = nearest_neighbor(xClick, yClick, tbase_array, bperp_array)
                print date8_list[idx]
                date12_idx.append(idx)
                if len(date12_idx) == 2:
                    date12_idx = sorted(date12_idx)
                    date12 = date6_list[date12_idx[0]]+'-'+date6_list[date12_idx[1]]
                    print date12
                    date12_to_rmv.append(date12)
                    ax1.plot(tbase_array[date12_idx], bperp_array[date12_idx], 'r', lw=4)
                    date12_idx = []
            plt.show()
        cid = fig1.canvas.mpl_connect('button_press_event', onclick)
        plt.show()

    # 4. drop duplicate date12 and sort in order
    date12_to_rmv = list(set(date12_to_rmv))
    date12_to_rmv = sorted(date12_to_rmv)

    ############################################################
    for File in inps.file:
        modify_file_date12_list(File, date12_to_rmv)
    print 'Done.'


########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])



