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
import re

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import pysar._datetime as ptime
import pysar._network as pnet
import pysar._pysar_utilities as ut
import pysar._readfile as readfile
from pysar._readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file


###########################  Sub Function  #############################
def nearest_neighbor(x,y, x_array, y_array):
    """ find nearest neighbour
    Input:
        x/y       : float
        x/y_array : numpy.array, temporal/perpendicular spatial baseline
    Output:
        idx : int, index of min distance - nearest neighbour
    """
    dist = np.sqrt((x_array -x)**2 + (y_array -y)**2)
    idx = np.argmin(dist)
    #idx = dist==np.min(dist)
    return idx


def manual_select_pairs_to_remove(File):
    '''Manually select interferograms to remove'''
    print '----------------------------------------------------------------------------'
    print 'Manually select interferograms to remove'
    print 'Click two dates - points - in the figure to select one pair of interferogram'
    print 'repeat until you select all pairs you would like to remove'
    print 'then close the figure to continue the program ...'
    print '----------------------------------------------------------------------------'
    # Display the network
    fig = plt.figure()
    ax = fig.add_subplot(111)

    pairs_idx = pnet.read_igram_pairs(File)
    bperp_list = ut.perp_baseline_ifgram2timeseries(File)[0].tolist()
    date8_list = ptime.ifgram_date_list(File)
    ax = pnet.plot_network(ax, pairs_idx, date8_list, bperp_list)
    print 'display the network of interferogram of file: '+File

    date12_orig = pnet.get_date12_list(File)
    date6_list = ptime.yymmdd(date8_list)
    dates_array = np.array(ptime.date_list2vector(date8_list)[0])
    dateNum_array = mdates.date2num(dates_array)
    bperp_array = np.array(bperp_list)

    date_click = []
    date12_click = []
    def onclick(event):
        xClick = event.xdata
        yClick = event.ydata
        idx = nearest_neighbor(xClick, yClick, dateNum_array, bperp_array)
        date6 = date6_list[idx]
        print 'click at '+date6
        date_click.append(date6)
        if len(date_click)%2 == 0 and date_click[-2] != date_click[-1]:
            [m_date, s_date] = sorted(date_click[-2:])
            m_idx = date6_list.index(m_date)
            s_idx = date6_list.index(s_date)
            date12 = m_date+'-'+s_date
            if date12 in date12_orig:
                print 'select date12: '+date12
                date12_click.append(date12)
                ax.plot([dateNum_array[m_idx],dateNum_array[s_idx]], [bperp_array[m_idx],bperp_array[s_idx]], 'r', lw=4)
            else:
                 print date12+' is not existed in input file'
        plt.draw()
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    return date12_click


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
    
    # Coherence-Based 
    if 'pysar.network.coherenceBase' in keyList:
        if not inps.coherence_file and template_dict['pysar.network.coherenceBase'].lower() in ['yes','y','auto']:
            # Search coherence file from input files
            k_list = [readfile.read_attribute(f)['FILE_TYPE'] for f in inps.file]
            try:  cohFileIdx = k_list.index('coherence')
            except:  sys.exit("ERROR: No coherence found in input files, cannot use coherence-based approach without it.")
            inps.coherence_file = inps.file[cohFileIdx]
            
            # Search mask file
            if not inps.mask_file and os.path.isfile('Mask.h5'):
                inps.mask_file = 'Mask.h5'
    
    return inps


def modify_file_date12_list(File, date12_to_rmv, mark_attribute=False, outFile=None):
    '''Update multiple group hdf5 file using date12 to remove
    Inputs:
        File          - multi_group HDF5 file, i.e. unwrapIfgram.h5, coherence.h5
        date12_to_rmv - list of string indicating interferograms in YYMMDD-YYMMDD format
        mark_attribute- bool, if True, change 'drop_ifgram' attribute only; otherwise, write
                        resutl to a new file
        outFile       - string, output file name
    Output:
        outFile       - string, output file name, if mark_attribute=True, outFile = File
    '''
    k = readfile.read_attribute(File)['FILE_TYPE']
    print '----------------------------------------------------------------------------'
    print 'file: '+File

    if mark_attribute:
        print "set drop_ifgram to 'yes' for all interferograms to remove, and 'no' for all the others."
        h5 = h5py.File(File,'r+')
        ifgram_list = sorted(h5[k].keys())
        for ifgram in ifgram_list:
            if h5[k][ifgram].attrs['DATE12'] in date12_to_rmv:
                h5[k][ifgram].attrs['drop_ifgram'] = 'yes'
            else:
                h5[k][ifgram].attrs['drop_ifgram'] = 'no'
        h5.close()
        outFile = File

    else:
        date12_orig = pnet.get_date12_list(File)
        date12_to_write = sorted(list(set(date12_orig) - set(date12_to_rmv)))
        print 'number of interferograms in file      : '+str(len(date12_orig))
        print 'number of interferograms to keep/write: '+str(len(date12_to_write))
        print 'list   of interferograms to keep/write: '
        print date12_to_write
        date12Num = len(date12_to_write)
    
        if not outFile:
            outFile = 'Modified_'+os.path.basename(File)
        print 'writing >>> '+outFile
        h5out = h5py.File(outFile, 'w')
        gg = h5out.create_group(k)

        h5 = h5py.File(File, 'r')
        igramList = sorted(h5[k].keys())
        for i in range(date12Num):
            date12 = date12_to_write[i]
            idx = date12_orig.index(date12)
            igram = igramList[idx]
            ut.print_progress(i+1, date12Num, prefix='', suffix=igram)
    
            data = h5[k][igram].get(igram)[:]
            group = gg.create_group(igram)
            dset = group.create_dataset(igram, data=data, compression='gzip')
            for key, value in h5[k][igram].attrs.iteritems():
                group.attrs[key] = value
        h5.close()
        h5out.close()
        print 'finished writing >>> '+outFile
    
    return outFile


###############################  Usage  ################################
EXAMPLE='''example:
  modify_network.py unwrapIfgram.h5 coherence.h5 --template KyushuT422F650AlosA.template
  modify_network.py unwrapIfgram.h5 coherence.h5 -t 365 -b 200
  modify_network.py unwrapIfgram.h5 coherence.h5 --coherence-base coherence.h5 --mask Mask.h5 --min-coherence 0.7
  modify_network.py unwrapIfgram.h5 -r Modified_coherence.h5
  modify_network.py unwrapIfgram.h5 --drop-date 20080520 20090816
  modify_network.py unwrapIfgram.h5 --drop-ifg-index 3:9 11 23
  modify_network.py unwrapIfgram.h5 --manual
'''

TEMPLATE='''
pysar.network.dropIfgramIndex = 7:9 15 25 26      #start from 1
pysar.network.dropDate        = 20080520 20090816
pysar.network.maxTempBaseline = 720
pysar.network.maxPerpBaseline = 2000
pysar.network.reference       = Modified_unwrapIfgram.h5
pysar.network.reference       = Paris.list
pysar.network.coherenceBase   = yes    #search and use input coherence file, set to no or comment the line to disable
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Modify the network of interferograms',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)
    parser.add_argument('file', nargs='+',\
                        help='Files to modify/drop network.\n'\
                             'i.e. unwrapIfgram.h5, wrapIfgram.h5, coherence.h5, ...')
    parser.add_argument('--mark-attribute', dest='mark_attribute', action='store_true',\
                        help='mark dropped interferograms in attribute only, do not write new file')
    
    parser.add_argument('-t', dest='max_temp_baseline', type=float, help='temporal baseline threshold/maximum in days')
    parser.add_argument('-b', dest='max_perp_baseline', type=float, help='perpendicular baseline threshold/maximum in meters')

    parser.add_argument('-r','--reference', dest='reference_file',\
                        help='Reference hdf5 / list file with network information.\n'\
                             'i.e. Modified_unwrapIfgram.h5, Pairs.list')
    parser.add_argument('--template', dest='template_file', help='Template file with input options:\n'+TEMPLATE+'\n')
    
    parser.add_argument('--drop-ifg-index', dest='drop_ifg_index', nargs='*',\
                        help='index of interferograms to remove/drop.\n1 as the first')
    parser.add_argument('--drop-date', dest='drop_date', nargs='*',\
                        help='date(s) to remove/drop, all interferograms included date(s) will be removed')
    parser.add_argument('--plot', action='store_true', help='plot and save the result to image files.')
    
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
    manualGroup.add_argument('--manual', dest='disp_network', action='store_true',\
                             help='display network to manually choose line/interferogram to remove')

    inps = parser.parse_args()
    return inps


#########################  Main Function  ##############################
def main(argv):
    ##### Read Inputs
    inps = cmdLineParse()
    inps.file = ut.get_file_list(inps.file)
    date12_orig = pnet.get_date12_list(inps.file[0])
    print 'input file(s) to be modified: '+str(inps.file)
    print 'number of interferograms: '+str(len(date12_orig))
    #print '\n****************** Network Modification ********************'

    if (not inps.reference_file and not inps.template_file and\
        not inps.max_temp_baseline and not inps.max_perp_baseline and\
        not inps.drop_ifg_index and not inps.drop_date and \
        not inps.coherence_file):
        # Display network for manually modification when there is no other modification input.
        print 'No input found to remove interferogram, continue by display the network to select it manually ...'
        inps.disp_network = True

    # Update inps if template is input
    if inps.template_file:
        inps = update_inps_with_template(inps, inps.template_file)
    
    # Convert index : input to continous index list
    if inps.drop_ifg_index:
        ifg_index = list(inps.drop_ifg_index)
        inps.drop_ifg_index = []
        for index in ifg_index:
            index_temp = [int(i) for i in index.split(':')]
            index_temp.sort()
            if len(index_temp)==2:
                for j in range(index_temp[0], index_temp[1]+1):
                    inps.drop_ifg_index.append(j)
            elif len(index_temp)==1:
                inps.drop_ifg_index.append(int(index))
            else:
                print 'Unrecoganized input: '+index
        inps.drop_ifg_index = sorted(inps.drop_ifg_index)
        if max(inps.drop_ifg_index) > len(date12_orig):
            raise Exception('Input index out of range!\n'+\
                            'input index:'+str(inps.drop_ifg_index)+'\n'+\
                            'index range of file: '+str(len(date12_orig)))

    ##### Get date12_to_rmv
    date12_to_rmv = []

    # 1. Update date12_to_rmv from reference file
    if inps.reference_file:
        date12_to_keep = pnet.get_date12_list(inps.reference_file)
        print '----------------------------------------------------------------------------'
        print 'use reference pairs info from file: '+inps.reference_file
        print 'number of interferograms in reference: '+str(len(date12_to_keep))
        print 'date12 not in reference file:'
        for date12 in date12_orig:
            if date12 not in date12_to_keep:
                date12_to_rmv.append(date12)
                print date12

    # 2.1 Update date12_to_rmv from coherence file
    if inps.coherence_file:
        print '----------------------------------------------------------------------------'
        print 'use coherence-based network modification from coherence file: '+inps.coherence_file
        # Calculate spatial average coherence
        if not inps.mask_file:
            mask = readfile.read(inps.mask_file)[0]
            print 'mask coherence with file: '+inps.mask_file
        else:
            mask = None
        cohTextFile = os.path.splitext(inps.coherence_file)[0]+'_spatialAverage.list'
        if os.path.isfile(cohTextFile):
            print 'average coherence in space has been calculated before and store in file: '+cohTextFile
            print 'read it directly, or delete it and re-run the script to re-calculate the list'
            cohTxt = np.loadtxt(cohTextFile, dtype=str)
            mean_coherence_list = [float(i) for i in cohTxt[:,1]]
            coh_date12_list = [i for i in cohTxt[:,0]]
        else:
            print 'calculating average coherence of each interferogram ...'
            mean_coherence_list = ut.spatial_average(inps.coherence_file, mask, saveList=True)
            coh_date12_list = pnet.get_date12_list(inps.coherence_file)
        print 'date12 with average coherence < '+str(inps.min_coherence)+': '
        for i in range(len(coh_date12_list)):
            if mean_coherence_list[i] < inps.min_coherence:
                date12 = coh_date12_list[i]
                date12_to_rmv.append(date12)
                print date12

    # 2.2 Update date12_to_rmv from perp baseline threshold
    if inps.max_perp_baseline:
        print '----------------------------------------------------------------------------'
        print 'Drop pairs with perpendicular spatial baseline > '+str(inps.max_perp_baseline)+' meters'
        ifg_bperp_list = pnet.igram_perp_baseline_list(inps.file[0])
        for i in range(len(ifg_bperp_list)):
            if ifg_bperp_list[i] > inps.max_perp_baseline:
                date12 = date12_orig[i]
                date12_to_rmv.append(date12)
                print date12

    # 2.3 Update date12_to_rmv from temp baseline threshold
    if inps.max_temp_baseline:
        print '----------------------------------------------------------------------------'
        print 'Drop pairs with temporal baseline > '+str(inps.max_temp_baseline)+' days'
        date8_list = ptime.ifgram_date_list(inps.file[0])
        date6_list = ptime.yymmdd(date8_list)
        tbase_list = ptime.date_list2tbase(date8_list)[0]
        for i in range(len(date12_orig)):
            date1, date2 = date12_orig[i].split('-')
            idx1 = date6_list.index(date1)
            idx2 = date6_list.index(date2)
            t_diff = tbase_list[idx2] - tbase_list[idx1]
            if t_diff > inps.max_temp_baseline:
                date12 = date12_orig[i]
                date12_to_rmv.append(date12)
                print date12

    # 2.4 Update date12_to_rmv from drop_ifg_index
    if inps.drop_ifg_index:
        print '----------------------------------------------------------------------------'
        print 'drop date12/pair with the following index number:'
        for index in inps.drop_ifg_index:
            date12 = date12_orig[index-1]
            date12_to_rmv.append(date12)
            print str(index)+'    '+date12

    # 2.5 Update date12_to_rmv from drop_date
    if inps.drop_date:
        inps.drop_date = ptime.yymmdd(inps.drop_date)
        print '----------------------------------------------------------------------------'
        print 'Drop pairs including the following dates: \n'+str(inps.drop_date)
        for i in range(len(date12_orig)):
            date1, date2 = date12_orig[i].split('-')
            if (date1 in inps.drop_date) or (date2 in inps.drop_date):
                date12 = date12_orig[i]
                date12_to_rmv.append(date12)
                print date12

    # 3. Manually drop pairs
    if inps.disp_network:
        date12_click = manual_select_pairs_to_remove(inps.file[0])
        for date12 in list(date12_click):
            if date12 not in date12_orig:
                date12_click.remove(date12)
        print 'date12 selected to remove:'
        print date12_click
        date12_to_rmv += date12_click

    # 4. drop duplicate date12 and sort in order
    date12_to_rmv = list(set(date12_to_rmv))
    date12_to_rmv = sorted(date12_to_rmv)
    print '----------------------------------------------------------------------------'
    print 'number of interferograms to remove: '+str(len(date12_to_rmv))
    print 'list   of interferograms to remove:'
    print date12_to_rmv

    # Check existing mark for --mark-attribute option
    if inps.mark_attribute:
        # Get list of date12 of interferograms already been marked.
        k = readfile.read_attribute(inps.file[0])['FILE_TYPE']
        atr = readfile.read_attribute(inps.file[0])
        h5 = h5py.File(inps.file[0], 'r')
        ifgram_list_all = sorted(h5[k].keys())
        ifgram_list_keep = ut.check_drop_ifgram(h5, atr, ifgram_list_all)
        ifgram_list_dropped = sorted(list(set(ifgram_list_all) - set(ifgram_list_keep)))
        date12_list_dropped = [str(re.findall('\d{6}-\d{6}', i)[0]) for i in ifgram_list_dropped]
        h5.close()
        
        if date12_to_rmv == date12_list_dropped:
            date12_to_rmv = []
            print 'calculated date12 to drop is the same as exsiting marked input file, set it empty.'

    if date12_to_rmv:
        ##### Update Input Files with date12_to_rmv
        Modified_CoherenceFile = 'Modified_coherence.h5'
        for File in inps.file:
            Modified_File = modify_file_date12_list(File, date12_to_rmv, inps.mark_attribute)

            k = readfile.read_attribute(File)['FILE_TYPE']
            # Update Mask File
            if k == 'interferograms':
                print 'update mask file for input '+k+' file based on '+Modified_File
                inps.mask_file = 'maskModified.h5'
                print 'writing >>> '+inps.mask_file
                ut.nonzero_mask(Modified_File, inps.mask_file)
            elif k == 'coherence':
                print 'update average spatial coherence for input '+k+' file based on: '+Modified_File
                outFile = 'averageSpatialCoherenceModified.h5'
                print 'writing >>> '+outFile
                ut.temporal_average(Modified_File, outFile)
                Modified_CoherenceFile = Modified_File
    
        # Plot result
        if inps.plot:
            print '\nplot modified network and save to file.'
            plotCmd = 'plot_network.py '+Modified_File+' --coherence '+Modified_CoherenceFile+' --nodisplay'
            if inps.mask_file:
                plotCmd += ' --mask '+inps.mask_file
            print plotCmd
            os.system(plotCmd)
        
        print 'Done.'
        return
    else:
        print 'No interferogram dropped, skip update.'
        return


########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])



