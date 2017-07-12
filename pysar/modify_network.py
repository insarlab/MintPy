#! /usr/bin/env python2
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

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import pysar._datetime as ptime
import pysar._network as pnet
import pysar._pysar_utilities as ut
import pysar._readfile as readfile
import pysar.subset as subset
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


def reset_pairs(File):
    '''Reset/restore all pairs within the input file by set all drop_ifgram=no'''
    print "set drop_ifgram to 'no' for all interferograms for file: "+File
    k = readfile.read_attribute(File)['FILE_TYPE']
    h5 = h5py.File(File,'r+')
    ifgram_list = sorted(h5[k].keys())
    for ifgram in ifgram_list:
        h5[k][ifgram].attrs['drop_ifgram'] = 'no'
    h5.close()
    return File


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

    date12_orig = pnet.get_date12_list(File)
    bperp_list = ut.perp_baseline_ifgram2timeseries(File)[0].tolist()
    date8_list = ptime.ifgram_date_list(File)
    ax = pnet.plot_network(ax, date12_orig, date8_list, bperp_list)
    print 'display the network of interferogram of file: '+File

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
        date12_list = ptime.list_ifgram2date12(igramList)
        prog_bar = ptime.progress_bar(maxValue=date12Num, prefix='writing: ')
        for i in range(date12Num):
            date12 = date12_to_write[i]
            idx = date12_orig.index(date12)
            igram = igramList[idx]
    
            data = h5[k][igram].get(igram)[:]
            group = gg.create_group(igram)
            dset = group.create_dataset(igram, data=data, compression='gzip')
            for key, value in h5[k][igram].attrs.iteritems():
                group.attrs[key] = value
            prog_bar.update(i+1, suffix=date12_list[i])
        prog_bar.close()
        h5.close()
        h5out.close()
        print 'finished writing >>> '+outFile
    
    return outFile


def read_template2inps(template_file, inps=None):
    '''Read input template options into Namespace inps'''
    if not inps:
        inps = cmdLineParse()

    template = readfile.read_template(inps.template_file)
    key_list = template.keys()

    # Coherence-based network modification
    prefix = 'pysar.network.'

    key = prefix+'coherenceBased'
    if key in key_list and template[key] in ['auto','yes']:
        inps.coherence_based = True

    key = prefix+'coherenceFile'
    if key in key_list:
        if template[key] == 'auto':
            inps.coherence_file = 'coherence.h5'
        else:
            inps.coherence_file = template[key]

    # find coherence file from input files if inps.coherence_file does not exists.
    if inps.coherence_based and not os.path.isfile(inps.coherence_file):
        k_list = [readfile.read_attribute(f)['FILE_TYPE'] for f in inps.file]
        try:
            coh_file_idx = k_list.index('coherence')
        except ValueError:
            print 'No coherence file found! Can not use coherence-based method without it.'
        inps.coherence_file = inps.file[coh_file_idx]

    key = prefix+'minCoherence'
    if key in key_list:
        if template[key] == 'auto':
            inps.min_coherence = 0.7
        else:
            inps.min_coherence = float(template[key])

    key = prefix+'maskFile'
    if key in key_list:
        value = template[key]
        if value == 'auto':
            inps.mask_file = 'mask.h5'
        elif value == 'no':
            inps.mask_file = None
        else:
            inps.mask_file = value

    key = prefix+'maskAoi.yx'
    if key in key_list:
        value = template[key]
        if value in ['auto','no']:
            inps.aoi_pix_box = None
        else:
            tmp = [i.strip() for i in value.split(',')]
            sub_y = sorted([int(i.strip()) for i in tmp[0].split(':')])
            sub_x = sorted([int(i.strip()) for i in tmp[1].split(':')])
            inps.aoi_pix_box = (sub_x[0], sub_y[0], sub_x[1], sub_y[1])

    key = prefix+'maskAoi.lalo'
    if key in key_list:
        value = template[key]
        if value in ['auto','no']:
            inps.aoi_geo_box = None
        else:
            tmp = [i.strip() for i in value.split(',')]
            sub_lat = sorted([float(i.strip()) for i in tmp[0].split(':')])
            sub_lon = sorted([float(i.strip()) for i in tmp[1].split(':')])
            inps.aoi_geo_box = (sub_lon[0], sub_lat[1], sub_lon[1], sub_lat[0])
            # Check trans file
            try:
                inps.trans_file = ut.get_file_list(inps.trans_file)[0]
            except:
                inps.trans_file = None
                print 'Warning: no mapping transformation file found! Can not use '+key+' option without it.'
                print 'skip this option.'
                inps.aoi_pix_box = None

    ## Network Modification based on thresholds
    key = prefix+'tempBaseMax'
    if key in key_list:
        value = template[key]
        if value not in ['auto','no']:
            inps.max_temp_baseline = float(value)

    key = prefix+'perpBaseMax'
    if key in key_list:
        value = template[key]
        if value not in ['auto','no']:
            inps.max_perp_baseline = float(value)

    key = prefix+'referenceFile'
    if key in key_list:
        value = template[key]
        if value in ['auto','no']:
            inps.reference_file = None
        else:
            inps.reference_file = value

    key = prefix+'excludeDate'
    if key in key_list:
        value = template[key]
        if value not in ['auto','no']:
            inps.exclude_date = [i for i in value.replace(',',' ').split()]

    key = prefix+'excludeIfgIndex'
    if key in key_list:
        value = template[key]
        if value not in ['auto','no']:
            inps.exclude_ifg_index = [i for i in value.replace(',',' ').split()]

    return inps


###############################  Usage  ################################
EXAMPLE='''example:
  modify_network.py unwrapIfgram.h5 coherence.h5 --template pysarApp_template.txt --trans geomap_4rlks.trans
  modify_network.py unwrapIfgram.h5 coherence.h5 -t 365 -b 200
  modify_network.py unwrapIfgram.h5 coherence.h5 --coherence-base coherence.h5 --mask Mask.h5 --min-coherence 0.7
  modify_network.py unwrapIfgram.h5 -r Modified_coherence.h5
  modify_network.py unwrapIfgram.h5 --start-date 20080520  --end-date 20110101
  modify_network.py unwrapIfgram.h5 --exclude-date 20080520 20090816
  modify_network.py unwrapIfgram.h5 --exclude-ifg-index 3:9 11 23
  modify_network.py unwrapIfgram.h5 --manual
'''

TEMPLATE='''
pysar.network.coherenceBased     = yes      #[yes / no], auto for yes
pysar.network.coherenceFile      = auto     #[filename], auto for coherence.h5
pysar.network.minCoherence       = auto     #[0.0-1.0], auto for 0.7
pysar.network.maskFile           = auto     #[file name, no], auto for mask.h5, no for all pixels
pysar.network.maskAoi.yx         = no       #[y0:y1,x0:x1 / no], area of interest for coherence calculation, auto for no
pysar.network.maskAoi.lalo       = no       #[lat0:lat1,lon0:lon1 / no], similar to maskAoi.yx but in lat/lon, auto for no

pysar.network.tempBaseMax        = 36500    #[1-inf], day, maximum temporal baseline, auto for 3.65e4
pysar.network.perpBaseMax        = 10e3     #[1-inf], meter, maximum perpendicular spatial baseline, auto for 10e3
pysar.network.referenceFile      = no       #[date12_list.txt / Modified_unwrapIfgram.h5 / no], auto for no
pysar.network.excludeDate        = no       #[20080520,20090817 / no], auto for no
pysar.network.excludeIfgIndex = no       #[1,3,25 / no], list of interferogram number starting from 1, auto for no
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Modify the network of interferograms',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)
    parser.add_argument('file', nargs='+',\
                        help='Files to modify/drop network.\n'\
                             'i.e. unwrapIfgram.h5, wrapIfgram.h5, coherence.h5, ...')
    parser.add_argument('--reset', action='store_true',\
                        help='restore all interferograms existed in the file, by marking all drop_ifgram=no')
    parser.add_argument('--write-file', dest='mark_attribute', action='store_false',\
                        help='mark dropped interferograms in attribute only, do not write new file')
    parser.add_argument('--plot', action='store_true', help='plot and save the result to image files.')

    parser.add_argument('-t', dest='max_temp_baseline', type=float, help='temporal baseline threshold/maximum in days')
    parser.add_argument('-b', dest='max_perp_baseline', type=float, help='perpendicular baseline threshold/maximum in meters')

    parser.add_argument('-r','--reference', dest='reference_file',\
                        help='Reference hdf5 / list file with network information.\n'\
                             'i.e. Modified_unwrapIfgram.h5, Pairs.list')
    parser.add_argument('--template', dest='template_file', help='Template file with input options:\n'+TEMPLATE+'\n')

    parser.add_argument('--exclude-ifg-index', dest='exclude_ifg_index', nargs='*',\
                        help='index of interferograms to remove/drop.\n1 as the first')
    parser.add_argument('--exclude-date', dest='exclude_date', nargs='*',\
                        help='date(s) to remove/drop, all interferograms included date(s) will be removed')
    parser.add_argument('--start-date','--min-date', dest='start_date',\
                        help='remove/drop interferograms with date earlier than start-date in YYMMDD or YYYYMMDD format')
    parser.add_argument('--end-date','--max-date', dest='end_date',\
                        help='remove/drop interferograms with date later than end-date in YYMMDD or YYYYMMDD format')

    # Coherence-based network
    cohBased = parser.add_argument_group('Coherence-based Network',\
                                               'Drop/modify network based on spatial coherence')
    cohBased.add_argument('--coherence-based', dest='coherence_based', action='store_true',\
                          help='Enable coherence-based network modification')
    cohBased.add_argument('--coherence', dest='coherence_file', default='coherence.h5',\
                          help='Coherence file used to calculate average value for each interferograms\n'+\
                               'Input coherence file should have the same network as input file(s)\n'+\
                               'default: coherence.h5')
    cohBased.add_argument('--mask', dest='mask_file',\
                          help='Mask file used to calculate the spatial coherence\n'\
                               'Will use the whole area if not assigned')
    cohBased.add_argument('--min-coherence', dest='min_coherence', type=float, default=0.7,\
                          help='Minimum coherence value')
    cohBased.add_argument('--trans', dest='trans_file', default='geomap*.trans',\
                          help='mapping transformation file for geo/radar coordinate conversion.\n'+\
                               'Needed for mask AOI in lalo')

    # Manually select network
    manual = parser.add_argument_group('Manual Network', 'Manually select/drop/modify network')
    manual.add_argument('--manual', dest='disp_network', action='store_true',\
                        help='display network to manually choose line/interferogram to remove')

    inps = parser.parse_args()
    inps.aoi_geo_box = None
    inps.aoi_pix_box = None
    return inps


#########################  Main Function  ##############################
def main(argv):
    ##### Read Inputs
    inps = cmdLineParse()
    inps.file = ut.get_file_list(inps.file)
    date12_orig = pnet.get_date12_list(inps.file[0])
    print 'input file(s) to be modified: '+str(inps.file)
    print 'number of interferograms: '+str(len(date12_orig))
    atr = readfile.read_attribute(inps.file[0])

    #print '\n****************** Network Modification ********************'

    if inps.reset:
        print '----------------------------------------------------------------------------'
        for file in inps.file:
            reset_pairs(file)
        return

    # Update inps if template is input
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    if all(not i for i in [inps.reference_file, inps.max_temp_baseline, inps.max_perp_baseline,\
                           inps.exclude_ifg_index, inps.exclude_date, inps.coherence_based, inps.start_date, inps.end_date]):
        # Display network for manually modification when there is no other modification input.
        print 'No input option found to remove interferogram'
        print 'To manually modify network, please use --manual option '
        return
    
    # Convert index : input to continous index list
    if inps.exclude_ifg_index:
        ifg_index = list(inps.exclude_ifg_index)
        inps.exclude_ifg_index = []
        for index in ifg_index:
            index_temp = [int(i) for i in index.split(':')]
            index_temp.sort()
            if len(index_temp)==2:
                for j in range(index_temp[0], index_temp[1]+1):
                    inps.exclude_ifg_index.append(j)
            elif len(index_temp)==1:
                inps.exclude_ifg_index.append(int(index))
            else:
                print 'Unrecoganized input: '+index
        inps.exclude_ifg_index = sorted(inps.exclude_ifg_index)
        if max(inps.exclude_ifg_index) > len(date12_orig):
            raise Exception('Input index out of range!\n'+\
                            'input index:'+str(inps.exclude_ifg_index)+'\n'+\
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
    if inps.coherence_based and os.path.isfile(inps.coherence_file):
        print '----------------------------------------------------------------------------'
        print 'use coherence-based network modification from coherence file: '+inps.coherence_file
        # check mask AOI in lalo
        if inps.aoi_geo_box and inps.trans_file:
            print 'input AOI in (lon0, lat1, lon1, lat0): '+str(inps.aoi_geo_box)
            inps.aoi_pix_box = subset.bbox_geo2radar(inps.aoi_geo_box, atr, inps.trans_file) 
        if inps.aoi_pix_box:
            print 'input AOI in (x0,y0,x1,y1): '+str(inps.aoi_pix_box)

        # Calculate spatial average coherence
        coh_list, coh_date12_list = ut.get_spatial_average(inps.coherence_file, inps.mask_file,\
                                                           inps.aoi_pix_box, saveList=True)

        # MST network
        print 'Get minimum spanning tree (MST) of interferograms with inverse of coherence.'
        mst_date12_list = pnet.threshold_coherence_based_mst(coh_date12_list, coh_list)

        print 'date12 with average coherence < '+str(inps.min_coherence)+' and not in MST: '
        for i in range(len(coh_date12_list)):
            date12 = coh_date12_list[i]
            if coh_list[i] < inps.min_coherence and date12 not in mst_date12_list:
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

    # 2.4 Update date12_to_rmv from exclude_ifg_index
    if inps.exclude_ifg_index:
        print '----------------------------------------------------------------------------'
        print 'drop date12/pair with the following index number:'
        for index in inps.exclude_ifg_index:
            date12 = date12_orig[index-1]
            date12_to_rmv.append(date12)
            print str(index)+'    '+date12

    # 2.5 Update date12_to_rmv from exclude_date
    if inps.exclude_date:
        inps.exclude_date = ptime.yymmdd(inps.exclude_date)
        print '----------------------------------------------------------------------------'
        print 'Drop pairs including the following dates: \n'+str(inps.exclude_date)
        for i in range(len(date12_orig)):
            date1, date2 = date12_orig[i].split('-')
            if (date1 in inps.exclude_date) or (date2 in inps.exclude_date):
                date12 = date12_orig[i]
                date12_to_rmv.append(date12)
                print date12

    # 2.6 Update date12_to_rmv from start_date
    if inps.start_date:
        inps.start_date = ptime.yymmdd(inps.start_date)
        print '----------------------------------------------------------------------------'
        print 'Drop pairs with date earlier than start-date: '+inps.start_date
        min_date = int(ptime.yyyymmdd(inps.start_date))
        for i in range(len(date12_orig)):
            date12 = date12_orig[i]
            if any(int(j) < min_date for j in ptime.yyyymmdd(date12.split('-'))):
                date12_to_rmv.append(date12)
                print date12

    # 2.7 Update date12_to_rmv from end_date
    if inps.end_date:
        inps.end_date = ptime.yymmdd(inps.end_date)
        print '----------------------------------------------------------------------------'
        print 'Drop pairs with date earlier than end-date: '+inps.end_date
        max_date = int(ptime.yyyymmdd(inps.end_date))
        for i in range(len(date12_orig)):
            date12 = date12_orig[i]
            if any(int(j) > max_date for j in ptime.yyyymmdd(date12.split('-'))):
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
        h5 = h5py.File(inps.file[0], 'r')
        ifgram_list_all = sorted(h5[k].keys())
        ifgram_list_keep = ut.check_drop_ifgram(h5, atr, ifgram_list_all, print_message=False)
        ifgram_list_dropped = sorted(list(set(ifgram_list_all) - set(ifgram_list_keep)))
        date12_list_dropped = ptime.list_ifgram2date12(ifgram_list_dropped)
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
                inps.mask_file = 'mask.h5'
                print 'writing >>> '+inps.mask_file
                ut.nonzero_mask(Modified_File, inps.mask_file)
            elif k == 'coherence':
                print 'update average spatial coherence for input '+k+' file based on: '+Modified_File
                outFile = 'averageSpatialCoherence.h5'
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



