############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import h5py
import numpy as np
from matplotlib import pyplot as plt, dates as mdates

from mintpy.objects import ifgramStack
from mintpy.utils import (
    ptime,
    readfile,
    utils as ut,
    network as pnet,
    plot as pp,
)


###############################  Usage  ################################
def read_input_index_list(idxList, stackFile=None):
    """Read ['2','3:5','10'] into ['2','3','4','5','10']"""
    idxListOut = []
    for idx in idxList:
        c = sorted([int(i) for i in idx.split(':')])
        if len(c) == 2:
            idxListOut += list(range(c[0], c[1]+1))
        elif len(c) == 1:
            idxListOut.append(c[0])
        else:
            print('Unrecoganized input: '+idx)
    idxListOut = sorted(set(idxListOut))

    if stackFile:
        obj = ifgramStack(stackFile)
        obj.open(print_msg=False)
        idxListOut = [i for i in idxListOut if i < obj.numIfgram]
        obj.close(print_msg=False)
    return idxListOut


def read_template2inps(template_file, inps):
    """Read input template options into Namespace inps"""
    inpsDict = vars(inps)
    print('read options from template file: '+os.path.basename(template_file))
    template = readfile.read_template(inps.template_file, skip_chars=['[', ']'])
    template = ut.check_template_auto_value(template)

    # Update inps if key existed in template file
    prefix = 'mintpy.network.'
    keyList = [i for i in list(inpsDict.keys()) if prefix+i in template.keys()]
    for key in keyList:
        value = template[prefix+key]
        if key in ['coherenceBased', 'areaRatioBased', 'keepMinSpanTree']:
            inpsDict[key] = value
        elif value:
            if key in ['minCoherence', 'minAreaRatio', 'tempBaseMax', 'perpBaseMax']:
                inpsDict[key] = float(value)
            elif key in ['connNumMax']:
                inpsDict[key] = int(value)
            elif key in ['maskFile', 'referenceFile']:
                inpsDict[key] = value
            elif key == 'aoiYX':
                tmp = [i.strip() for i in value.split(',')]
                sub_y = sorted([int(i.strip()) for i in tmp[0].split(':')])
                sub_x = sorted([int(i.strip()) for i in tmp[1].split(':')])
                inps.aoi_pix_box = (sub_x[0], sub_y[0], sub_x[1], sub_y[1])
            elif key == 'aoiLALO':
                tmp = [i.strip() for i in value.split(',')]
                sub_lat = sorted([float(i.strip()) for i in tmp[0].split(':')])
                sub_lon = sorted([float(i.strip()) for i in tmp[1].split(':')])
                inps.aoi_geo_box = (sub_lon[0], sub_lat[1], sub_lon[1], sub_lat[0])
                # Check lookup file
                if not inps.lookupFile:
                    print('Warning: no lookup table file found! Can not use '+key+' option without it.')
                    print('skip this option.')
                    inps.aoi_pix_box = None
            elif key in ['startDate', 'endDate']:
                inpsDict[key] = ptime.yyyymmdd(value)
            elif key == 'excludeDate':
                inpsDict[key] = ptime.yyyymmdd(value.split(','))
            elif key == 'excludeIfgIndex':
                inpsDict[key] += value.split(',')
                inpsDict[key] = read_input_index_list(inpsDict[key], stackFile=inps.file)

    # Turn reset on if 1) no input options found to drop ifgram AND 2) there is template input
    if all(not i for i in [inps.referenceFile, inps.tempBaseMax, inps.perpBaseMax, inps.connNumMax,
                           inps.excludeIfgIndex, inps.excludeDate, inps.coherenceBased, inps.areaRatioBased,
                           inps.startDate, inps.endDate, inps.reset, inps.manual]):
        print('No input option found to remove interferogram')
        print('Keep all interferograms by enable --reset option')
        inps.reset = True
    return inps


###########################  Sub Function  #############################
def reset_network(stackFile):
    """Reset/restore all pairs within the input file by set all DROP_IFGRAM=no"""
    print("reset dataset 'dropIfgram' to True for all interferograms for file: "+stackFile)
    obj = ifgramStack(stackFile)
    obj.open(print_msg=False)
    if np.all(obj.dropIfgram):
        print('All dropIfgram are already True, no need to reset.')
    else:
        with h5py.File(stackFile, 'r+') as f:
            f['dropIfgram'][:] = True
        ut.touch('coherenceSpatialAvg.txt')
    return stackFile


def nearest_neighbor(x, y, x_array, y_array):
    """ find nearest neighbour
    Input:
        x/y       : float
        x/y_array : numpy.array, temporal/perpendicular spatial baseline
    Output:
        idx : int, index of min distance - nearest neighbour
    """
    dist = np.sqrt((x_array - x)**2 + (y_array - y)**2)
    idx = np.argmin(dist)
    #idx = dist==np.min(dist)
    return idx


def manual_select_pairs_to_remove(stackFile):
    """Manually select interferograms to remove"""
    print('\n-------------------------------------------------------------')
    print('Manually select interferograms to remove')
    print('1) click two dates/points to select one pair of interferogram')
    print('2) repeat until you select all pairs you would like to remove')
    print('3) close the figure to continue the program ...')
    print('-------------------------------------------------------------\n')
    obj = ifgramStack(stackFile)
    obj.open()
    date12ListAll = obj.date12List
    pbase = obj.get_perp_baseline_timeseries(dropIfgram=False)
    dateList = obj.dateList
    datesNum = mdates.date2num(np.array(ptime.date_list2vector(dateList)[0]))

    date12ListKept = obj.get_date12_list(dropIfgram=True)
    date12ListDropped = sorted(list(set(date12ListAll) - set(date12ListKept)))

    # Display the network
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax = pp.plot_network(ax, date12ListAll, dateList, list(pbase), date12List_drop=date12ListDropped)
    print('display the network of interferogram of file: '+stackFile)
    date_click = []
    date12_click = []

    def onclick(event):
        idx = nearest_neighbor(event.xdata, event.ydata, datesNum, pbase)
        print('click at '+dateList[idx])
        date_click.append(dateList[idx])
        if len(date_click) % 2 == 0 and date_click[-2] != date_click[-1]:
            [mDate, sDate] = sorted(date_click[-2:])
            mIdx = dateList.index(mDate)
            sIdx = dateList.index(sDate)
            date12 = mDate+'_'+sDate
            if date12 in date12ListAll:
                print('select date12: '+date12)
                date12_click.append(date12)
                ax.plot([datesNum[mIdx], datesNum[sIdx]], [pbase[mIdx], pbase[sIdx]], 'r', lw=4)
            else:
                print(date12+' is not existed in input file')
        plt.draw()
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()

    if not ut.yes_or_no('Proceed to drop the ifgrams/date12?'):
        date12_click = None

    return date12_click


def get_aoi_pix_box(meta, lookup_file, pix_box, geo_box):
    """Get pix_box for AOI."""
    coord = ut.coordinate(meta, lookup_file=lookup_file)

    # geo_box -> pix_box
    if geo_box and lookup_file:
        print('input AOI in (lon0, lat1, lon1, lat0): {}'.format(geo_box))
        pix_box = coord.bbox_geo2radar(geo_box)

    # check pix_box
    if pix_box:
        pix_box = coord.check_box_within_data_coverage(pix_box)
        print('input AOI in (x0,y0,x1,y1): {}'.format(pix_box))

    return pix_box


def get_mst_date12(keep_mst, par_list_all, date12_list_all, date12_to_drop, min_par, par_name='average coherence'):
    """Get the date12_list of the MST network for the given parameter."""
    if keep_mst:
        print(f'Get minimum spanning tree (MST) of interferograms with inverse of {par_name}.')
        msg = ('Drop ifgrams with '
               '1) {} < {} AND '
               '2) not in MST network: '.format(par_name, min_par))

        # get the current remaining network (after all the above criteria and before data-driven)
        date12_to_keep = sorted(list(set(date12_list_all) - set(date12_to_drop)))
        par_to_keep = [par for par, date12 in zip(par_list_all, date12_list_all)
                       if date12 in date12_to_keep]

        # get MST from the current remaining network
        mst_date12_list = pnet.threshold_coherence_based_mst(date12_to_keep, par_to_keep)
        mst_date12_list = ptime.yyyymmdd_date12(mst_date12_list)

    else:
        msg = 'Drop ifgrams with {} < {}: '.format(par_name, min_par)
        mst_date12_list = []

    return mst_date12_list, msg


def get_date12_to_drop(inps):
    """Get date12 list to dropped
    Return [] if no ifgram to drop, thus KEEP ALL ifgrams;
           None if nothing to change, exit without doing anything.
    """
    obj = ifgramStack(inps.file)
    obj.open()
    date12ListAll = obj.date12List
    dateList = obj.dateList
    print('number of interferograms: {}'.format(len(date12ListAll)))

    # Get date12_to_drop
    date12_to_drop = []

    # reference file
    if inps.referenceFile:
        date12_to_keep = pnet.get_date12_list(inps.referenceFile, dropIfgram=True)
        print('--------------------------------------------------')
        print('use reference pairs info from file: {}'.format(inps.referenceFile))
        print('number of interferograms in reference: {}'.format(len(date12_to_keep)))
        tempList = sorted(list(set(date12ListAll) - set(date12_to_keep)))
        date12_to_drop += tempList
        print('date12 not in reference file: ({})\n{}'.format(len(tempList), tempList))

    # temp baseline threshold
    if inps.tempBaseMax:
        tempIndex = np.abs(obj.tbaseIfgram) > inps.tempBaseMax
        tempList = list(np.array(date12ListAll)[tempIndex])
        date12_to_drop += tempList
        print('--------------------------------------------------')
        print('Drop ifgrams with temporal baseline > {} days: ({})\n{}'.format(
            inps.tempBaseMax, len(tempList), tempList))

    # perp baseline threshold
    if inps.perpBaseMax:
        tempIndex = np.abs(obj.pbaseIfgram) > inps.perpBaseMax
        tempList = list(np.array(date12ListAll)[tempIndex])
        date12_to_drop += tempList
        print('--------------------------------------------------')
        print('Drop ifgrams with perp baseline > {} meters: ({})\n{}'.format(
            inps.perpBaseMax, len(tempList), tempList))

    # connection number threshold
    if inps.connNumMax:
        seq_date12_list = pnet.select_pairs_sequential(dateList, inps.connNumMax)
        seq_date12_list = ptime.yyyymmdd_date12(seq_date12_list)
        tempList = [i for i in date12ListAll if i not in seq_date12_list]
        date12_to_drop += tempList
        print('--------------------------------------------------')
        msg = 'Drop ifgrams with temporal baseline beyond {} neighbors: ({})'.format(
            inps.connNumMax, len(tempList))
        if len(tempList) <= 200:
            msg += '\n{}'.format(tempList)
        print(msg)

    # excludeIfgIndex
    if inps.excludeIfgIndex:
        tempList = [date12ListAll[i] for i in inps.excludeIfgIndex]
        date12_to_drop += tempList
        print('--------------------------------------------------')
        print('Drop ifgrams with the following index number: {}'.format(len(tempList)))
        for i, date12 in enumerate(tempList):
            print('{} : {}'.format(i, date12))

    # excludeDate
    if inps.excludeDate:
        tempList = [i for i in date12ListAll if any(j in inps.excludeDate for j in i.split('_'))]
        date12_to_drop += tempList
        print('-'*50+'\nDrop ifgrams including the following dates: ({})\n{}'.format(
            len(tempList), inps.excludeDate))
        print('-'*30+'\n{}'.format(tempList))

    # startDate
    if inps.startDate:
        minDate = int(inps.startDate)
        tempList = [i for i in date12ListAll if any(int(j) < minDate for j in i.split('_'))]
        date12_to_drop += tempList
        print('--------------------------------------------------')
        print('Drop ifgrams with date earlier than: {} ({})\n{}'.format(
            inps.startDate, len(tempList), tempList))

    # endDate
    if inps.endDate:
        maxDate = int(inps.endDate)
        tempList = [i for i in date12ListAll if any(int(j) > maxDate for j in i.split('_'))]
        date12_to_drop += tempList
        print('--------------------------------------------------')
        print('Drop ifgrams with date later than: {} ({})\n{}'.format(
            inps.endDate, len(tempList), tempList))

    # coherence file
    if inps.coherenceBased:
        print('--------------------------------------------------')
        print('use coherence-based network modification')

        # get area of interest for coherence calculation
        pix_box = get_aoi_pix_box(obj.metadata, inps.lookupFile, inps.aoi_pix_box, inps.aoi_geo_box)

        # calculate spatial average coherence
        cohList = ut.spatial_average(inps.file,
                                     datasetName='coherence',
                                     maskFile=inps.maskFile,
                                     box=pix_box,
                                     saveList=True)[0]

        # get coherence-based network
        coh_date12_list = list(np.array(date12ListAll)[np.array(cohList) >= inps.minCoherence])

        # get MST network
        mst_date12_list, msg = get_mst_date12(inps.keepMinSpanTree, cohList, date12ListAll, date12_to_drop,
                                              min_par=inps.minCoherence,
                                              par_name='average coherence')

        # drop all dates (below cohh thresh AND not in MST)
        tempList = sorted(list(set(date12ListAll) - set(coh_date12_list + mst_date12_list)))
        date12_to_drop += tempList

        msg += '({})'.format(len(tempList))
        if len(tempList) <= 200:
            msg += '\n{}'.format(tempList)
        print(msg)

    # area ratio file
    if inps.areaRatioBased:
        print('--------------------------------------------------')
        print('use area-ratio-based network modification')

        # get area of interest for coherence calculation
        pix_box = get_aoi_pix_box(obj.metadata, inps.lookupFile, inps.aoi_pix_box, inps.aoi_geo_box)

        # calculate average coherence in masked out areas as threshold
        meanMaskCoh = np.nanmean(ut.spatial_average(inps.file,
                                                    datasetName='coherence',
                                                    maskFile=inps.maskFile,
                                                    saveList=True,
                                                    reverseMask=True)[0])
        print(f'Average coherence of {inps.maskFile} reverse is {meanMaskCoh:.2f}')

        # calculate area-ratio with pixels greater than meanMaskCoh
        areaRatioList = ut.spatial_average(inps.file,
                                           datasetName='coherence',
                                           maskFile=inps.maskFile,
                                           box=pix_box,
                                           saveList=True,
                                           checkAoi=True,
                                           threshold=meanMaskCoh)[0]

        # get area-ratio-based network
        area_ratio_date12_list = list(np.array(date12ListAll)[np.array(areaRatioList) >= inps.minAreaRatio])

        # get MST network
        mst_date12_list, msg = get_mst_date12(inps.keepMinSpanTree, areaRatioList, date12ListAll, date12_to_drop,
                                              min_par=inps.minAreaRatio,
                                              par_name='coherent area ratio')

        # drop all dates (below area-ratio thresh AND not in MST)
        tempList = sorted(list(set(date12ListAll) - set(area_ratio_date12_list + mst_date12_list)))
        date12_to_drop += tempList

        msg += '({})'.format(len(tempList))
        if len(tempList) <= 200:
            msg += '\n{}'.format(tempList)
        print(msg)

    # Manually drop pairs
    if inps.manual:
        tempList = manual_select_pairs_to_remove(inps.file)
        if tempList is None:
            return None
        tempList = [i for i in tempList if i in date12ListAll]
        print('date12 selected to remove: ({})\n{}'.format(len(tempList), tempList))
        date12_to_drop += tempList

    ## summary
    # drop duplicate date12 and sort in order
    date12_to_drop = sorted(list(set(date12_to_drop)))
    date12_to_keep = sorted(list(set(date12ListAll) - set(date12_to_drop)))
    print('--------------------------------------------------')
    print('number of interferograms to remove: {}'.format(len(date12_to_drop)))
    print('number of interferograms to keep  : {}'.format(len(date12_to_keep)))

    # print list of date to drop
    date_to_keep = [d for date12 in date12_to_keep for d in date12.split('_')]
    date_to_keep = sorted(list(set(date_to_keep)))
    date_to_drop = sorted(list(set(dateList) - set(date_to_keep)))
    if len(date_to_drop) > 0:
        print('number of acquisitions to remove: {}\n{}'.format(len(date_to_drop), date_to_drop))

    # checking:
    # 1) no new date12 to drop against existing file
    # 2) no date12 left after dropping
    date12ListKept = obj.get_date12_list(dropIfgram=True)
    date12ListDropped = sorted(list(set(date12ListAll) - set(date12ListKept)))
    if date12_to_drop == date12ListDropped:
        print('Calculated date12 to drop is the same as exsiting marked input file, skip updating file.')
        date12_to_drop = None
    elif date12_to_drop == date12ListAll:
        raise Exception('Zero interferogram left! Please adjust your setting and try again.')

    return date12_to_drop
