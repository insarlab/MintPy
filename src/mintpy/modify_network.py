############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import h5py
import numpy as np
from matplotlib import dates as mdates, pyplot as plt

from mintpy.objects import ifgramStack
from mintpy.utils import network as pnet, plot as pp, ptime, utils as ut


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
                print(date12+' does not exist in input file')
        plt.draw()

    fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()

    if not ut.yes_or_no('Proceed to drop the ifgrams/date12?'):
        date12_click = None

    return date12_click


def get_aoi_pix_box(meta, lookup_file, pix_box, geo_box):
    """Get pix_box for AOI."""
    coord = ut.coordinate(meta, lookup_file=lookup_file)

    # geo_box -> pix_box
    if geo_box and lookup_file:
        print(f'input AOI in (lon0, lat1, lon1, lat0): {geo_box}')
        pix_box = coord.bbox_geo2radar(geo_box)

    # check pix_box
    if pix_box:
        pix_box = coord.check_box_within_data_coverage(pix_box)
        print(f'input AOI in (x0,y0,x1,y1): {pix_box}')

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
        msg = f'Drop ifgrams with {par_name} < {min_par}: '
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
    print(f'number of interferograms: {len(date12ListAll)}')

    # Get date12_to_drop
    date12_to_drop = []

    # reference file
    if inps.referenceFile:
        date12_to_keep = pnet.get_date12_list(inps.referenceFile, dropIfgram=True)
        print('--------------------------------------------------')
        print(f'use reference pairs info from file: {inps.referenceFile}')
        print(f'number of interferograms in reference: {len(date12_to_keep)}')
        tempList = sorted(list(set(date12ListAll) - set(date12_to_keep)))
        date12_to_drop += tempList
        print(f'date12 not in reference file: ({len(tempList)})\n{tempList}')

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
            msg += f'\n{tempList}'
        print(msg)

    # excludeIfgIndex
    if inps.excludeIfgIndex:
        tempList = [date12ListAll[i] for i in inps.excludeIfgIndex]
        date12_to_drop += tempList
        print('--------------------------------------------------')
        print(f'Drop ifgrams with the following index number: {len(tempList)}')
        for ifg_idx, date12 in zip(inps.excludeIfgIndex, tempList):
            print(f'{ifg_idx} : {date12}')

    # excludeDate
    if inps.excludeDate:
        tempList = [i for i in date12ListAll if any(j in inps.excludeDate for j in i.split('_'))]
        date12_to_drop += tempList
        print('-'*50+'\nDrop ifgrams including the following dates: ({})\n{}'.format(
            len(tempList), inps.excludeDate))
        print('-'*30+f'\n{tempList}')

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
        pix_box = get_aoi_pix_box(obj.metadata, inps.lookupFile, inps.aoiYX, inps.aoiLALO)

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

        msg += f'({len(tempList)})'
        if len(tempList) <= 200:
            msg += f'\n{tempList}'
        print(msg)

    # area ratio file
    if inps.areaRatioBased:
        print('--------------------------------------------------')
        print('use area-ratio-based network modification')

        # get area of interest for coherence calculation
        pix_box = get_aoi_pix_box(obj.metadata, inps.lookupFile, inps.aoiYX, inps.aoiLALO)

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

        msg += f'({len(tempList)})'
        if len(tempList) <= 200:
            msg += f'\n{tempList}'
        print(msg)

    # Manually drop pairs
    if inps.manual:
        tempList = manual_select_pairs_to_remove(inps.file)
        if tempList is None:
            return None
        tempList = [i for i in tempList if i in date12ListAll]
        print(f'date12 selected to remove: ({len(tempList)})\n{tempList}')
        date12_to_drop += tempList

    ## summary
    # drop duplicate date12 and sort in order
    date12_to_drop = sorted(list(set(date12_to_drop)))
    date12_to_keep = sorted(list(set(date12ListAll) - set(date12_to_drop)))
    print('--------------------------------------------------')
    print(f'number of interferograms to remove: {len(date12_to_drop)}')
    print(f'number of interferograms to keep  : {len(date12_to_keep)}')

    # print list of date to drop
    date_to_keep = [d for date12 in date12_to_keep for d in date12.split('_')]
    date_to_keep = sorted(list(set(date_to_keep)))
    date_to_drop = sorted(list(set(dateList) - set(date_to_keep)))
    if len(date_to_drop) > 0:
        print(f'number of acquisitions to remove: {len(date_to_drop)}\n{date_to_drop}')

    # checking:
    # 1) no new date12 to drop against existing file
    # 2) no date12 left after dropping
    date12ListKept = obj.get_date12_list(dropIfgram=True)
    date12ListDropped = sorted(list(set(date12ListAll) - set(date12ListKept)))
    if date12_to_drop == date12ListDropped:
        print('Calculated date12 to drop is the same as existing marked input file, skip updating file.')
        date12_to_drop = None
    elif date12_to_drop == date12ListAll:
        raise Exception('Zero interferogram left! Please adjust your setting and try again.')

    return date12_to_drop


########################################################################
def modify_network(inps):
    """Run network modification."""

    if inps.reset:
        print('--------------------------------------------------')
        reset_network(inps.file)
        return

    # get the list of date12 to drop/exclude
    inps.date12_to_drop = get_date12_to_drop(inps)

    # update date dataset in the ifgram stack h5 file
    if inps.date12_to_drop is not None:
        ifgramStack(inps.file).update_drop_ifgram(inps.date12_to_drop)
        ut.touch('coherenceSpatialAvg.txt')
        print('Done.')

    return
