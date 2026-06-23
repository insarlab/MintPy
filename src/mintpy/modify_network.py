############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


# h5py 用来直接读写 HDF5 文件；这里会直接修改 ifgramStack.h5 里的 dropIfgram 数据集。
import h5py
# numpy 用来做数组计算和布尔筛选，例如按阈值选出要删除的干涉图。
import numpy as np
# matplotlib 用来交互式画图；manual_select_pairs_to_remove() 里用户可以点击网络图手动删边。
from matplotlib import dates as mdates, pyplot as plt

# ifgramStack 是 MintPy 表示干涉图栈 HDF5 文件的对象。
from mintpy.objects import ifgramStack
# network(pnet) 负责干涉网算法；plot(pp) 负责绘图；ptime 处理日期；ut 是通用工具函数。
from mintpy.utils import network as pnet, plot as pp, ptime, utils as ut


###########################  Sub Function  #############################
def reset_network(stackFile):
    """Reset/restore all pairs within the input file by set all DROP_IFGRAM=no"""
    # 这个函数用于“恢复网络”：把所有干涉图都标记为保留。
    # stackFile 通常是 ./inputs/ifgramStack.h5。
    print("reset dataset 'dropIfgram' to True for all interferograms for file: "+stackFile)
    # 创建 ifgramStack 对象，并打开 HDF5 文件读取元数据和 dropIfgram 数组。
    obj = ifgramStack(stackFile)
    obj.open(print_msg=False)
    if np.all(obj.dropIfgram):
        print('All dropIfgram are already True, no need to reset.')
    else:
        # h5py.File(..., 'r+') 表示以可读写模式打开 HDF5 文件。
        with h5py.File(stackFile, 'r+') as f:
            # dropIfgram 是布尔数组；True 表示该干涉图参与后续反演，False 表示丢弃。
            f['dropIfgram'][:] = True
        # touch() 更新 coherenceSpatialAvg.txt 的修改时间，让后续步骤知道网络状态变了。
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
    # 计算鼠标点击点 (x, y) 和所有网络节点 (x_array, y_array) 的欧氏距离。
    dist = np.sqrt((x_array - x)**2 + (y_array - y)**2)
    # np.argmin() 返回最小距离所在位置，也就是离鼠标点击最近的日期节点。
    idx = np.argmin(dist)
    #idx = dist==np.min(dist)
    return idx


def manual_select_pairs_to_remove(stackFile):
    """Manually select interferograms to remove"""
    # 这个函数会打开交互式网络图，让用户用鼠标点击两个日期来选择一条要删除的干涉边。
    print('\n-------------------------------------------------------------')
    print('Manually select interferograms to remove')
    print('1) click two dates/points to select one pair of interferogram')
    print('2) repeat until you select all pairs you would like to remove')
    print('3) close the figure to continue the program ...')
    print('-------------------------------------------------------------\n')
    obj = ifgramStack(stackFile)
    obj.open()
    # date12ListAll 是所有干涉图日期对，例如 20180101_20180113。
    date12ListAll = obj.date12List
    # pbase 是每个日期对应的垂直基线时间序列，用来绘制干涉网络纵轴。
    pbase = obj.get_perp_baseline_timeseries(dropIfgram=False)
    dateList = obj.dateList
    # matplotlib 内部用浮点数表示日期，date2num() 把日期转换成可画图的数字。
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
        # onclick 是鼠标点击事件回调函数：用户每点一次图，matplotlib 会调用它。
        idx = nearest_neighbor(event.xdata, event.ydata, datesNum, pbase)
        print('click at '+dateList[idx])
        date_click.append(dateList[idx])
        if len(date_click) % 2 == 0 and date_click[-2] != date_click[-1]:
            # 每点击两个不同日期，就组成一个干涉图日期对。
            [mDate, sDate] = sorted(date_click[-2:])
            mIdx = dateList.index(mDate)
            sIdx = dateList.index(sDate)
            date12 = mDate+'_'+sDate
            if date12 in date12ListAll:
                print('select date12: '+date12)
                date12_click.append(date12)
                # 在图上用红色粗线标出用户选择删除的干涉边。
                ax.plot([datesNum[mIdx], datesNum[sIdx]], [pbase[mIdx], pbase[sIdx]], 'r', lw=4)
            else:
                print(date12+' does not exist in input file')
        plt.draw()

    fig.canvas.mpl_connect('button_press_event', onclick)
    # plt.show() 打开图形窗口；用户关闭窗口后代码继续往下执行。
    plt.show()

    # yes_or_no() 会询问用户是否确认删除刚才选中的干涉图。
    if not ut.yes_or_no('Proceed to drop the ifgrams/date12?'):
        date12_click = None

    return date12_click


def get_aoi_pix_box(meta, lookup_file, pix_box, geo_box):
    """Get pix_box for AOI."""
    # AOI 是 area of interest，感兴趣区域。这个函数把用户输入的区域统一转换成像素框。
    coord = ut.coordinate(meta, lookup_file=lookup_file)

    # geo_box -> pix_box
    if geo_box and lookup_file:
        # 如果用户输入的是经纬度范围，并且有查找表，就转换成雷达坐标像素范围。
        print(f'input AOI in (lon0, lat1, lon1, lat0): {geo_box}')
        pix_box = coord.bbox_geo2radar(geo_box)

    # check pix_box
    if pix_box:
        # check_box_within_data_coverage() 保证像素框不会超出数据覆盖范围。
        pix_box = coord.check_box_within_data_coverage(pix_box)
        print(f'input AOI in (x0,y0,x1,y1): {pix_box}')

    return pix_box


def get_mst_date12(keep_mst, par_list_all, date12_list_all, date12_to_drop, min_par, par_name='average coherence'):
    """Get the date12_list of the MST network for the given parameter."""
    # MST 是 minimum spanning tree，最小生成树。
    # 在干涉网里保留 MST 可以避免网络被删断，保证时间序列反演仍有基本连通性。
    if keep_mst:
        print(f'Get minimum spanning tree (MST) of interferograms with inverse of {par_name}.')
        msg = ('Drop ifgrams with '
               '1) {} < {} AND '
               '2) not in MST network: '.format(par_name, min_par))

        # get the current remaining network (after all the above criteria and before data-driven)
        # 先从全部日期对里去掉前面规则已经决定删除的日期对。
        date12_to_keep = sorted(list(set(date12_list_all) - set(date12_to_drop)))
        par_to_keep = [par for par, date12 in zip(par_list_all, date12_list_all)
                       if date12 in date12_to_keep]

        # get MST from the current remaining network
        # threshold_coherence_based_mst() 根据参数反比构建 MST，返回需要保留的日期对。
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
    # 这是本文件最核心的函数：综合所有筛选规则，计算哪些 date12 干涉图要删除。
    obj = ifgramStack(inps.file)
    obj.open()
    date12ListAll = obj.date12List
    dateList = obj.dateList
    print(f'number of interferograms: {len(date12ListAll)}')

    # Get date12_to_drop
    # date12_to_drop 是待删除干涉图日期对列表，后面每条规则都会往里面追加。
    date12_to_drop = []

    # reference file
    if inps.referenceFile:
        # 如果提供参考文件，就只保留参考文件中仍然保留的干涉图。
        date12_to_keep = pnet.get_date12_list(inps.referenceFile, dropIfgram=True)
        print('--------------------------------------------------')
        print(f'use reference pairs info from file: {inps.referenceFile}')
        print(f'number of interferograms in reference: {len(date12_to_keep)}')
        tempList = sorted(list(set(date12ListAll) - set(date12_to_keep)))
        date12_to_drop += tempList
        print(f'date12 not in reference file: ({len(tempList)})\n{tempList}')

    # temp baseline threshold
    if inps.tempBaseMax:
        # tbaseIfgram 是时间基线；超过阈值的干涉图会被删除。
        tempIndex = np.abs(obj.tbaseIfgram) > inps.tempBaseMax
        tempList = list(np.array(date12ListAll)[tempIndex])
        date12_to_drop += tempList
        print('--------------------------------------------------')
        print('Drop ifgrams with temporal baseline > {} days: ({})\n{}'.format(
            inps.tempBaseMax, len(tempList), tempList))

    # perp baseline threshold
    if inps.perpBaseMax:
        # pbaseIfgram 是垂直空间基线；超过阈值的干涉图会被删除。
        tempIndex = np.abs(obj.pbaseIfgram) > inps.perpBaseMax
        tempList = list(np.array(date12ListAll)[tempIndex])
        date12_to_drop += tempList
        print('--------------------------------------------------')
        print('Drop ifgrams with perp baseline > {} meters: ({})\n{}'.format(
            inps.perpBaseMax, len(tempList), tempList))

    # connection number threshold
    if inps.connNumMax:
        # connNumMax 控制每个日期最多连接多少个相邻日期，常用于构建较稀疏的顺序网络。
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
        # excludeIfgIndex 按干涉图在文件里的索引编号删除。
        tempList = [date12ListAll[i] for i in inps.excludeIfgIndex]
        date12_to_drop += tempList
        print('--------------------------------------------------')
        print(f'Drop ifgrams with the following index number: {len(tempList)}')
        for ifg_idx, date12 in zip(inps.excludeIfgIndex, tempList):
            print(f'{ifg_idx} : {date12}')

    # excludeDate12
    if inps.excludeDate12:
        # excludeDate12 直接按日期对字符串删除，例如 20180101_20180113。
        tempList = [i for i in inps.excludeDate12 if i in date12ListAll]
        date12_to_drop += tempList
        print('--------------------------------------------------')
        print(f'Drop ifgrams with the following date12: {len(tempList)}')
        for date12 in tempList:
            print(date12)

    # excludeDate
    if inps.excludeDate:
        # excludeDate 删除所有包含指定单景日期的干涉图。
        tempList = [i for i in date12ListAll if any(j in inps.excludeDate for j in i.split('_'))]
        date12_to_drop += tempList
        print('-'*50+'\nDrop ifgrams including the following dates: ({})\n{}'.format(
            len(tempList), inps.excludeDate))
        print('-'*30+f'\n{tempList}')

    # startDate
    if inps.startDate:
        # startDate 删除早于指定日期的干涉图。
        minDate = int(inps.startDate)
        tempList = [i for i in date12ListAll if any(int(j) < minDate for j in i.split('_'))]
        date12_to_drop += tempList
        print('--------------------------------------------------')
        print('Drop ifgrams with date earlier than: {} ({})\n{}'.format(
            inps.startDate, len(tempList), tempList))

    # endDate
    if inps.endDate:
        # endDate 删除晚于指定日期的干涉图。
        maxDate = int(inps.endDate)
        tempList = [i for i in date12ListAll if any(int(j) > maxDate for j in i.split('_'))]
        date12_to_drop += tempList
        print('--------------------------------------------------')
        print('Drop ifgrams with date later than: {} ({})\n{}'.format(
            inps.endDate, len(tempList), tempList))

    # coherence file
    if inps.coherenceBased:
        # coherenceBased 根据平均相干性筛选低质量干涉图。
        print('--------------------------------------------------')
        print('use coherence-based network modification')

        # get area of interest for coherence calculation
        pix_box = get_aoi_pix_box(obj.metadata, inps.lookupFile, inps.aoiYX, inps.aoiLALO)

        # calculate spatial average coherence
        # ut.spatial_average() 计算每个干涉图在 AOI 或掩膜区域内的空间平均值。
        cohList = ut.spatial_average(inps.file,
                                     datasetName='coherence',
                                     maskFile=inps.maskFile,
                                     box=pix_box,
                                     saveList=True)[0]

        # get coherence-based network
        # 只保留平均相干性大于等于阈值的日期对。
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
        # areaRatioBased 根据“高相干像元面积比例”筛选干涉图。
        print('--------------------------------------------------')
        print('use area-ratio-based network modification')

        # get area of interest for coherence calculation
        pix_box = get_aoi_pix_box(obj.metadata, inps.lookupFile, inps.aoiYX, inps.aoiLALO)

        # calculate average coherence in masked out areas as threshold
        # reverseMask=True 表示在掩膜外区域计算平均相干性，作为动态阈值参考。
        meanMaskCoh = np.nanmean(ut.spatial_average(inps.file,
                                                    datasetName='coherence',
                                                    maskFile=inps.maskFile,
                                                    saveList=True,
                                                    reverseMask=True)[0])
        print(f'Average coherence of {inps.maskFile} reverse is {meanMaskCoh:.2f}')

        # calculate area-ratio with pixels greater than meanMaskCoh
        # threshold=meanMaskCoh 表示统计大于这个阈值的像元比例。
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
        # manual=True 时进入交互式手动选边删除流程。
        tempList = manual_select_pairs_to_remove(inps.file)
        if tempList is None:
            return None
        tempList = [i for i in tempList if i in date12ListAll]
        print(f'date12 selected to remove: ({len(tempList)})\n{tempList}')
        date12_to_drop += tempList

    ## summary
    # drop duplicate date12 and sort in order
    # set() 去重，sorted() 排序，得到最终稳定的删除列表。
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
    # 如果计算结果和当前文件已有 dropIfgram 状态完全一致，就不需要重复写文件。
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
    # modify_network() 是核心入口：计算要删除的干涉图，并写回 ifgramStack.h5 的 dropIfgram。

    if inps.reset:
        # reset 模式直接恢复所有干涉图，不执行其它筛选规则。
        print('--------------------------------------------------')
        reset_network(inps.file)
        return

    # get the list of date12 to drop/exclude
    # get_date12_to_drop() 返回 None 表示没有变化；返回 [] 表示保留全部；返回列表表示要删除这些日期对。
    inps.date12_to_drop = get_date12_to_drop(inps)

    # update date dataset in the ifgram stack h5 file
    if inps.date12_to_drop is not None:
        # update_drop_ifgram() 会把指定日期对的 dropIfgram 标记改为 False。
        ifgramStack(inps.file).update_drop_ifgram(inps.date12_to_drop)
        # 更新辅助文件时间戳，提示后续网络图/平均相干性可能需要重新生成。
        ut.touch('coherenceSpatialAvg.txt')
        print('Done.')

    return
