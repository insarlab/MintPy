############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


# glob 用来按通配符搜索文件，例如 '/path/*.unw' 会找到所有 .unw 文件。
import glob
# importlib 用来按字符串动态导入模块；这里根据 processor 名称导入 prep_isce/prep_aria 等模块。
import importlib
# os 用来处理路径、当前目录、判断文件是否存在等操作。
import os
# time 用来统计 load_data 整体运行耗时。
import time
# warnings 用来打印警告信息；警告不会立刻中断程序，和 raise Exception 不同。
import warnings

# subset 模块负责从模板中读取裁剪范围，例如按像素范围或经纬度范围裁剪数据。
from mintpy import subset
# auto_path 模块保存不同 InSAR 处理软件的默认路径规则，用于把模板里的 auto 转换成真实路径。
from mintpy.defaults import auto_path
# 这些对象和常量定义了 MintPy HDF5 文件中常见的数据集名称，以及几何/干涉图栈的数据结构。
from mintpy.objects import (
    GEOMETRY_DSET_NAMES,
    IFGRAM_DSET_NAMES,
    geometry,
    ifgramStack,
    sensor,
)
# geometryDict、ifgramDict、ifgramStackDict 是临时字典对象：
# 它们先收集“外部文件路径 + 元数据”，最后统一写成 MintPy 标准 HDF5 文件。
from mintpy.objects.stackDict import geometryDict, ifgramDict, ifgramStackDict
# ptime 处理日期格式；readfile 读取模板/属性/数据；ut 是 MintPy 通用工具函数集合。
from mintpy.utils import ptime, readfile, utils as ut

#################################################################
# PROCESSOR_LIST 列出 load_data 支持的上游处理软件名称。
# 用户模板中的 mintpy.load.processor 必须是这些名字之一。
PROCESSOR_LIST = ['isce', 'aria', 'hyp3', 'gmtsar', 'snap', 'gamma', 'roipac', 'cosicorr', 'nisar']

# primary observation dataset names
# OBS_DSET_NAMES 是主要观测数据集名称：相位、距离向偏移、方位向偏移。
OBS_DSET_NAMES = ['unwrapPhase', 'rangeOffset', 'azimuthOffset']

# 下面几个 *_DSET_NAME2TEMPLATE_KEY 字典是“数据集名称 -> 模板配置项”的对应表。
# 例如 unwrapPhase 在 HDF5 里叫 unwrapPhase，但它的输入文件路径来自模板项 mintpy.load.unwFile。
IFG_DSET_NAME2TEMPLATE_KEY = {
    'unwrapPhase'     : 'mintpy.load.unwFile',
    'coherence'       : 'mintpy.load.corFile',
    'connectComponent': 'mintpy.load.connCompFile',
    'wrapPhase'       : 'mintpy.load.intFile',
    'magnitude'       : 'mintpy.load.magFile',
}

ION_DSET_NAME2TEMPLATE_KEY = {
    'unwrapPhase'     : 'mintpy.load.ionUnwFile',
    'coherence'       : 'mintpy.load.ionCorFile',
    'connectComponent': 'mintpy.load.ionConnCompFile',
}

OFF_DSET_NAME2TEMPLATE_KEY = {
    'azimuthOffset'   : 'mintpy.load.azOffFile',
    'azimuthOffsetStd': 'mintpy.load.azOffStdFile',
    'rangeOffset'     : 'mintpy.load.rgOffFile',
    'rangeOffsetStd'  : 'mintpy.load.rgOffStdFile',
    'offsetSNR'       : 'mintpy.load.offSnrFile',
}

GEOM_DSET_NAME2TEMPLATE_KEY = {
    'height'          : 'mintpy.load.demFile',
    'latitude'        : 'mintpy.load.lookupYFile',
    'longitude'       : 'mintpy.load.lookupXFile',
    'azimuthCoord'    : 'mintpy.load.lookupYFile',
    'rangeCoord'      : 'mintpy.load.lookupXFile',
    'incidenceAngle'  : 'mintpy.load.incAngleFile',
    'azimuthAngle'    : 'mintpy.load.azAngleFile',
    'shadowMask'      : 'mintpy.load.shadowMaskFile',
    'waterMask'       : 'mintpy.load.waterMaskFile',
    'bperp'           : 'mintpy.load.bperpFile',
}


#################################################################
def read_inps2dict(inps):
    """Read input Namespace object info into iDict.

    It grab the following contents into iDict
    1. inps & all template files
    2. configurations: processor, autoPath, updateMode, compression, x/ystep
    3. extra metadata: PLATFORM, PROJECT_NAME,
    4. translate autoPath

    Parameters: inps  - namespace, input arguments from command line & template file
    Returns:    iDict - dict,      input arguments from command line & template file
    """
    # Read input info into iDict
    # vars(inps) 把 argparse.Namespace 对象转换成普通字典，方便用 iDict['key'] 读写。
    iDict = vars(inps)
    # 先给一些关键字段设置默认值，后面模板或项目名识别会覆盖它们。
    iDict['PLATFORM'] = None
    iDict['processor'] = 'isce'

    # Read template file
    template = {}
    for fname in inps.template_file:
        # read_template() 把 cfg 文件读成字典；check_template_auto_value() 会处理 yes/no/auto 等值。
        temp = readfile.read_template(fname)
        temp = ut.check_template_auto_value(temp)
        # update() 会把当前模板内容合并到 template；后读的模板可以覆盖先读的同名配置。
        template.update(temp)
    for key, value in template.items():
        iDict[key] = value

    # group - load
    prefix = 'mintpy.load.'
    # 找出所有以 mintpy.load. 开头的模板配置，并去掉前缀得到短 key。
    key_list = [i.split(prefix)[1] for i in template.keys() if i.startswith(prefix)]
    for key in key_list:
        value = template[prefix+key]
        if key in ['processor', 'autoPath', 'updateMode', 'compression']:
            # 这些配置在后面经常使用，所以同时保存成短 key，例如 iDict['processor']。
            iDict[key] = template[prefix+key]
        elif value:
            # 其它 load 配置保留完整 key，避免和其它配置冲突。
            iDict[prefix+key] = template[prefix+key]
    print('processor : {}'.format(iDict['processor']))

    # group - multilook
    prefix = 'mintpy.multilook.'
    key_list = [i.split(prefix)[1] for i in template.keys() if i.startswith(prefix)]
    for key in key_list:
        value = template[prefix+key]
        if key in ['xstep', 'ystep', 'method']:
            iDict[key] = template[prefix+key]

    iDict['xstep']  = int(iDict.get('xstep', 1))
    iDict['ystep']  = int(iDict.get('ystep', 1))
    iDict['method'] = str(iDict.get('method', 'nearest'))

    # PROJECT_NAME --> PLATFORM
    if not iDict['PROJECT_NAME']:
        # 如果命令行没有给项目名，就尝试从自定义模板文件名中推断项目名和卫星平台。
        cfile = [i for i in list(inps.template_file) if os.path.basename(i) != 'smallbaselineApp.cfg']
        iDict['PROJECT_NAME'] = sensor.project_name2sensor_name(cfile)[1]

    msg = 'SAR platform/sensor : '
    sensor_name = sensor.project_name2sensor_name(str(iDict['PROJECT_NAME']))[0]
    if sensor_name:
        msg += str(sensor_name)
        iDict['PLATFORM'] = str(sensor_name)
    else:
        msg += 'unknown from project name "{}"'.format(iDict['PROJECT_NAME'])
    print(msg)

    # update file path with auto
    if iDict.get('autoPath', False):
        print('use auto path defined in mintpy.defaults.auto_path for options in auto')
        # auto_path.get_auto_path() 会把模板里的 auto 路径替换成当前 processor 的默认文件路径模式。
        iDict = auto_path.get_auto_path(processor=iDict['processor'],
                                        work_dir=os.getcwd(),
                                        template=iDict)

    return iDict


def read_subset_box(iDict):
    """read the following items:
    geocoded - bool, if the stack of observations geocoded or not
    box      - tuple of 4 int, pixel box for stackObj and geomRadarObj, for obs in geo & radar coordinates
    box4geo  - tuple of 4 int, pixel box for geomGeoObj, box4geo is the same as box, except for:
               obs in radar coordinate with lookup table [for gamma and roipac], where box4geo is
               the geo bounding box of the box above.
    """
    # Read subset info from template
    # box 是用于普通观测/雷达几何文件的裁剪范围；box4geo 是用于地理坐标查找表的裁剪范围。
    iDict['box'] = None
    iDict['box4geo'] = None
    # subset.read_subset_template2box() 从模板读取裁剪设置，返回像素框 pix_box 和经纬度框 geo_box。
    pix_box, geo_box = subset.read_subset_template2box(iDict['template_file'][0])

    # Grab required info to read input geo_box into pix_box
    lookup_y_files = glob.glob(str(iDict['mintpy.load.lookupYFile']))
    lookup_x_files = glob.glob(str(iDict['mintpy.load.lookupXFile']))
    if len(lookup_y_files) > 0 and len(lookup_x_files) > 0:
        lookup_file = [lookup_y_files[0], lookup_x_files[0]]
    else:
        lookup_file = None

    # use DEM file attribute as reference, because
    # 1) it is required AND
    # 2) it is in the same coordinate type as observation files
    dem_files = glob.glob(iDict['mintpy.load.demFile'])
    if len(dem_files) > 0:
        # 用 DEM 的属性判断输入数据是地理坐标还是雷达坐标；Y_FIRST 是地理坐标文件常见属性。
        atr = readfile.read_attribute(dem_files[0])
    else:
        atr = dict()

    geocoded = True if 'Y_FIRST' in atr.keys() else False
    iDict['geocoded'] = geocoded

    # Check conflict
    if geo_box and not geocoded and lookup_file is None:
        geo_box = None
        print('WARNING: mintpy.subset.lalo is not supported'
              ' if 1) no lookup file AND'
              '    2) radar/unknown coded dataset')
        print('\tignore it and continue.')

    if not geo_box and not pix_box:
        # adjust for the size inconsistency problem in SNAP geocoded products
        # ONLY IF there is no input subset
        # Use the min bbox if files size are different
        if iDict['processor'] == 'snap':
            fnames = ut.get_file_list(iDict['mintpy.load.unwFile'])
            pix_box = update_box4files_with_inconsistent_size(fnames)

        if not pix_box:
            return iDict

    # geo_box --> pix_box
    # ut.coordinate() 创建坐标转换对象，可以在经纬度框和像素框之间转换。
    coord = ut.coordinate(atr, lookup_file=lookup_file)
    if geo_box is not None:
        # bbox_geo2radar() 把经纬度范围转换为雷达坐标下的像素范围。
        pix_box = coord.bbox_geo2radar(geo_box)
        pix_box = coord.check_box_within_data_coverage(pix_box)
        print(f'input bounding box of interest in lat/lon: {geo_box}')
    print(f'box to read for datasets in y/x: {pix_box}')

    # Get box for geocoded lookup table (for gamma/roipac)
    box4geo_lut = None
    if lookup_file is not None:
        atrLut = readfile.read_attribute(lookup_file[0])
        if not geocoded and 'Y_FIRST' in atrLut.keys():
            # 对于 Gamma/ROI_PAC 等情况，查找表可能是地理坐标，需要单独计算它的裁剪框。
            geo_box = coord.bbox_radar2geo(pix_box)
            box4geo_lut = ut.coordinate(atrLut).bbox_geo2radar(geo_box)
            print(f'box to read for geocoded lookup file in y/x: {box4geo_lut}')

    iDict['box'] = pix_box
    iDict['box4geo'] = box4geo_lut if box4geo_lut else pix_box
    return iDict


#################################################################
def update_box4files_with_inconsistent_size(fnames):
    """Check the size (row / column number) of a list of files
    For SNAP geocoded products has one line missing in some interferograms, Andre, 2019-07-16
    Parameters: fnames  - list of path for interferogram files
    Returns:    pix_box - None if all files are in same size
                          (0, 0, min_width, min_length) if not.
    """
    # 读取每个输入干涉图的 LENGTH/WIDTH 属性，检查行列数是否一致。
    atr_list = [readfile.read_attribute(fname) for fname in fnames]
    length_list = [int(atr['LENGTH']) for atr in atr_list]
    width_list = [int(atr['WIDTH']) for atr in atr_list]
    if any(len(set(i)) for i in [length_list, width_list]):
        # 如果尺寸不一致，选择所有文件都共有的最小范围，避免读取越界。
        min_length = min(length_list)
        min_width = min(width_list)
        pix_box = (0, 0, min_width, min_length)

        # print out warning message
        msg = '\n'+'*'*80
        msg += '\nWARNING: NOT all input unwrapped interferograms have the same row/column number!'
        msg += f'\nMinimum size is: ({min_length}, {min_width})'
        msg += '\n'+'-'*30
        msg += '\nThe following dates have different size:'

        for i in range(len(fnames)):
            if length_list[i] != min_length or width_list[i] != min_width:
                msg += '\n\t{}\t({}, {})'.format(atr_list[i]['DATE12'], length_list[i], width_list[i])

        msg += '\n'+'-'*30
        msg += '\nAssuming the interferograms above have:'
        msg += '\n\textra line(s) at the bottom OR'
        msg += '\n\textra column(s) at the right'
        msg += '\nContinue to load data using subset of the minimum size.'
        msg += '\n'+'*'*80+'\n'
        print(msg)
    else:
        pix_box = None
    return pix_box


def skip_files_with_inconsistent_size(dsPathDict, pix_box=None, dsName='unwrapPhase'):
    """Skip files by removing the file path from the input dsPathDict."""
    # dsPathDict 保存每类数据集对应的文件路径列表，例如 {'unwrapPhase': [...], 'coherence': [...]}。
    atr_list = [readfile.read_attribute(fname) for fname in dsPathDict[dsName]]
    length_list = [int(atr['LENGTH']) for atr in atr_list]
    width_list = [int(atr['WIDTH']) for atr in atr_list]

    # Check size requirements
    drop_inconsistent_files = False
    if any(len(set(size_list)) > 1 for size_list in [length_list, width_list]):
        # 如果没有指定裁剪框，尺寸不同的文件无法放进同一个 HDF5 栈，只能跳过异常尺寸文件。
        if pix_box is None:
            drop_inconsistent_files = True
        else:
            # if input subset is within the min file sizes: do NOT drop
            max_box_width, max_box_length = pix_box[2:4]
            if max_box_length > min(length_list) or max_box_width > min(width_list):
                drop_inconsistent_files = True

    # update dsPathDict
    if drop_inconsistent_files:
        # most_common() 找出最常见的行列数；非这个尺寸的干涉图会被移除。
        common_length = ut.most_common(length_list)
        common_width = ut.most_common(width_list)

        # print out warning message
        msg = '\n'+'*'*80
        msg += '\nWARNING: NOT all input unwrapped interferograms have the same row/column number!'
        msg += f'\nThe most common size is: ({common_length}, {common_width})'
        msg += '\n'+'-'*30
        msg += '\nThe following dates have different size:'

        dsNames = list(dsPathDict.keys())
        date12_list = [atr['DATE12'] for atr in atr_list]
        num_drop = 0
        for date12, length, width in zip(date12_list, length_list, width_list):
            if length != common_length or width != common_width:
                dates = ptime.yyyymmdd(date12.split('-'))
                # update file list for all datasets
                # 对某一个异常 date12，要从 unwrapPhase/coherence/connectComponent 等所有数据集中一起移除。
                for dsName in dsNames:
                    fnames = [i for i in dsPathDict[dsName]
                              if all(d[2:8] in i for d in dates)]
                    if len(fnames) > 0:
                        dsPathDict[dsName].remove(fnames[0])
                msg += f'\n\t{date12}\t({length}, {width})'
                num_drop += 1

        msg += '\n'+'-'*30
        msg += f'\nSkip loading the above interferograms ({num_drop}).'
        msg += f'\nContinue to load the rest interferograms ({len(date12_list) - num_drop}).'
        msg += '\n'+'*'*80+'\n'
        print(msg)
    return dsPathDict


def read_inps_dict2ifgram_stack_dict_object(iDict, ds_name2template_key):
    """Read input arguments into ifgramStackDict object.

    Parameters: iDict                - dict, input arguments from command line & template file
                ds_name2template_key - dict, to relate the HDF5 dataset name to the template key
    Returns:    stackObj             - ifgramStackDict object or None
    """
    if iDict['only_load_geometry']:
        # 用户指定 --geom 时，只加载几何文件，不创建干涉图/电离层/偏移栈对象。
        return None

    # 根据传入的映射表判断这次要加载的是普通干涉图栈、电离层栈还是偏移栈。
    if 'mintpy.load.unwFile' in ds_name2template_key.values():
        obs_type = 'interferogram'
    elif 'mintpy.load.ionUnwFile' in ds_name2template_key.values():
        obs_type = 'ionosphere'
    elif 'mintpy.load.azOffFile' in ds_name2template_key.values():
        obs_type = 'offset'

    # iDict --> dsPathDict
    print('-'*50)
    print(f'searching {obs_type} pairs info')
    print('input data files:')
    max_digit = max(len(i) for i in list(ds_name2template_key.keys()))
    dsPathDict = {}
    for dsName in [i for i in IFGRAM_DSET_NAMES if i in ds_name2template_key.keys()]:
        key = ds_name2template_key[dsName]
        if key in iDict.keys():
            # glob.glob() 按模板路径搜索真实文件；sorted() 保证顺序稳定。
            files = sorted(glob.glob(str(iDict[key])))
            if len(files) > 0:
                dsPathDict[dsName] = files
                print(f'{dsName:<{max_digit}}: {iDict[key]}')

    # Check 1: required dataset
    # 普通干涉图必须有 unwrapPhase；偏移栈必须至少有 rangeOffset/azimuthOffset 之一。
    dsName0s = [x for x in OBS_DSET_NAMES if x in ds_name2template_key.keys()]
    dsName0 = [i for i in dsName0s if i in dsPathDict.keys()]
    if len(dsName0) == 0:
        print(f'WARNING: No data files found for the required dataset: {dsName0s}! Skip loading for {obs_type} stack.')
        return None
    else:
        dsName0 = dsName0[0]

    # Check 2: data dimension for unwrapPhase files
    dsPathDict = skip_files_with_inconsistent_size(
        dsPathDict=dsPathDict,
        pix_box=iDict['box'],
        dsName=dsName0)

    # Check 3: number of files for all dataset types
    # dsPathDict --> dsNumDict
    dsNumDict = {}
    for key in dsPathDict.keys():
        num_file = len(dsPathDict[key])
        dsNumDict[key] = num_file
        print(f'number of {key:<{max_digit}}: {num_file}')

    dsNumList = list(dsNumDict.values())
    if any(i != dsNumList[0] for i in dsNumList):
        # 不同数据集数量不一致时继续运行，但后面会按 date12 匹配，跳过缺文件的干涉对。
        msg = 'WARNING: NOT all types of dataset have the same number of files.'
        msg += ' -> skip interferograms with missing files and continue.'
        print(msg)
        #raise Exception(msg)

    # dsPathDict --> pairsDict --> stackObj
    dsNameList = list(dsPathDict.keys())

    #####################################
    # A dictionary of data file paths for a list of pairs, e.g.:
    # pairsDict = {
    #     ('date1', 'date2') : ifgramPathDict1,
    #     ('date1', 'date3') : ifgramPathDict2,
    #     ...,
    # }

    pairsDict = {}
    for i, dsPath0 in enumerate(dsPathDict[dsName0]):
        # date string used in the file/dir path
        # YYYYDDD       for gmtsar [day of the year - 1]
        # YYYYMMDDTHHMM for uavsar
        # YYYYMMDD      for all the others
        date6s = readfile.read_attribute(dsPath0)['DATE12'].replace('_','-').split('-')
        if iDict['processor'] == 'gmtsar':
            date12MJD = os.path.basename(os.path.dirname(dsPath0))
        else:
            date12MJD = None

        #####################################
        # A dictionary of data file paths for a given pair.
        # One pair may have several types of dataset, e.g.:
        # ifgramPathDict1 = {
        #     'unwrapPhase': /dirPathToFile/filt_fine.unw,
        #     'coherence'  : /dirPathToFile/filt_fine.cor,
        #     ...
        # }
        # All path of data file must contain the reference and secondary date, in file/dir name.

        ifgramPathDict = {}
        for dsName in dsNameList:
            # search the matching data file for the given date12
            # 1st guess: file in the same order as the one for dsName0
            dsPath1 = dsPathDict[dsName][i]
            if (all(d6 in dsPath1 for d6 in date6s)
                    or (date12MJD and date12MJD in dsPath1)):
                ifgramPathDict[dsName] = dsPath1

            else:
                # 2nd guess: any file in the list
                # 如果同一索引位置不是对应日期，就在整个文件列表里搜索包含同一对日期的文件。
                dsPath2 = [p for p in dsPathDict[dsName]
                           if (all(d6 in p for d6 in date6s)
                                   or (date12MJD and date12MJD in dsPath1))]

                if len(dsPath2) > 0:
                    ifgramPathDict[dsName] = dsPath2[0]
                else:
                    print(f'WARNING: {dsName:>18} file missing for pair {date6s}')

        # initiate ifgramDict object
        # ifgramDict 表示一个干涉对的所有数据文件路径集合。
        ifgramObj = ifgramDict(datasetDict=ifgramPathDict)

        # update pairsDict object
        # pairsDict 的 key 是日期对，例如 ('20180101', '20180113')。
        date8s = ptime.yyyymmdd(date6s)
        pairsDict[tuple(date8s)] = ifgramObj

    if len(pairsDict) > 0:
        # ifgramStackDict 表示整个干涉图栈，后面可以一次性写入 ifgramStack.h5。
        stackObj = ifgramStackDict(pairsDict=pairsDict, dsName0=dsName0)
    else:
        stackObj = None
    return stackObj


def read_inps_dict2geometry_dict_object(iDict, dset_name2template_key):
    """Read input arguments into geometryDict object(s).

    Parameters: iDict        - dict, input arguments from command line & template file
    Returns:    geomGeoObj   - geometryDict object in geo   coordinates or None
                geomRadarObj - geometryDict object in radar coordinates or None
    """

    # eliminate lookup table dsName for input files in radar-coordinates
    # 不同上游处理软件输出的查找表坐标系不同，因此需要去掉不适用的数据集名称。
    if iDict['processor'] in ['isce', 'doris']:
        # for processors with lookup table in radar-coordinates, remove azimuth/rangeCoord
        dset_name2template_key.pop('azimuthCoord')
        dset_name2template_key.pop('rangeCoord')
    elif iDict['processor'] in ['roipac', 'gamma']:
        # for processors with lookup table in geo-coordinates, remove latitude/longitude
        dset_name2template_key.pop('latitude')
        dset_name2template_key.pop('longitude')
    elif iDict['processor'] in ['aria', 'gmtsar', 'hyp3', 'snap', 'cosicorr']:
        # for processors with geocoded products support only, do nothing for now.
        # check again when adding products support in radar-coordiantes
        pass
    else:
        print('Un-recognized InSAR processor: {}'.format(iDict['processor']))

    # iDict --> dsPathDict
    print('-'*50)
    print('searching geometry files info')
    print('input data files:')
    max_digit = max(len(i) for i in list(dset_name2template_key.keys()))
    dsPathDict = {}
    for dsName in [i for i in GEOMETRY_DSET_NAMES if i in dset_name2template_key.keys()]:
        key = dset_name2template_key[dsName]
        if key in iDict.keys():
            # 对每个几何数据集，按模板路径找到真实输入文件。
            files = sorted(glob.glob(str(iDict[key])))
            if len(files) > 0:
                if dsName == 'bperp':
                    # bperp 是垂直基线，通常每个日期一个文件，因此保存成 date -> file 的字典。
                    bperpDict = {}
                    for file in files:
                        date = ptime.yyyymmdd(os.path.basename(os.path.dirname(file)))
                        bperpDict[date] = file
                    dsPathDict[dsName] = bperpDict
                    print(f'{dsName:<{max_digit}}: {iDict[key]}')
                    print(f'number of bperp files: {len(list(bperpDict.keys()))}')
                else:
                    # 其它几何量通常只需要一个文件，例如 DEM、高程角、阴影掩膜。
                    dsPathDict[dsName] = files[0]
                    print(f'{dsName:<{max_digit}}: {files[0]}')

    # Check required dataset
    dsName0 = GEOMETRY_DSET_NAMES[0]
    if dsName0 not in dsPathDict.keys():
        print(f'WARNING: No reqired {dsName0} data files found!')

    # extra metadata from observations
    # e.g. EARTH_RADIUS, HEIGHT, etc.
    # 几何文件可能缺少部分元数据，因此从观测文件中补充一些通用元数据。
    obsMetaGeo = None
    obsMetaRadar = None
    for obsName in OBS_DSET_NAMES:
        obsFiles = sorted(glob.glob(iDict[dset_name2template_key[obsName]]))
        if len(obsFiles) > 0:
            atr = readfile.read_attribute(obsFiles[0])
            if 'Y_FIRST' in atr.keys():
                obsMetaGeo = atr.copy()
            else:
                obsMetaRadar = atr.copy()
            break

    # dsPathDict --> dsGeoPathDict + dsRadarPathDict
    # 按文件属性中是否存在 Y_FIRST，把几何文件分成地理坐标和雷达坐标两类。
    dsNameList = list(dsPathDict.keys())
    dsGeoPathDict = {}
    dsRadarPathDict = {}
    for dsName in dsNameList:
        if dsName == 'bperp':
            atr = readfile.read_attribute(next(iter(dsPathDict[dsName].values())))
        else:
            atr = readfile.read_attribute(dsPathDict[dsName])
        if 'Y_FIRST' in atr.keys():
            dsGeoPathDict[dsName] = dsPathDict[dsName]
        else:
            dsRadarPathDict[dsName] = dsPathDict[dsName]

    geomGeoObj = None
    geomRadarObj = None
    if len(dsGeoPathDict) > 0:
        # geometryDict 是“待写入 geometryGeo.h5 的几何数据集合”。
        geomGeoObj = geometryDict(
            processor=iDict['processor'],
            datasetDict=dsGeoPathDict,
            extraMetadata=obsMetaGeo)
    if len(dsRadarPathDict) > 0:
        # geometryDict 是“待写入 geometryRadar.h5 的几何数据集合”。
        geomRadarObj = geometryDict(
            processor=iDict['processor'],
            datasetDict=dsRadarPathDict,
            extraMetadata=obsMetaRadar)

    return geomGeoObj, geomRadarObj


#################################################################
def run_or_skip(outFile, inObj, box, updateMode=True, xstep=1, ystep=1, geom_obj=None):
    """Check if re-writing is necessary.

    Do not write HDF5 file if ALL the following meet:
        1. HDF5 file exists and is readable,
        2. HDF5 file contains all the datasets and in the same size
        3. For ifgramStackDict, HDF5 file contains all date12.

    Parameters: outFile    - str, path to the output HDF5 file
                inObj      - ifgramStackDict or geometryDict, object to write
                box        - tuple of int, bounding box in (x0, y0, x1, y1)
                updateMode - bool
                x/ystep    - int
                geom_obj   - geometryDict object or None, for ionosphere only
    Returns:    flag       - str, run or skip
    """

    flag = 'run'

    # skip if there is no dict object to write
    if not inObj:
        # 没有输入对象，说明没有对应数据需要写入。
        flag = 'skip'
        return flag

    # run if not in update mode
    if not updateMode:
        # updateMode 关闭时，总是重新写入。
        return flag

    if ut.run_or_skip(outFile, readable=True) == 'skip':
        # 只有输出文件已经存在且可读时，才进一步检查里面的数据集是否完整。
        kwargs = dict(box=box, xstep=xstep, ystep=ystep)

        if inObj.name == 'ifgramStack':
            # 对干涉图栈，要同时检查尺寸、数据集名称和 date12 日期对是否完整。
            in_size = inObj.get_size(geom_obj=geom_obj, **kwargs)[1:]
            in_dset_list = inObj.get_dataset_list()
            in_date12_list = inObj.get_date12_list()

            outObj = ifgramStack(outFile)
            outObj.open(print_msg=False)
            out_size = (outObj.length, outObj.width)
            out_dset_list = outObj.datasetNames
            out_date12_list = outObj.date12List

            if (out_size[1:] == in_size[1:]
                    and set(in_dset_list).issubset(set(out_dset_list))
                    and set(in_date12_list).issubset(set(out_date12_list))):
                print('All date12   exists in file {} with same size as required,'
                      ' no need to re-load.'.format(os.path.basename(outFile)))
                flag = 'skip'

        elif inObj.name == 'geometry':
            # 对几何文件，只需要检查尺寸和数据集名称是否完整。
            in_size = inObj.get_size(**kwargs)
            in_dset_list = inObj.get_dataset_list()

            outObj = geometry(outFile)
            outObj.open(print_msg=False)
            out_size = (outObj.length, outObj.width)
            out_dset_list = outObj.datasetNames

            if (out_size == in_size
                    and set(in_dset_list).issubset(set(out_dset_list))):
                print('All datasets exists in file {} with same size as required,'
                      ' no need to re-load.'.format(os.path.basename(outFile)))
                flag = 'skip'

    return flag


def prepare_metadata(iDict):
    """Prepare metadata via prep_{processor}.py scripts."""
    # metadata 是描述数据的重要属性，例如日期、轨道、波长、坐标范围等。
    # 不同处理器的原始产品格式不同，所以这里先调用对应 prep_xxx.py 标准化元数据。
    processor = iDict['processor']
    script_name = f'prep_{processor}.py'
    print('-'*50)
    print(f'prepare metadata files for {processor} products')

    if processor not in PROCESSOR_LIST:
        msg = f'un-recognized InSAR processor: {processor}'
        msg += f'\nsupported processors: {PROCESSOR_LIST}'
        raise ValueError(msg)

    # import prep_{processor}
    # importlib.import_module() 通过字符串动态导入模块，例如 processor='isce' 时导入 mintpy.cli.prep_isce。
    prep_module = importlib.import_module(f'mintpy.cli.prep_{processor}')

    if processor in ['gamma', 'hyp3', 'roipac', 'snap', 'cosicorr']:
        # run prep_module
        # 这些 processor 通常需要对每类输入文件逐个补充/转换元数据。
        for key in [i for i in iDict.keys()
                    if (i.startswith('mintpy.load.')
                        and i.endswith('File')
                        and i != 'mintpy.load.metaFile')]:
            if len(glob.glob(str(iDict[key]))) > 0:
                # print command line
                # iargs 模拟命令行参数，传给 prep_xxx.py 的 main()。
                iargs = [iDict[key]]
                if processor == 'gamma':
                    if iDict['PLATFORM']:
                        iargs += ['--sensor', iDict['PLATFORM'].lower()]
                    # add DEM file to faciliate the checking and metadata extraction for geocoded datasets
                    iargs += ['--dem', iDict['mintpy.load.demFile']]

                elif processor == 'cosicorr':
                    iargs += ['--metadata', iDict['mintpy.load.metaFile']]
                ut.print_command_line(script_name, iargs)
                # run
                prep_module.main(iargs)

    elif processor == 'nisar':
        # NISAR 的 GUNW 产品需要 DEM、频率等信息，且这里要求 DEM 必须真实存在。
        dem_file = iDict['mintpy.load.demFile']
        gunw_files = iDict['mintpy.load.unwFile']
        water_mask = iDict['mintpy.load.waterMaskFile']
        frequency = iDict.get('mintpy.load.frequency', 'auto')

        if str(dem_file).lower() in ['auto', 'none', 'no', '']:
            raise ValueError(
                'mintpy.load.demFile is required for processor=nisar. '
                'Please set it to a real DEM path in the template.'
            )
        dem_file = os.path.expanduser(str(dem_file))
        if not os.path.isfile(dem_file):
            raise FileNotFoundError(
                f'No DEM file found for mintpy.load.demFile: {dem_file}'
            )

        if len(glob.glob(str(gunw_files))) == 0:
            raise FileNotFoundError(
                f'No input GUNW files found for mintpy.load.unwFile: {gunw_files}'
            )

        # run prep_*.py
        # prep_nisar.py 会读取 GUNW 文件并准备 MintPy 后续加载所需的元数据。
        iargs = ['-i', gunw_files, '-d', dem_file, '--frequency', frequency]

        if str(water_mask).lower() not in ['auto', 'none', 'no', ''] and os.path.exists(water_mask):
            iargs = iargs + ['--mask', water_mask]

        if iDict['mintpy.subset.yx']:
            warnings.warn('Subset in Y/X is not implemented for NISAR. \n'
                          'There might be shift in the coordinates of different products. \n'
                          'Use lat/lon instead.')
        if iDict['mintpy.subset.lalo']:
            lalo = iDict['mintpy.subset.lalo'].split(',')
            lats = lalo[0].split(':')
            lons = lalo[1].split(':')
            iargs = iargs + ['--sub-lat', lats[0], lats[1], '--sub-lon', lons[0], lons[1]]

        ut.print_command_line(script_name, iargs)
        try:
            prep_module.main(iargs)
        except:
            # 这里捕获所有异常并给警告，是为了兼容“元数据已经存在”的情况继续往下走。
            warnings.warn('prep_nisar.py failed. Assuming its result exists and continue...')

    elif processor == 'isce':
        # ISCE 产品需要 meta 文件、baseline 目录、geometry 目录以及观测文件路径。
        from mintpy.utils import isce_utils, s1_utils

        # --meta-file
        meta_files = sorted(glob.glob(iDict['mintpy.load.metaFile']))
        if len(meta_files) > 0:
            meta_file = meta_files[0]
        else:
            warnings.warn('No input metadata file found: {}'.format(iDict['mintpy.load.metaFile']))
            meta_file = 'auto'

        # --baseline-dir / --geometry-dir
        baseline_dir = iDict['mintpy.load.baselineDir']
        geom_dir = os.path.dirname(iDict['mintpy.load.demFile'])
        geom_dir = os.path.abspath(geom_dir)

        # --dset-dir / --file-pattern
        obs_keys = [
            'mintpy.load.unwFile',
            'mintpy.load.ionUnwFile',
            'mintpy.load.rgOffFile',
            'mintpy.load.azOffFile',
        ]
        obs_paths = [iDict[key] for key in obs_keys if iDict[key].lower() != 'auto']
        obs_paths = [x for x in obs_paths if len(glob.glob(x)) > 0]

        # --geom-files for the basenames only
        geom_names = ['dem', 'lookupY', 'lookupX', 'incAngle', 'azAngle', 'shadowMask', 'waterMask']
        geom_keys = [f'mintpy.load.{i}File' for i in geom_names]
        geom_files = [os.path.basename(iDict[key]) for key in geom_keys
                      if iDict.get(key, 'auto') not in ['auto', 'None', 'no',  None, False]]

        # compose list of input arguments
        iargs = ['-m', meta_file, '-g', geom_dir]
        if baseline_dir:
            iargs += ['-b', baseline_dir]
        if len(obs_paths) > 0:
            iargs += ['-f'] + obs_paths
        if geom_files:
            iargs += ['--geom-files'] + geom_files

        # run module
        ut.print_command_line(script_name, iargs)
        try:
            prep_module.main(iargs)
        except:
            warnings.warn('prep_isce.py failed. Assuming its result exists and continue...')

        # [optional] for topsStack: SAFE_files.txt --> S1A/B_date.txt
        # Sentinel-1 TOPS 数据有时需要 SAFE_files.txt 来生成 S1A/B 日期列表，用于后续范围偏差修正等。
        if os.path.isfile(meta_file) and isce_utils.get_processor(meta_file) == 'topsStack':
            safe_list_file = os.path.join(os.path.dirname(os.path.dirname(meta_file)), 'SAFE_files.txt')
            if os.path.isfile(safe_list_file):
                s1_utils.get_s1ab_date_list_file(
                    mintpy_dir=os.getcwd(),
                    safe_list_file=safe_list_file,
                    print_msg=True)

    elif processor == 'aria':
        ## compose input arguments
        # use the default template file if exists & input
        # ARIA 的 prep 阶段会直接生成部分 MintPy 输入文件，因此后面 load_data() 会跳过重复写入。
        default_temp_files = [fname for fname in iDict['template_file']
                              if fname.endswith('smallbaselineApp.cfg')]
        if len(default_temp_files) > 0:
            temp_file = default_temp_files[0]
        else:
            temp_file = iDict['template_file'][0]
        iargs = ['--template', temp_file]

        # file name/dir/path
        # ARG2OPT_DICT 把 prep_aria.py 的命令行参数名映射到 MintPy 模板配置项。
        ARG2OPT_DICT = {
            '--stack-dir'           : 'mintpy.load.unwFile',
            '--unwrap-stack-name'   : 'mintpy.load.unwFile',
            '--coherence-stack-name': 'mintpy.load.corFile',
            '--conn-comp-stack-name': 'mintpy.load.connCompFile',
            '--dem'                 : 'mintpy.load.demFile',
            '--incidence-angle'     : 'mintpy.load.incAngleFile',
            '--azimuth-angle'       : 'mintpy.load.azAngleFile',
            '--water-mask'          : 'mintpy.load.waterMaskFile',
            '--iono'                : 'mintpy.load.ionUnwFile',
        }

        for arg_name, opt_name in ARG2OPT_DICT.items():
            arg_value = iDict.get(opt_name, 'auto')
            if arg_value.lower() not in ['auto', 'no', 'none']:
                # 对目录、文件名、普通文件路径三类参数分别处理。
                if arg_name.endswith('dir'):
                    iargs += [arg_name, os.path.dirname(arg_value)]
                elif arg_name.endswith('name'):
                    iargs += [arg_name, os.path.basename(arg_value)]
                else:
                    iargs += [arg_name, arg_value]

        # configurations
        if iDict['compression']:
            iargs += ['--compression', iDict['compression']]
        if iDict['updateMode']:
            iargs += ['--update']

        ## run
        ut.print_command_line(script_name, iargs)
        prep_module.main(iargs)

    elif processor == 'gmtsar':
        # use the custom template file if exists & input
        # GMTSAR 的 prep 脚本依赖自定义模板；只有默认模板时信息不够。
        custom_temp_files = [fname for fname in iDict['template_file']
                             if not fname.endswith('smallbaselineApp.cfg')]
        if len(custom_temp_files) == 0:
            raise FileExistsError('Custom template file NOT found and is required for GMTSAR!')

        # run prep_*.py
        iargs = [custom_temp_files[0]]
        ut.print_command_line(script_name, iargs)
        try:
            prep_module.main(iargs)
        except:
            warnings.warn('prep_gmtsar.py failed. Assuming its result exists and continue...')

    return


def get_extra_metadata(iDict):
    """Extra metadata with key names in MACRO_CASE to be written into stack file.

    Parameters: iDict     - dict, input arguments from command lines & template file
                extraDict - dict, extra metadata from template file:
                            E.g. PROJECT_NAME, PLATFORM, ORBIT_DIRECTION, SUBSET_X/YMIN, etc.
    """
    extraDict = {}
    # all keys in MACRO_CASE
    # 约定：全大写 key 通常是要写入输出 HDF5 的额外元数据。
    upper_keys = [i for i in iDict.keys() if i.isupper()]
    for key in upper_keys:
        value  = iDict[key]
        if key in ['PROJECT_NAME', 'PLATFORM']:
            if value:
                extraDict[key] = value
        else:
            extraDict[key] = value
    return extraDict


#################################################################
def load_data(inps):
    """load data into HDF5 files."""

    ## 0. read input
    start_time = time.time()
    # 第一步：把命令行参数和模板文件合并成一个总字典 iDict。
    iDict = read_inps2dict(inps)

    ## 1. prepare metadata
    # 第二步：调用 prep_xxx.py 为不同处理器产品准备/补全元数据。
    prepare_metadata(iDict)
    extraDict = get_extra_metadata(iDict)

    # skip data writing as it is included in prep_aria/nisar
    if iDict['processor'] in ['aria', 'nisar']:
        # ARIA/NISAR 的准备脚本已经生成 MintPy 可用文件，因此这里直接返回。
        return

    ## 2. search & write data files
    print('-'*50)
    print('updateMode : {}'.format(iDict['updateMode']))
    print('compression: {}'.format(iDict['compression']))
    print('multilook x/ystep: {}/{}'.format(iDict['xstep'], iDict['ystep']))
    print('multilook method : {}'.format(iDict['method']))
    kwargs = dict(updateMode=iDict['updateMode'], xstep=iDict['xstep'], ystep=iDict['ystep'])

    # read subset info [need the metadata from above]
    # 第三步：读取裁剪范围，决定后面只加载全图还是子区域。
    iDict = read_subset_box(iDict)

    # geometry in geo / radar coordinates
    # 第四步：准备几何文件。这里把几何数据集和观测数据集的模板映射合并，是为了能从观测文件补充尺寸/坐标信息。
    geom_dset_name2template_key = {
        **GEOM_DSET_NAME2TEMPLATE_KEY,
        **IFG_DSET_NAME2TEMPLATE_KEY,
        **OFF_DSET_NAME2TEMPLATE_KEY,
    }
    geom_geo_obj, geom_radar_obj = read_inps_dict2geometry_dict_object(iDict, geom_dset_name2template_key)
    # 几何输出统一写到 ./inputs 目录。
    geom_geo_file = os.path.abspath('./inputs/geometryGeo.h5')
    geom_radar_file = os.path.abspath('./inputs/geometryRadar.h5')

    # use 'lzf' HDF compression for a significantly smaller geometry file size, w/o much impact on the performance
    compression = 'lzf' if iDict['compression'] == 'default' else iDict['compression']

    if run_or_skip(geom_geo_file, geom_geo_obj, iDict['box4geo'], **kwargs) == 'run':
        # write2hdf5() 是 geometryDict 对象的方法，会把外部几何文件读取并写成标准 HDF5。
        geom_geo_obj.write2hdf5(
            outputFile=geom_geo_file,
            access_mode='w',
            box=iDict['box4geo'],
            xstep=iDict['xstep'],
            ystep=iDict['ystep'],
            compression=compression)

    if run_or_skip(geom_radar_file, geom_radar_obj, iDict['box'], **kwargs) == 'run':
        # radar 几何文件会额外写入 PROJECT_NAME、PLATFORM 等元数据。
        geom_radar_obj.write2hdf5(
            outputFile=geom_radar_file,
            access_mode='w',
            box=iDict['box'],
            xstep=iDict['xstep'],
            ystep=iDict['ystep'],
            compression=compression,
            extra_metadata=extraDict)

    # observations: ifgram, ion or offset
    # loop over obs stacks
    # 第五步：依次加载普通干涉图栈、电离层栈、偏移栈。
    stack_ds_name2tmpl_key_list = [
        IFG_DSET_NAME2TEMPLATE_KEY,
        ION_DSET_NAME2TEMPLATE_KEY,
        OFF_DSET_NAME2TEMPLATE_KEY,
    ]
    stack_files = ['ifgramStack.h5', 'ionStack.h5', 'offsetStack.h5']
    stack_files = [os.path.abspath(os.path.join('./inputs', x)) for x in stack_files]
    compression = None if iDict['compression'] == 'default' else iDict['compression']

    for ds_name2tmpl_opt, stack_file in zip(stack_ds_name2tmpl_key_list, stack_files):
        # initiate dict objects
        # 根据模板路径搜索文件，并构建 ifgramStackDict 对象。
        stack_obj = read_inps_dict2ifgram_stack_dict_object(iDict, ds_name2tmpl_opt)

        # use geom_obj as size reference while loading ionosphere
        geom_obj = None
        if os.path.basename(stack_file).startswith('ion'):
            # 电离层栈的尺寸需要参考几何对象，避免和主观测栈坐标/尺寸不一致。
            geom_obj = geom_geo_obj if iDict['geocoded'] else geom_radar_obj

        # write dict objects to HDF5 files
        if run_or_skip(stack_file, stack_obj, iDict['box'], geom_obj=geom_obj, **kwargs) == 'run':
            # write2hdf5() 会真正读取所有外部干涉图文件，写入 ifgramStack.h5/ionStack.h5/offsetStack.h5。
            stack_obj.write2hdf5(
                outputFile=stack_file,
                access_mode='w',
                box=iDict['box'],
                xstep=iDict['xstep'],
                ystep=iDict['ystep'],
                mli_method=iDict['method'],
                compression=compression,
                extra_metadata=extraDict,
                geom_obj=geom_obj)

    # used time
    m, s = divmod(time.time()-start_time, 60)
    print(f'time used: {m:02.0f} mins {s:02.1f} secs.\n')

    return
