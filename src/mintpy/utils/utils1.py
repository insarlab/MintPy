"""Miscellaneous utilities - dependent on utils0."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################
# Recommend import:
#   from mintpy.utils import utils as ut


import glob
import os
import re
import shutil
import time

import h5py
import numpy as np

import mintpy
from mintpy.objects import GEOMETRY_DSET_NAMES, deramp, ifgramStack, timeseries
from mintpy.utils import ptime, readfile, writefile
from mintpy.utils.utils0 import *


#################################### Geometry #########################################
def get_center_lat_lon(geom_file, box=None):
    """Get the lat/lon of the scene center"""
    meta = readfile.read_attribute(geom_file)
    if box is None:
        box = (0, 0, int(meta['WIDTH']), int(meta['LENGTH']))

    col_c = int((box[0] + box[2]) / 2)
    row_c = int((box[1] + box[3]) / 2)
    if 'Y_FIRST' in meta.keys():
        lat0 = float(meta['Y_FIRST'])
        lon0 = float(meta['X_FIRST'])
        lat_step = float(meta['Y_STEP'])
        lon_step = float(meta['X_STEP'])
        lat_c = lat0 + lat_step * row_c
        lon_c = lon0 + lon_step * col_c
    else:
        box_c = (col_c, row_c, col_c+1, row_c+1)
        lat_c = float(readfile.read(geom_file, datasetName='latitude', box=box_c)[0])
        lon_c = float(readfile.read(geom_file, datasetName='longitude', box=box_c)[0])
    return lat_c, lon_c


#################################### Data Operation ###################################
def get_residual_std(timeseries_resid_file, mask_file='maskTempCoh.h5', ramp_type='quadratic'):
    """Calculate deramped standard deviation in space for each epoch of input timeseries file.
    Parameters: timeseries_resid_file - string, timeseries HDF5 file,
                    e.g. timeseries_ERA5_demErrInvResid.h5
                mask_file - string, mask file, e.g. maskTempCoh.h5
                ramp_type - string, ramp type, e.g. linear, quadratic, no for do not remove ramp
    Returns:    std_list  - list of float, standard deviation of deramped input timeseries file
                date_list - list of string in YYYYMMDD format, corresponding dates
                std_file  - string, text file with std and date info.
    Example:    import mintpy.utils.utils as ut
                std_list, date_list = ut.get_residual_std('timeseries_ERA5_demErrInvResid.h5',
                                                          'maskTempCoh.h5')[:2]
    """
    # Intermediate files name
    # ramp_type can sometimes be False, thus, should be treated the same as "no"
    if not ramp_type or ramp_type == 'no':
        print('No ramp removal')
        deramped_file = timeseries_resid_file
    else:
        deramped_file = f'{os.path.splitext(timeseries_resid_file)[0]}_ramp.h5'
    std_file = os.path.splitext(deramped_file)[0]+'_std.txt'

    # Get residual std text file
    if run_or_skip(out_file=std_file, in_file=[timeseries_resid_file, mask_file], readable=False) == 'run':
        if run_or_skip(out_file=deramped_file, in_file=timeseries_resid_file) == 'run':
            if not os.path.isfile(timeseries_resid_file):
                msg = 'Can not find input timeseries residual file: '+timeseries_resid_file
                msg += '\nRe-run dem_error.py to generate it.'
                raise Exception(msg)
            else:
                #print('removing a {} ramp from file: {}'.format(ramp_type, timeseries_resid_file))
                deramped_file = run_deramp(
                    timeseries_resid_file,
                    ramp_type=ramp_type,
                    mask_file=mask_file,
                    out_file=deramped_file,
                )
        print('calculating residual standard deviation for each epoch from file: '+deramped_file)
        std_file = timeseries(deramped_file).timeseries_std(maskFile=mask_file, outFile=std_file)

    # Read residual std text file
    print('read timeseries RMS from file: '+std_file)
    fc = np.loadtxt(std_file, dtype=bytes).astype(str)
    std_list = fc[:, 1].astype(np.float32).tolist()
    date_list = list(fc[:, 0])
    return std_list, date_list, std_file


def get_residual_rms(timeseries_resid_file, mask_file='maskTempCoh.h5', ramp_type='quadratic'):
    """Calculate deramped Root Mean Square in space for each epoch of input timeseries file.
    Parameters: timeseries_resid_file : string,
                    timeseries HDF5 file, e.g. timeseries_ERA5_demErrInvResid.h5
                mask_file : string,
                    mask file, e.g. maskTempCoh.h5
                ramp_type : string,
                    ramp type, e.g. linear, quadratic, no for do not remove ramp
    Returns:    rms_list : list of float,
                    Root Mean Square of deramped input timeseries file
                date_list : list of string in YYYYMMDD format,
                    corresponding dates
                rms_file : string, text file with rms and date info.
    Example:
        import mintpy.utils.utils as ut
        rms_list, date_list = ut.get_residual_rms('timeseriesResidual.h5', 'maskTempCoh.h5')
    """
    # Intermediate files name
    # ramp_type can sometimes be False, thus, should be treated the same as "no"
    if not ramp_type or ramp_type == 'no':
        print('No ramp removal')
        deramped_file = timeseries_resid_file
    else:
        deramped_file = f'{os.path.splitext(timeseries_resid_file)[0]}_ramp.h5'
    fdir = os.path.dirname(os.path.abspath(deramped_file))
    fbase = os.path.splitext(os.path.basename(deramped_file))[0]
    rms_file = os.path.join(fdir, f'rms_{fbase}.txt')

    # Get residual RMS text file
    if run_or_skip(out_file=rms_file, in_file=[timeseries_resid_file, mask_file], readable=False) == 'run':
        if run_or_skip(out_file=deramped_file, in_file=timeseries_resid_file) == 'run':
            if not os.path.isfile(timeseries_resid_file):
                msg = 'Can not find input timeseries residual file: '+timeseries_resid_file
                msg += '\nRe-run dem_error.py to generate it.'
                raise Exception(msg)
            else:
                #print('remove {} ramp from file: {}'.format(ramp_type, timeseries_resid_file))
                deramped_file = run_deramp(
                    timeseries_resid_file,
                    ramp_type=ramp_type,
                    mask_file=mask_file,
                    out_file=deramped_file,
                )

        print('\ncalculating residual RMS for each epoch from file: '+deramped_file)
        rms_file = timeseries(deramped_file).timeseries_rms(
            maskFile=mask_file,
            outFile=rms_file,
        )

    # Read residual RMS text file
    print('read timeseries residual RMS from file: '+rms_file)
    fc = np.loadtxt(rms_file, dtype=bytes).astype(str)
    rms_list = fc[:, 1].astype(np.float32).tolist()
    date_list = list(fc[:, 0])
    return rms_list, date_list, rms_file


def nonzero_mask(File, out_file='maskConnComp.h5', datasetName=None):
    """Generate mask file for non-zero value of input multi-group hdf5 file"""
    atr = readfile.read_attribute(File)
    k = atr['FILE_TYPE']
    if k == 'ifgramStack':
        mask = ifgramStack(File).nonzero_mask(datasetName=datasetName)
    else:
        print('Only ifgramStack file is supported for now, input is '+k)
        return None

    atr['FILE_TYPE'] = 'mask'
    writefile.write(mask, out_file=out_file, metadata=atr)
    return out_file


def spatial_average(fname, datasetName='coherence', maskFile=None, box=None,
                    saveList=False, checkAoi=True, reverseMask=False, threshold=None):
    """Read/Calculate Spatial Average of input file.

    If input file is text file, read it directly;
    If input file is data matrix file:
        If corresponding text file exists with the same mask file/AOI info, read it directly;
        Otherwise, calculate it from data file.

        Only non-nan pixel is considered.
    Parameters: fname       - string, path of input file
                maskFile    - string, path of mask file, e.g. maskTempCoh.h5
                box         - 4-tuple defining the left, upper, right, and lower pixel coordinate
                saveList    - bool, save (list of) mean value into text file
                reverseMask - bool, perform analysis within masked regions instead of outside of them
                threshold   - float, calculate area ratio above threshold instead of spatial average
    Returns:    meanList    - list(float) or float, average value in space for each epoch of input file
                dateList    - list(str) or str, for date info
                              date12_list, e.g. 101120-110220, for interferograms/coherence
                              date8_list, e.g. 20101120, for timeseries
                              file name, e.g. velocity.h5, for all the other file types
    Example:    meanList = spatial_average('inputs/ifgramStack.h5')[0]
                meanList, date12_list = spatial_average('inputs/ifgramStack.h5',
                                                        maskFile='maskTempCoh.h5',
                                                        saveList=True)
    """
    def read_text_file(fname):
        txtContent = np.loadtxt(fname, dtype=bytes).astype(str)
        meanList = [float(i) for i in txtContent[:, 1]]
        dateList = [i for i in txtContent[:, 0]]
        return meanList, dateList

    # Baic File Info
    atr = readfile.read_attribute(fname)
    k = atr['FILE_TYPE']
    if not box:
        box = (0, 0, int(atr['WIDTH']), int(atr['LENGTH']))

    # default output filename
    prefix = datasetName if k == 'ifgramStack' else os.path.splitext(os.path.basename(fname))[0]
    suffix = 'SpatialAvg' if threshold is None else 'AreaRatio'
    suffix += 'RevMsk' if reverseMask else ''
    txtFile = prefix + suffix + '.txt'

    # If input is text file
    if fname.endswith(suffix):
        print('Input file is spatial average txt already, read it directly')
        meanList, dateList = read_text_file(fname)
        return meanList, dateList

    # Read existing txt file only if 1) data file is older AND 2) same AOI
    file_line  = f'# Data file: {os.path.basename(fname)}\n'
    mask_line  = f'# Mask file: {maskFile}\n'
    aoi_line   = f'# AOI box: {box}\n'
    thres_line = f'# Threshold: {threshold}\n'

    try:
        # Read AOI line from existing txt file
        fl = open(txtFile)
        lines = fl.readlines()
        fl.close()
        # 1. aoi
        if checkAoi:
            try:
                aoi_line_orig = [i for i in lines if '# AOI box:' in i][0]
            except:
                aoi_line_orig = ''
        else:
            aoi_line_orig = aoi_line
        # 2. mask file
        try:
            mask_line_orig = [i for i in lines if '# Mask file:' in i][0]
        except:
            mask_line_orig = ''
        # 3. mask file - modification time
        update_mask_file = run_or_skip(out_file=txtFile, in_file=[maskFile], readable=False)
        # 4. data file - modification time
        if k == 'ifgramStack':
            with h5py.File(fname, 'r') as f:
                ti = float(f[datasetName].attrs.get('MODIFICATION_TIME', os.path.getmtime(fname)))
        else:
            ti = os.path.getmtime(fname)
        to = os.path.getmtime(txtFile)

        if (aoi_line_orig == aoi_line
                and mask_line_orig == mask_line
                and update_mask_file == 'skip'
                and ti <= to):
            print(txtFile+' already exists, read it directly')
            meanList, dateList = read_text_file(txtFile)
            return meanList, dateList
    except:
        pass

    # use median instead of mean for offset measurement
    if datasetName and 'offset' in datasetName:
        useMedian = True
    else:
        useMedian = False

    # Calculate mean coherence or area ratio list
    if k == 'ifgramStack':
        obj = ifgramStack(fname)
        obj.open(print_msg=False)
        meanList, dateList = obj.spatial_average(
            datasetName=datasetName,
            maskFile=maskFile,
            box=box,
            useMedian=useMedian,
            reverseMask=reverseMask,
            threshold=threshold,
        )
        pbase = obj.pbaseIfgram
        tbase = obj.tbaseIfgram
        obj.close()

    elif k == 'timeseries':
        meanList, dateList = timeseries(fname).spatial_average(
            maskFile=maskFile,
            box=box,
            reverseMask=reverseMask,
            threshold=threshold,
        )

    else:
        data = readfile.read(fname, box=box)[0]
        if maskFile and os.path.isfile(maskFile):
            print('mask from file: '+maskFile)
            mask = readfile.read(maskFile, datasetName='mask', box=box)[0]
            data[mask == int(reverseMask)] = np.nan

        # calculate area ratio if threshold is specified
        # percentage of pixels with value above the threshold
        if threshold is not None:
            data[data > threshold] = 1
            data[data <= threshold] = 0

        meanList = np.nanmean(data)
        dateList = [os.path.basename(fname)]

    # Write mean coherence list into text file
    if saveList:
        print('write average value in space into text file: '+txtFile)
        fl = open(txtFile, 'w')
        # Write comments
        fl.write(file_line+mask_line+aoi_line+thres_line)
        # Write data list
        numLine = len(dateList)
        if k == 'ifgramStack':
            fl.write('#\tDATE12\t\tMean\tBtemp/days\tBperp/m\t\tNum\n')
            for i in range(numLine):
                fl.write('%s\t%.4f\t%8.0f\t%8.1f\t%d\n' %
                         (dateList[i], meanList[i], tbase[i], pbase[i], i))
        else:
            fl.write('#\tDATE12\t\tMean\n')
            for i in range(numLine):
                fl.write(f'{dateList[i]}\t{meanList[i]:.4f}\n')
        fl.close()

    # read from text file (in 1e-4 precision)
    # to ensure output value consistency
    meanList, dateList = read_text_file(txtFile)

    if len(meanList) == 1:
        meanList = meanList[0]
        dateList = dateList[0]

    return meanList, dateList


def temporal_average(fname, datasetName='coherence', updateMode=False, outFile=None):
    """Calculate temporal average of multi-temporal dataset, equivalent to stacking
    For ifgramStack/unwrapPhase, return average phase velocity

    Parameters: fname       - str, file to be averaged in time
                datasetName - str, dataset to be read from input file, for multiple
                              datasets file - ifgramStack - only
                              e.g.: coherence, unwrapPhase
                updateMode  - bool
                outFile     - str, output filename
                              None for auto output filename
                              False for do not save as output file
    Returns:    dataMean    - 2D np.ndarray
                outFile     - str, output file name
    Examples:   avgPhaseVel = ut.temporal_average('ifgramStack.h5', datasetName='unwrapPhase')[0]
                ut.temporal_average('ifgramStack.h5', datasetName='coherence',
                                    outFile='avgSpatialCoh.h5', updateMode=True)
    """
    atr = readfile.read_attribute(fname, datasetName=datasetName)
    k = atr['FILE_TYPE']
    if k not in ['ifgramStack', 'timeseries']:
        print(f'WARNING: input file is not multi-temporal file: {fname}, return itself.')
        data = readfile.read(fname)[0]
        return data, fname

    # Default output filename
    if outFile is None:
        ext = os.path.splitext(fname)[1]
        if not outFile:
            if k == 'ifgramStack':
                if datasetName == 'coherence':
                    outFile = 'avgSpatialCoh.h5'
                elif 'unwrapPhase' in datasetName:
                    outFile = 'avgPhaseVelocity.h5'
                else:
                    outFile = f'avg{datasetName}.h5'
            elif k == 'timeseries':
                if k in fname:
                    processMark = os.path.basename(fname).split('timeseries')[1].split(ext)[0]
                    outFile = f'avgDisplacement{processMark}.h5'
            else:
                outFile = f'avg{fname}.h5'

    if updateMode and os.path.isfile(outFile):
        dataMean = readfile.read(outFile)[0]
        return dataMean, outFile

    # Calculate temporal average
    if k == 'ifgramStack':
        dataMean = ifgramStack(fname).temporal_average(datasetName=datasetName)
        if 'unwrapPhase' in datasetName:
            atr['FILE_TYPE'] = 'velocity'
            atr['UNIT'] = 'm/year'
        else:
            atr['FILE_TYPE'] = datasetName
    elif k == 'timeseries':
        dataMean = timeseries(fname).temporal_average()
        atr['FILE_TYPE'] = 'displacement'

    if outFile:
        writefile.write(dataMean, out_file=outFile, metadata=atr)
    return dataMean, outFile



#################################### File IO ##########################################
def get_file_list(file_list, abspath=False, coord=None):
    """Get all existed files matching the input list of file pattern
    Parameters: file_list - string or list of string, input file/directory pattern
                abspath - bool, return absolute path or not
                coord - string, return files with specific coordinate type: geo or radar
                    if none, skip the checking and return all files
    Returns:    file_list_out - list of string, existed file path/name, [] if not existed
    Example:    file_list = get_file_list(['*velocity*.h5','timeseries*.h5'])
                file_list = get_file_list('timeseries*.h5')
    """
    if not file_list:
        return []

    if isinstance(file_list, str):
        file_list = [file_list]

    # Get rid of None element
    file_list = [x for x in file_list if x is not None]
    file_list_out = []
    for fname in file_list:
        fnames = glob.glob(fname)
        file_list_out += sorted(list(set(fnames) - set(file_list_out)))

    if abspath:
        file_list_out = [os.path.abspath(i) for i in file_list_out]

    if coord is not None:
        for fname in list(file_list_out):
            atr = readfile.read_attribute(fname)
            if coord in ['geo']:
                if 'Y_FIRST' not in atr.keys():
                    file_list_out.remove(fname)
            elif coord in ['radar', 'rdr', 'rdc']:
                if 'Y_FIRST' in atr.keys():
                    file_list_out.remove(fname)
            else:
                msg = f'un-recognized input coord type: {coord}'
                raise ValueError(msg)
    return file_list_out


def get_lookup_file(filePattern=None, abspath=False, print_msg=True):
    """Find lookup table file with/without input file pattern
    Parameters: filePattern - list of str
                abspath     - bool, return absolute path or not
                print_msg   - bool, printout message or not
    Returns:    outFile     - str, path of the lookup file
    """
    # Search Existing Files
    if not filePattern:
        fileList = ['geometryRadar.h5',
                    'geometryGeo_tight.h5', 'geometryGeo.h5',
                    'geomap*lks_tight.trans', 'geomap*lks.trans',
                    'sim*_tight.UTM_TO_RDC', 'sim*.UTM_TO_RDC']
        dirList = ['inputs', '', '../inputs']

        # file/dirList --> filePattern
        filePattern = []
        for dirname in dirList:
            filePattern += [os.path.join(dirname, fname) for fname in fileList]

    existFiles = []
    try:
        existFiles = get_file_list(filePattern)
    except:
        if print_msg:
            print('ERROR: No geometry / lookup table file found!')
            print('It should be like:')
            print(filePattern)
        return None

    # Check Files Info
    outFile = None
    for fname in existFiles:
        readfile.read_attribute(fname)
        for dsName in ['longitude', 'rangeCoord']:
            try:
                readfile.read(fname, datasetName=dsName, print_msg=False)
                outFile = fname
                break
            except:
                pass

    if not outFile:
        if print_msg:
            print('No lookup table (longitude or rangeCoord) found in files.')
        return None

    # Path Format
    if abspath:
        outFile = os.path.abspath(outFile)
    return outFile


def get_geometry_file(dset_list, work_dir=None, coord='geo', abspath=True, print_msg=True):
    """Find geometry file containing input specific dataset"""
    if isinstance(dset_list, str):
        dset_list = [dset_list]
    for dset in dset_list:
        if dset not in GEOMETRY_DSET_NAMES:
            raise ValueError(f'unrecognized geometry dataset name: {dset}')

    if not work_dir:
        work_dir = os.getcwd()

    # search *geometry*.h5 files
    fname_list = [os.path.join(work_dir, i) for i in ['*geometry*.h5', '*/*geometry*.h5', '../*/geometry*.h5']]
    fname_list = get_file_list(fname_list, coord=coord)
    if len(fname_list) == 0:
        if print_msg:
            print('No geometry file found.')
        return None

    # check dset in the existing h5 files
    for fname in list(fname_list):   #use list() as temp copy to handle modifying list during the loop
        if any(dset not in readfile.get_dataset_list(fname) for dset in dset_list):
            fname_list.remove(fname)
    if len(fname_list) == 0:
        if print_msg:
            print(f'No geometry file with dataset {dset_list} found')
        return None

    geom_file = fname_list[0]
    if abspath:
        geom_file = os.path.abspath(geom_file)
    return geom_file


def update_template_file(template_file, extra_dict, delimiter='='):
    """Update option value in template_file with value from input extra_dict"""
    # Compare and skip updating template_file if no new option value found.
    update = False
    orig_dict = readfile.read_template(template_file)
    for key, value in orig_dict.items():
        if key in extra_dict.keys() and extra_dict[key] != value:
            update = True
    if not update:
        print('No new option value found, skip updating '+template_file)
        return template_file

    # Update template_file with new value from extra_dict
    tmp_file = template_file+'.tmp'
    f_tmp = open(tmp_file, 'w')
    for line in open(template_file):
        c = [i.strip() for i in line.strip().split(delimiter, 1)]
        if not line.startswith(('%', '#')) and len(c) > 1:
            key = c[0]
            value = str.replace(c[1], '\n', '').split("#")[0].strip()

            if key in extra_dict.keys() and extra_dict[key] != value:
                # prepare value string to search & replace following "re" expression syntax
                # link: https://docs.python.org/3/library/re.html
                value2search = value
                # 1. interpret special symbols as characters
                for symbol in ['*', '[', ']', '(', ')']:
                    value2search = value2search.replace(symbol, fr"\{symbol}")
                # 2. use "= {OLD_VALUE}" for search/replace to be more robust
                # against the scenario when key name contains {OLD_VALUE}
                # i.e. mintpy.load.autoPath
                value2search = delimiter+r'[\s]*'+value2search

                old_value_str = re.findall(value2search, line)[0]
                new_value_str = old_value_str.replace(value, extra_dict[key])
                line = line.replace(old_value_str, new_value_str, 1)
                print(f'    {key}: {value} --> {extra_dict[key]}')

        f_tmp.write(line)
    f_tmp.close()

    # Overwrite existing original template file
    shutil.move(tmp_file, template_file)
    return template_file


def add_attribute(fname, atr_new=dict(), print_msg=False):
    """Add/update input attribute of the give file.

    Parameters: fname   - string, path/name of file
                atr_new - dict, attributes to be added/updated
                          if value is None, delete the item from input file attributes
    Returns:    fname   - string, path/name of updated file
    """
    vprint = print if print_msg else lambda *args, **kwargs: None

    # read existing attributes
    atr = readfile.read_attribute(fname)
    key_list = list(atr.keys())

    # compare new attributes with existing ones
    update = update_attribute_or_not(atr_new, atr)
    if not update:
        vprint('All updated (removed) attributes already exists (do not exists)'
               ' and have the same value, skip update.')
        return fname

    # update attributes in the input data file
    fext = os.path.splitext(fname)[1]
    if fext in ['.h5', '.he5']:
        with h5py.File(fname, 'r+') as f:
            for key, value in iter(atr_new.items()):

                if value == 'None' or value is None:
                    # delete the item for invalid input (None)
                    if key in key_list:
                        f.attrs.pop(key)
                        vprint(f'remove {key}')
                else:
                    # update the item for valid input
                    f.attrs[key] = str(value)
                    vprint(f'add/update {key} = {str(value)}')

    else:
        for key, value in iter(atr_new.items()):

            if value == 'None' or value is None:
                # delete the item for invalid input (None)
                if key in key_list:
                    atr.pop(key)
                    vprint(f'remove {key}')
            else:
                # update the item for valid input
                atr[key] = str(value)
                vprint(f'add/update {key} = {str(value)}')

        # write to RSC file
        writefile.write_roipac_rsc(atr, fname+'.rsc', print_msg=print_msg)

    return fname


def check_file_size(fname_list, mode_width=None, mode_length=None):
    """Check file size in the list of files, and drop those not in the same size with majority."""
    # If input file list is empty
    if not fname_list:
        return fname_list, None, None

    # Read Width/Length list
    width_list = []
    length_list = []
    for fname in fname_list:
        atr = readfile.read_attribute(fname)
        width_list.append(atr['WIDTH'])
        length_list.append(atr['LENGTH'])

    # Mode of Width and Length
    mode_width = mode_width if mode_width else most_common(width_list)
    mode_length = mode_length if mode_length else most_common(length_list)

    # Update Input List
    fname_list_out = list(fname_list)
    if (width_list.count(mode_width) != len(width_list)
            or length_list.count(mode_length) != len(length_list)):
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('WARNING: Some files may have the wrong dimensions!')
        print('All files should have the same size.')
        print('The width and length of the majority of files are: %s, %s' %
              (mode_width, mode_length))
        print('But the following files have different dimensions and thus will not be loaded:')
        for fname, length, width in zip(fname_list, length_list, width_list):
            if width != mode_width or length != mode_length:
                print(f'{fname}    width: {width}  length: {length}')
                fname_list_out.remove(fname)
        print('\nNumber of files left: '+str(len(fname_list_out)))
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    return fname_list_out, mode_width, mode_length



#################################### Interaction ##########################################
def is_file_exist(file_list, abspath=True):
    """Check if any file in the file list 1) exists and 2) readable
    Parameters: file_list : str or list(str), file name with/without wildcards
                abspath   : bool, return absolute file name/path or not
    Returns:    file_path : string, found file name/path; None if not.
    """
    try:
        file = get_file_list(file_list, abspath=abspath)[0]
        readfile.read_attribute(file)
    except:
        file = None
    return file


def run_or_skip(out_file, in_file=None, readable=True, print_msg=True):
    """Check whether to update out_file or not.
    return run if any of the following meets:
        1. out_file is empty, e.g. None, []
        2. out_file is not existed
        3. out_file is not readable by readfile.read_attribute() when readable=True
        4. out_file is older than in_file, if in_file is not None
    Otherwise, return skip.

    If in_file=None and out_file exists and readable, return skip

    Parameters: out_file  - string or list of string, output file(s)
                in_file   - string or list of string, input file(s)
                readable  - bool, check if the 1st output file has attribute 'WIDTH'
                print_msg - bool, print message
    Returns:    run/skip  - str, whether to update output file or not
    Example:    if ut.run_or_skip(out_file='timeseries_ERA5_demErr.h5', in_file='timeseries_ERA5.h5'):
                if ut.run_or_skip(out_file='exclude_date.txt',
                                  in_file=['timeseries_ERA5_demErrInvResid.h5',
                                           'maskTempCoh.h5',
                                           'smallbaselineApp.cfg'],
                                  readable=False):
    """
    # 1 - check existence of output files
    if not out_file:
        return 'run'
    else:
        if isinstance(out_file, str):
            out_file = [out_file]
        if not all(os.path.isfile(i) for i in out_file):
            return 'run'

    # 2 - check readability of output files
    if readable:
        try:
            readfile.read_attribute(out_file[0])['WIDTH']
        except:
            if print_msg:
                print(f'{out_file[0]} exists, but can not read, remove it.')
            os.remove(out_file[0])
            return 'run'

    # 3 - check modification time of output and input files
    if in_file:
        in_file = get_file_list(in_file)
        # Check modification time
        if in_file:
            t_in  = max(os.path.getmtime(i) for i in in_file)
            t_out = min(os.path.getmtime(i) for i in out_file)
            if t_in > t_out:
                return 'run'
            elif print_msg:
                print(f'{out_file} exists and is newer than {in_file} --> skip.')
    return 'skip'


def check_template_auto_value(templateDict, auto_file='defaults/smallbaselineApp_auto.cfg'):
    """Replace auto value based on the input auto config file."""
    ## Read default template value and turn yes/no to True/False
    templateAutoFile = os.path.join(os.path.dirname(mintpy.__file__), auto_file)
    templateAutoDict = readfile.read_template(templateAutoFile)

    # if cluster != local, change auto value of numWorker
    cluster_key = 'mintpy.compute.cluster'
    cluster = templateDict.get(cluster_key, 'auto').lower()
    if cluster == 'auto':
        cluster = templateAutoDict[cluster_key]

    if cluster != 'local':
        templateAutoDict['mintpy.compute.numWorker'] = '40'

    ## Update auto value of input template dict
    for key, value in templateDict.items():
        if value == 'auto' and key in templateAutoDict.keys():
            templateDict[key] = templateAutoDict[key]

    # Change yes --> True, no --> False and none --> None
    special_values = {
        'yes'  : True,
        'true' : True,
        'no'   : False,
        'false': False,
        'none' : None,
    }
    for key, value in templateDict.items():
        value = value.lower()
        if value in special_values.keys():
            templateDict[key] = special_values[value]

    return templateDict


def run_deramp(fname, ramp_type, mask_file=None, out_file=None, datasetName=None,
               save_ramp_coeff=False, extra_meta=None):
    """ Remove ramp from each 2D matrix of input file
    Parameters: fname           - str, data file to be deramped
                ramp_type       - str, name of ramp to be estimated.
                mask_file       - str, file of mask of pixels used for ramp estimation
                out_file        - str, output file name
                datasetName     - str, output dataset name, for ifgramStack file type only
                save_ramp_coeff - bool, save the estimated ramp coefficients to text file
                extra_meta      - dict, extra metadata to add to the output file
    Returns:    out_file        - str, output file name
    """

    start_time = time.time()

    # file/dir
    fdir = os.path.dirname(fname)
    fbase, fext = os.path.splitext(os.path.basename(fname))

    # metadata
    atr = readfile.read_attribute(fname)
    ftype = atr['FILE_TYPE']
    length = int(atr['LENGTH'])
    width = int(atr['WIDTH'])

    print(f'remove {ramp_type} ramp from file: {fname}')
    out_file = out_file if out_file else os.path.join(fdir, f'{fbase}_ramp{fext}')
    # ignore out_file for ifgramStack (write back to the same HDF5 file)
    if ftype == 'ifgramStack':
        out_file = fname

    # mask
    if os.path.isfile(mask_file):
        mask = readfile.read(mask_file)[0]
        print('read mask file: '+mask_file)
    else:
        mask = np.ones((length, width), dtype=np.bool_)
        print('use mask of the whole area')

    # write coefficient of specified surface function fit
    coeff_file = None
    if save_ramp_coeff:
        coeff_file = os.path.join(fdir, f'rampCoeff_{fbase}.txt')
        with open(coeff_file, 'w') as f:
            f.write(f'# input  file: {fname}\n')
            f.write(f'# output file: {out_file}\n')
            f.write(f'# ramp type: {ramp_type}\n')

    # deramping
    if ftype == 'timeseries':
        # write HDF5 file with defined metadata and (empty) dataset structure
        writefile.layout_hdf5(out_file, ref_file=fname, print_msg=True)

        print('estimating phase ramp one date at a time ...')
        date_list = timeseries(fname).get_date_list()
        num_date = len(date_list)
        prog_bar = ptime.progressBar(maxValue=num_date)
        for i in range(num_date):
            if coeff_file:
                # prepend epoch name to line of coefficients
                with open(coeff_file, 'a') as f:
                    f.write(f'{(date_list[i])}    ')

            # read
            data = readfile.read(fname, datasetName=date_list[i])[0]

            # deramp
            data = deramp(
                data,
                mask,
                ramp_type=ramp_type,
                metadata=atr,
                coeff_file=coeff_file,
            )[0]

            # write
            writefile.write_hdf5_block(
                out_file,
                data=data,
                datasetName='timeseries',
                block=[i, i+1, 0, length, 0, width],
                print_msg=False,
            )

            prog_bar.update(i+1, suffix=f'{i+1}/{num_date}')
        prog_bar.close()
        print(f'finished writing to file: {out_file}')

    elif ftype == 'ifgramStack':
        obj = ifgramStack(fname)
        obj.open(print_msg=False)
        if not datasetName:
            datasetName = 'unwrapPhase'

        with h5py.File(fname, 'a') as f:
            ds = f[datasetName]
            dsNameOut = f'{datasetName}_ramp'
            if dsNameOut in f.keys():
                dsOut = f[dsNameOut]
                print(f'access HDF5 dataset /{dsNameOut}')
            else:
                dsOut = f.create_dataset(
                    dsNameOut,
                    shape=(obj.numIfgram, length, width),
                    dtype=np.float32,
                    chunks=True,
                    compression=None)
                print(f'create HDF5 dataset /{dsNameOut}')

            prog_bar = ptime.progressBar(maxValue=obj.numIfgram)
            for i in range(obj.numIfgram):
                if coeff_file:
                    # prepend IFG date12 to line of coefficients
                    with open(coeff_file, 'a') as f:
                        f.write(f'{str(obj.date12List[i])}    ')

                # read
                data = ds[i, :, :]

                # deramp
                data = deramp(
                    data,
                    mask,
                    ramp_type=ramp_type,
                    metadata=atr,
                    coeff_file=coeff_file,
                )[0]

                # write
                dsOut[i, :, :] = data

                prog_bar.update(i+1, suffix=f'{i+1}/{obj.numIfgram}')
            prog_bar.close()
            print(f'finished writing to file: {fname}')

    # Single Dataset File
    else:
        if coeff_file:
            # prepend file-type to line of coefficients
            with open(coeff_file, 'a') as f:
                f.write('{}    '.format(atr['FILE_TYPE']))

        # read
        if not datasetName and ftype == 'velocity':
            datasetName = 'velocity'
        data = readfile.read(fname, datasetName=datasetName)[0]

        # deramp
        data = deramp(
            data,
            mask,
            ramp_type=ramp_type,
            metadata=atr,
            coeff_file=coeff_file,
        )[0]

        # write
        print(f'writing >>> {out_file}')
        writefile.write(data, out_file=out_file, ref_file=fname)

    # add extra_meta to the output file
    if extra_meta:
        print('add/update the following metadata to file:')
        add_attribute(out_file, extra_meta, print_msg=True)

    # used time
    m, s = divmod(time.time()-start_time, 60)
    print(f'time used: {m:02.0f} mins {s:02.1f} secs.')

    return out_file
