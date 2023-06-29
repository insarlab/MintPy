############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Sara Mirzaee, Zhang Yunjun, Bhuvan Varugu, 2018  #
############################################################


import os
import re

import h5py
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from skimage.transform import resize

import mintpy.cli.diff
from mintpy.objects import timeseries
from mintpy.utils import ptime, readfile, utils as ut, writefile


############################################################################
def get_delay_geo(ztd_file, atr, cos_inc_angle):
    """calc single path tropo delay in line-of-sight direction

    Parameters: ztd_file      - str, path of zenith delay file
                atr           - dict, dictionary of attribute for output file
                cos_inc_angle - 2D np.ndarray in float32, cos(inc_angle)
    Returns:    delay         - 2D np.ndarray in float32, LOS delay
    """

    # get geo_box from ts_file attributes
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    geo_box = ut.coordinate(atr).box_pixel2geo(pixel_box=(0, 0, width, length))

    # geo_box --> pix_box in ztd file
    atr_ztd = readfile.read_attribute(ztd_file)
    pix_box = ut.coordinate(atr_ztd).box_geo2pixel(geo_box)

    # read ztd file
    delay = readfile.read(ztd_file, box=pix_box)[0]

    # interpolate/resample into the same resolution as ts_file
    delay = resize(
        delay, (length, width),
        order=1,
        mode='constant',
        anti_aliasing=True,
        preserve_range=True,
    )

    # project from zenith to line-of-sight
    delay /= cos_inc_angle

    # reverse the sign for consistency between different phase correction steps/methods
    delay *= -1

    return delay


def get_delay_radar(ztd_file, cos_inc_angle, pts_new):
    """calc single path tropo delay in line-of-sight direction

    Parameters: ztd_file      - str, path of zenith delay file
                cos_inc_angle - 2D np.ndarray in (len, wid) in float32, cos(inc_angle)
                pts_new       - 2D np.ndarray in (len*wid, 2) in float32
    Returns:    delay         - 2D np.ndarray in float32, LOS delay
    """
    # read ztd file
    delay_ztd, atr_ztd = readfile.read(ztd_file)
    # flip to be consistent with the reversed lats
    delay_ztd = np.flipud(delay_ztd)

    # pixel coordinates in ztd file
    lats, lons = ut.get_lat_lon(atr_ztd, dimension=1)
    # set lats in ascending order as required by RegularGridInterpolator
    lats = np.flipud(lats)
    pts_ztd = ((lats.flatten(),
                lons.flatten()))

    # resample in pts_new coordinates
    interp_func = RegularGridInterpolator(
        pts_ztd,
        delay_ztd,
        method='nearest',
        bounds_error=False,
        fill_value=0,
    )
    delay = interp_func(pts_new)
    delay = delay.reshape(cos_inc_angle.shape)

    # project from zenith to line-of-sight
    delay /= cos_inc_angle

    # reverse the sign for consistency between different phase correction steps/methods
    delay *= -1

    return delay


############################################################################
def calculate_delay_timeseries(tropo_file, dis_file, geom_file, gacos_dir):
    """calculate delay time-series and write to HDF5 file"""

    ## get list of dates
    atr = readfile.read_attribute(dis_file)
    ftype = atr['FILE_TYPE']
    if ftype == 'timeseries':
        date_list = timeseries(dis_file).get_date_list()

    elif ftype == '.unw':
        date12 = readfile.read_attribute(dis_file)['DATE12']
        date_list = ptime.yyyymmdd(date12.split('-'))

    else:
        raise ValueError(f'un-supported displacement file type: {ftype}')

    # list of dates --> list of ztd files
    ztd_files = []
    flag = np.ones(len(date_list), dtype=np.bool_)
    for i, date_str in enumerate(date_list):
        fnames = [os.path.join(gacos_dir, f'{date_str}{fext}') for fext in ['.ztd', '.ztd.tif']]
        fnames = [f for f in fnames if os.path.exists(f)]
        if len(fnames) > 0:
            ztd_files.append(fnames[0])
        else:
            print(f'WARNING: NO ztd file found for {date_str}! Continue without it.')
            flag[i] = False

    # update date_list to be consistent with ztd_files
    if np.any(flag == 0):
        date_list = np.array(date_list)[flag].tolist()

    ## update_mode
    def get_dataset_size(fname):
        atr = readfile.read_attribute(fname)
        return (atr['LENGTH'], atr['WIDTH'])

    def run_or_skip(ztd_files, tropo_file, geom_file):
        print('update mode: ON')
        print(f'output file: {tropo_file}')
        flag = 'skip'

        # check existence and modification time
        if ut.run_or_skip(out_file=tropo_file, in_file=ztd_files, print_msg=False) == 'run':
            flag = 'run'
            print('1) output file either do NOT exist or is NOT newer than all ZTD files.')

        else:
            print('1) output file exists and is newer than all ZTD files.')

            # check dataset size in space / time
            date_list = [str(re.findall(r'\d{8}', i)[0]) for i in ztd_files]
            if (get_dataset_size(tropo_file) != get_dataset_size(geom_file)
                    or any(i not in timeseries(tropo_file).get_date_list() for i in date_list)):
                flag = 'run'
                print(('2) output file does NOT have the same len/wid as the geometry file {}'
                       ' or does NOT contain all dates').format(geom_file))
            else:
                print('2) output file has the same len/wid as the geometry file and contains all dates')

                # check if output file is fully written
                with h5py.File(tropo_file, 'r') as f:
                    if np.all(f['timeseries'][-1,:,:] == 0):
                        flag = 'run'
                        print('3) output file is NOT fully written.')
                    else:
                        print('3) output file is fully written.')

        # result
        print(f'run or skip: {flag}')
        return flag

    if run_or_skip(ztd_files, tropo_file, geom_file) == 'skip':
        return


    ## prepare output file

    # metadata
    atr['FILE_TYPE'] = 'timeseries'
    atr['UNIT'] = 'm'

    # remove metadata related with double reference
    # because absolute delay is calculated and saved
    for key in ['REF_DATE', 'REF_X', 'REF_Y', 'REF_LAT', 'REF_LON']:
        if key in atr.keys():
            atr.pop(key)

    # instantiate time-series
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    num_date = len(date_list)
    dates = np.array(date_list, dtype=np.string_)
    ds_name_dict = {
        "date"       : [dates.dtype, (num_date,), dates],
        "timeseries" : [np.float32,  (num_date, length, width), None],
    }
    writefile.layout_hdf5(tropo_file, ds_name_dict, metadata=atr)


    ## calculate phase delay

    # read geometry
    print(f'read incidenceAngle from file: {geom_file}')
    inc_angle = readfile.read(geom_file, datasetName='incidenceAngle')[0]
    cos_inc_angle = np.cos(inc_angle * np.pi / 180.0)

    if 'Y_FIRST' in atr.keys():
        # No need for data in geo-coordinates
        pts_new = None

    else:
        # Get pixel lat/lon for data in radar-coordinates
        print('get pixel coordinates in geometry file')
        lats, lons = ut.get_lat_lon(atr, geom_file)
        pts_new = np.hstack((lats.reshape(-1, 1),
                             lons.reshape(-1, 1)))

    # loop for date-by-date IO
    prog_bar = ptime.progressBar(maxValue=num_date)
    for i in range(num_date):
        date_str = date_list[i]
        ztd_file = ztd_files[i]

        # calc delay
        if 'Y_FIRST' in atr.keys():
            delay = get_delay_geo(ztd_file, atr, cos_inc_angle)

        else:
            delay = get_delay_radar(ztd_file, cos_inc_angle, pts_new)

        # write delay to file
        block = [i, i+1, 0, length, 0, width]
        writefile.write_hdf5_block(
            tropo_file,
            data=delay,
            datasetName='timeseries',
            block=block,
            print_msg=False,
        )

        prog_bar.update(i + 1, suffix=os.path.basename(ztd_file))
    prog_bar.close()

    return tropo_file


############################################################################
def run_tropo_gacos(inps):

    # calculate tropo delay and savee to h5 file
    calculate_delay_timeseries(
        tropo_file=inps.tropo_file,
        dis_file=inps.dis_file,
        geom_file=inps.geom_file,
        gacos_dir=inps.gacos_dir)

    # correct tropo delay (using diff.py)
    # diff.py can handle different reference in space and time
    # e.g. the absolute delay and the double referenced time-series
    print('correcting delay for using diff.py')
    iargs = [inps.dis_file, inps.tropo_file, '-o', inps.cor_dis_file, '--force']
    print('diff.py', ' '.join(iargs))
    mintpy.cli.diff.main(iargs)

    return
