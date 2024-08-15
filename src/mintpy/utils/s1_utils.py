"""Utilities for Sentinel-1."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Aug 2021                           #
############################################################
# Recommend import:
#   from mintpy.utils import s1_utils


import os
import re

import numpy as np

from mintpy.objects import timeseries
from mintpy.utils import ptime, time_func


def estimate_s1ab_bias(mintpy_dir, ts_dis, safe_list_file=None):
    """Estimate the bias between Sentinel-1 A and B.
    Parameters: mintpy_dir     - str, path of the mintpy working directory
                ts_dis         - 2D np.ndarray in size of (num_date, num_pixel) in float32
                safe_list_file - str, path of the SAFE_files.txt file
    Returns:    bias           - 1D np.ndarray in size of (num_pixel) in float32
                flagA/B        - 1D np.ndarray in size of (num_date) in bool
                dates_fit      - list of datetime.datetime objects
                ts_fitA/B      - 2D np.ndarray in size of (num_date_fit, num_pixel) in float32
    """

    # dates/flags for S1A/B
    (s1a_date_list_file,
     s1b_date_list_file) = get_s1ab_date_list_file(mintpy_dir, safe_list_file, print_msg=False)
    date_listA = np.loadtxt(s1a_date_list_file, dtype=str).tolist()
    date_listB = np.loadtxt(s1b_date_list_file, dtype=str).tolist()
    date_list = sorted(date_listA + date_listB)
    num_date = len(date_list)

    min_date = date_listB[0]
    flagA = np.array([x in date_listA and x >= min_date for x in date_list], dtype=np.bool_)
    flagB = np.array([x in date_listB and x >= min_date for x in date_list], dtype=np.bool_)
    # update date_list to the shared time period only
    date_listA = np.array(date_list)[flagA].tolist()
    date_listB = np.array(date_list)[flagB].tolist()

    if not date_listA or not date_listB:
        sname = 'S1A' if not date_listA else 'S1B'
        msg = f'WARNING: NO {sname} acquisitions in the time series, thus,'
        msg += ' can NOT estimate S1A/B bias from it.'
        print(msg)
        return None, flagA, flagB, None, None, None

    # fit
    ts_dis = ts_dis.reshape(num_date, -1)
    model = dict(polynomial=1)
    mA = time_func.estimate_time_func(model, date_listA, ts_dis[flagA, :], ref_date=date_listA[0])[1]
    mB = time_func.estimate_time_func(model, date_listB, ts_dis[flagB, :], ref_date=date_listB[0])[1]

    # grab bias/offset from the fitting time-series
    date_list_fit = ptime.get_date_range(min_date, date_list[-1], dstep=1)
    dates_fit = ptime.date_list2vector(date_list_fit)[0]
    GA_fit = time_func.get_design_matrix4time_func(date_list_fit, model, ref_date=date_listA[0])
    GB_fit = time_func.get_design_matrix4time_func(date_list_fit, model, ref_date=date_listB[0])
    ts_fitA = np.matmul(GA_fit, mA)
    ts_fitB = np.matmul(GB_fit, mB)
    bias = np.median(ts_fitB - ts_fitA, axis=0)

    # ignore zero bias values
    bias[bias == 0] = np.nan

    return bias, flagA, flagB, dates_fit, ts_fitA, ts_fitB


def get_s1ab_date_list_file(mintpy_dir, safe_list_file=None, print_msg=True):
    """Get (and generate if not exist) the date list file of S1A/B.
    Parameters: mintpy_dir           - str, path of mintpy working directory
                safe_list_file       - str, path of SAFE_files.txt file
    Returns:    s1a/b_date_list_file - str, path of S1A/B_date.txt file
    """
    vprint = print if print_msg else lambda *args, **kwargs: None
    mintpy_dir = os.path.abspath(mintpy_dir)
    s1a_date_list_file = os.path.join(mintpy_dir, 'S1A_date.txt')
    s1b_date_list_file = os.path.join(mintpy_dir, 'S1B_date.txt')
    if not os.path.isfile(s1a_date_list_file):
        # get SAFE list filename
        if not safe_list_file:
            safe_list_file = os.path.join(os.path.dirname(mintpy_dir), 'SAFE_files.txt')
        if not os.path.isfile(safe_list_file):
            msg = f'Required file NOT found in: {safe_list_file}!'
            msg += '\nIt can be generated as: "ls ./SLC > SAFE_files.txt".'
            raise FileNotFoundError(msg)

        # get date/sensor_list
        vprint('\nread sensor info from file:', safe_list_file)
        ts_files = [os.path.join(mintpy_dir, f'timeseries{x}.h5') for x in ['', 'Rg', 'Az']]
        ts_file = [x for x in ts_files if os.path.isfile(x)][0]
        date_list = timeseries(ts_file).get_date_list()
        sensor_list = safe_list_file2sensor_list(safe_list_file, date_list, print_msg=False)[0]

        # write to text file for easy access by other scripts
        s1a_date_list = [i for i, j in zip(date_list, sensor_list) if j == 'S1A']
        s1b_date_list = [i for i, j in zip(date_list, sensor_list) if j == 'S1B']
        np.savetxt(s1a_date_list_file, np.array(s1a_date_list).reshape(-1,1), fmt='%s')
        vprint(f'write file: {s1a_date_list_file}')
        if len(s1b_date_list) > 0:
            np.savetxt(s1b_date_list_file, np.array(s1b_date_list).reshape(-1,1), fmt='%s')
            vprint(f'write file: {s1b_date_list_file}')
    else:
        vprint(f'S1A/B_date.txt files exist in: {mintpy_dir}.')

    return s1a_date_list_file, s1b_date_list_file


def safe_list_file2sensor_list(safe_list_file, date_list=None, print_msg=True):
    """Get list of Sentinel-1 sensor names from txt file with SAFE file names.

    Parameters: safe_list_file - str, path of the text file with Sentinel-1 SAFE file path
                                 E.g. SAFE_files.txt
                date_list      - list of str in YYYYMMDD format, reference list of dates
    Returns:    sensor_list    - list of str in S1A or S1B
                date_list      - list of str in YYYYMMDD format
    Example:
        date_list = timeseries('timeseries.h5').get_date_list()
        sensor_list = safe_list_file2sensor_list('../SAFE_files.txt',
                                                 date_list=date_list,
                                                 print_msg=False)[0]
        s1b_dates = [i for i, j in zip(date_list, sensor_list) if j == 'S1B']
        np.savetxt('S1B_date.txt', np.array(s1b_dates).reshape(-1,1), fmt='%s')
    """
    # read txt file
    fc = np.loadtxt(safe_list_file, dtype=str).astype(str).tolist()
    safe_fnames = [os.path.basename(i) for i in fc]

    # get date_list
    date_list_out = [re.findall(r'_\d{8}T', i)[0][1:-1] for i in safe_fnames]
    date_list_out = sorted(list(set(date_list_out)))

    # get sensor_list
    sensor_list = []
    for d in date_list_out:
        safe_fname = [i for i in safe_fnames if d in i][0]
        sensor = safe_fname.split('_')[0]
        sensor_list.append(sensor)

    # update against date_list_out
    if date_list is not None:
        # check possible missing dates
        dates_missing = [i for i in date_list if i not in date_list_out]
        if dates_missing:
            raise ValueError(f'The following dates are missing:\n{dates_missing}')

        # prune dates not-needed
        flag = np.array([i in date_list for i in date_list_out], dtype=np.bool_)
        if np.sum(flag) > 0:
            sensor_list = np.array(sensor_list)[flag].tolist()
            dates_removed = np.array(date_list_out)[~flag].tolist()
            date_list_out = np.array(date_list_out)[flag].tolist()
            if print_msg:
                print(f'The following dates are not needed and removed:\n{dates_removed}')

    return sensor_list, date_list


def get_subswath_masks(flag, cut_overlap_in_half=False):
    """Get the 3 masks for each of the Sentinel-1 subswath.
    Parameters: flag                - 2D np.ndarray in size of (length, width) in bool for valid observations
                cut_overlap_in_half - bool, turn on for offset estimated with very large chip size
    Returns:    mask1/2/3           - 2D np.ndarray in size of (length, width) in bool
                box1/2/3            - list of (x0, y0, x1, y1) in int
    Examples:   flag = readfile.read('inputs/geometryRadar.h5', datasetName='height')[0] != 0
                mask1, mask2, mask3 = s1_utils.get_subswath_masks(flag)[:3]
    """
    length, width = flag.shape
    iw1_x0 = 0
    iw3_x1 = width

    # get ymin/max based on the rough center colomn number for each subswath
    x1, x2, x3 = int(width/6), int(width*3/6), int(width*5/6)
    yind = np.where(flag[:, x1])[0];  iw1_y0, iw1_y1 = yind[0], yind[-1]
    yind = np.where(flag[:, x2])[0];  iw2_y0, iw2_y1 = yind[0], yind[-1]
    yind = np.where(flag[:, x3])[0];  iw3_y0, iw3_y1 = yind[0], yind[-1]

    # get xmin/max based on the non-overlap rows
    y0 = int((iw1_y0 + iw2_y0) / 2)
    y1 = int((iw1_y1 + iw2_y1) / 2)
    xs = [np.where(np.diff(flag[y0, :]))[0][0],
          np.where(np.diff(flag[y1, :]))[0][0]]
    iw2_x0, iw1_x1 = min(xs), max(xs)

    y0 = int((iw2_y0 + iw3_y0) / 2)
    y1 = int((iw2_y1 + iw3_y1) / 2)
    xs = [np.where(np.diff(flag[y0, :]))[0][-1],
          np.where(np.diff(flag[y1, :]))[0][-1]]
    iw3_x0, iw2_x1 = min(xs), max(xs)

    # iw1/2/3_x0/y0/x1/y1 --> box1/2/3
    box1 = [iw1_x0, iw1_y0, iw1_x1, iw1_y1]
    box2 = [iw2_x0, iw2_y0, iw2_x1, iw2_y1]
    box3 = [iw3_x0, iw3_y0, iw3_x1, iw3_y1]

    # adjust subswath overlap in X
    if cut_overlap_in_half:
        box2[0] += int((box1[2] - box2[0]) / 2)
        box3[0] += int((box2[2] - box3[0]) / 2)

    # initiate mask
    mask1 = np.zeros((length, width), dtype=np.bool_)
    mask2 = np.zeros((length, width), dtype=np.bool_)
    mask3 = np.zeros((length, width), dtype=np.bool_)

    # assign mask for each subswath
    mask1[box1[1]:box1[3], box1[0]:box1[2]] = 1
    mask2[box2[1]:box2[3], box2[0]:box2[2]] = 1
    mask3[box3[1]:box3[3], box3[0]:box3[2]] = 1
    mask1[mask2==1] = 0
    mask2[mask3==1] = 0

    return mask1, mask2, mask3, box1, box2, box3
