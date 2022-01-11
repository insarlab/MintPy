############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Aug 2021                           #
############################################################
# Recommand import:
#   from mintpy.utils import s1_utils


import os
import numpy as np
from mintpy.utils import ptime, time_func


def estimate_S1AB_bias(mintpy_dir, dates, ts_dis):
    """Estimate the bias between Sentinel-1 A and B.
    Parameters: mintpy_dir - str, path of the mintpy working directory
                dates      - list of datetime.datetime objects
                ts_dis     - 2D np.ndarray in size of (num_date, num_pixel) in float32
    Returns:    bias       - 1D np.ndarray in size of (num_pixel) in float32
                flagA/B    - 1D np.ndarray in size of (num_date) in bool
                dates_fit  - list of datetime.datetime objects
                ts_fitA/B  - 1D np.ndarray in size of (num_date_fit) in float32
    """
    num_date = len(dates)
    ts_dis = ts_dis.reshape(num_date, -1)

    # dates/flags for S1A/B
    date_listA = np.loadtxt(os.path.join(mintpy_dir, 'S1A_date.txt'), dtype=str).tolist()
    date_listB = np.loadtxt(os.path.join(mintpy_dir, 'S1B_date.txt'), dtype=str).tolist()
    date_list = sorted(date_listA + date_listB)
    min_date = date_listB[0]
    flagA = np.array([x in date_listA and x >= min_date for x in date_list], dtype=np.bool_)
    flagB = np.array([x in date_listB and x >= min_date for x in date_list], dtype=np.bool_)
    # update date_list to the shared time period only
    date_listA = np.array(date_list)[flagA].tolist()
    date_listB = np.array(date_list)[flagB].tolist()

    # fit
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

    return bias, flagA, flagB, dates_fit, ts_fitA, ts_fitB


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




