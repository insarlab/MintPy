############################################################
# Program is part of MintPy                                #
# Copyright(c) 2017-2019, Zhang Yunjun                     #
# Author:  Zhang Yunjun                                    #
############################################################
# Recommend import:
#   from mintpy.utils import variance as pvar


import numpy as np
import pyproj
import random
import matplotlib.pyplot as plt
from mintpy.utils import ptime


def sample_data(lat, lon, mask=None, num_sample=500):
    ''''''
    # Flatten input data
    for i in [lat, lon]:
        if len(i.shape) != 1:
            i = i.flatten()

    # Check number of samples and number of pixels
    num_pixel = len(lat)
    if num_sample > num_pixel:
        print('Number of samples > number of pixels, fix number of samples to number of pixels.')
        num_sample = num_pixel

    # Check input mask
    if mask is None:
        mask = np.ones((num_pixel))

    # Random select samples
    idx = np.arange(num_pixel)
    idx_sample = random.sample(idx[mask.flatten() == 1.0], int(num_sample))

    lat_sample = lat[idx_sample]
    lon_sample = lon[idx_sample]

    return idx_sample, lat_sample, lon_sample


def get_distance(lat, lon, i):
    '''Return the distance of all points in lat/lon from its ith point'''
    lat1 = lat[i]*np.ones(lat.shape)
    lon1 = lon[i]*np.ones(lon.shape)

    g = pyproj.Geod(ellps='WGS84')
    dist = g.inv(lon1, lat1, lon, lat)[2]
    return dist


def structure_function(data, lat, lon, step=5e3, min_pair_num=100e3, print_msg=True):
    num_sample = len(data)
    distance = np.zeros((num_sample**2))
    variance = np.zeros((num_sample**2))
    if print_msg:
        prog_bar = ptime.progressBar(maxValue=num_sample)
    for i in range(num_sample):
        distance[i*num_sample:(i+1)*num_sample] = get_distance(lat, lon, i)
        variance[i*num_sample:(i+1)*num_sample] = np.square(data - data[i])
        if print_msg:
            prog_bar.update(i+1, every=10)
    if print_msg:
        prog_bar.close()

    dist, std, stdStd = bin_variance(
        distance, variance, step=step, min_pair_num=min_pair_num, print_msg=print_msg)
    return dist, std, stdStd


def bin_variance(distance, variance, step=5e3, min_pair_num=100e3, print_msg=True):
    x_steps = np.arange(0, np.max(distance), step)
    num_step = len(x_steps)
    std = np.zeros(x_steps.shape)
    stdStd = np.zeros(std.shape)
    p_num = np.zeros(x_steps.shape)

    if print_msg:
        prog_bar = ptime.progressBar(maxValue=num_step)
    for i in range(num_step):
        x = x_steps[i]
        idx = (distance > max(0, x-step/2.)) * (distance < x+step/2.)
        p_num[i] = np.sum(idx)
        std[i] = np.mean(np.sqrt(variance[idx]))
        stdStd[i] = np.std(np.sqrt(variance[idx]))
        if print_msg:
            prog_bar.update(i+1, every=10)
    if print_msg:
        prog_bar.close()

    max_step_idx = int(max(np.argwhere(p_num > min_pair_num)))
    return x_steps[0:max_step_idx], std[0:max_step_idx], stdStd[0:max_step_idx]
