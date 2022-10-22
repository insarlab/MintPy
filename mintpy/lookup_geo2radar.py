############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Yunmeng Cao, Jul 2019                            #
############################################################


import os
import sys

import h5py
import numpy as np
from scipy.interpolate import griddata

try:
    from tqdm import tqdm
except ImportError:
    raise ImportError('Can not import tqdm!')

try:
    from concurrent.futures import ProcessPoolExecutor, as_completed
except ImportError:
    raise ImportError('Can not import concurrent!')

from mintpy.utils import readfile


################################################################################
def write_h5(datasetDict, out_file, metadata=None, ref_file=None, compression=None):

    if os.path.isfile(out_file):
        print(f'delete exsited file: {out_file}')
        os.remove(out_file)

    with h5py.File(out_file, 'w') as f:
        for dsName in datasetDict.keys():
            data = datasetDict[dsName]
            f.create_dataset(dsName,data=data,compression=compression)

        for key, value in metadata.items():
            f.attrs[key] = str(value)
            #print(key + ': ' +  value)
    print(f'finished writing to {out_file}')

    return out_file


def get_dataNames(FILE):
    with h5py.File(FILE, 'r') as f:
        dataNames = []
        for k0 in f.keys():
            dataNames.append(k0)
    return dataNames


def parallel_process(array, function, n_jobs=16, use_kwargs=False, front_num=1):
    """A parallel version of the map function with a progress bar.

        Args:
            array (array-like): An array to iterate over.
            function (function): A python function to apply to the elements of array
            n_jobs (int, default=16): The number of cores to use
            use_kwargs (boolean, default=False): Whether to consider the elements of array as dictionaries of
                keyword arguments to function
            front_num (int, default=3): The number of iterations to run serially before kicking off the parallel job.
                Useful for catching bugs
        Returns:
            [function(array[0]), function(array[1]), ...]
    """
    #We run the first few iterations serially to catch bugs
    if front_num > 0:
        front = [function(**a) if use_kwargs else function(a) for a in array[:front_num]]
    #If we set n_jobs to 1, just run a list comprehension. This is useful for benchmarking and debugging.
    if n_jobs==1:
        return front + [function(**a) if use_kwargs else function(a) for a in tqdm(array[front_num:])]
    #Assemble the workers
    with ProcessPoolExecutor(max_workers=n_jobs) as pool:
        #Pass the elements of array into function
        if use_kwargs:
            futures = [pool.submit(function, **a) for a in array[front_num:]]
        else:
            futures = [pool.submit(function, a) for a in array[front_num:]]
        kwargs = {
            'total': len(futures),
            'unit': 'it',
            'unit_scale': True,
            'leave': True
        }
        #Print out the progress as tasks complete
        for f in tqdm(as_completed(futures), **kwargs):
            del f
            #pass
    out = []
    #Get the results from the futures.
    for i, future in tqdm(enumerate(futures)):
        del i
        try:
            out.append(future.result())
        except Exception as e:
            out.append(e)
    return front + out


def split_range(N, M):
    #list0 = np.arange(0,N)
    dx = round(N/M)
    list00 = []
    for i in range(M):
        a0 = i*dx
        b0 = (i+1)*dx

        if b0 > N:
            b0 = N

        l0 = np.arange(a0,b0)
        list00.append(l0)

    return list00


def split_box(data,row_sample,col_sample):
    data_split= []
    row,col = data.shape
    list_row = split_range(row, row_sample)
    list_col = split_range(col, col_sample)

    for i in range(row_sample):
        for j in range(col_sample):
            y0 = min(list_row[i])
            y1 = max(list_row[i])
            x0 = min(list_col[j])
            x1 = max(list_col[j])

            data0 = data[y0:y1+1,x0:x1+1]
            data_split.append(data0)

    return data_split


def function(data0):
    points, zz1, zz2, grid_x0, grid_y0 = data0
    grid_lat0 = griddata(points, zz1, (grid_x0, grid_y0), method='nearest')
    grid_lon0 = griddata(points, zz2, (grid_x0, grid_y0), method='nearest')
    return grid_lat0, grid_lon0


def run_lookup_geo2radar(inps):
    """Convert the lookup table in geo-coordinates (from roipac, gamma) into radar-coordinates (from isce)."""

    rangeCoord = readfile.read(inps.geom_geo_file, datasetName = 'rangeCoord')[0].astype(np.float64)
    azimuthCoord = readfile.read(inps.geom_geo_file, datasetName = 'azimuthCoord')[0].astype(np.float64)
    #CPX_lt = complex(rangeCoord + '+' + azimuthCoord+'j')
    #CPX_lt = rangeCoord  + 1j *azimuthCoord

    meta_geo = readfile.read_attribute(inps.geom_geo_file)
    post_Lat = meta_geo['Y_STEP']
    post_Lon = meta_geo['X_STEP']
    Corner_LAT = meta_geo['Y_FIRST']
    Corner_LON = meta_geo['X_FIRST']

    if inps.write:
        meta = readfile.read_attribute(inps.write)
    elif inps.reference:
        meta = readfile.read_attribute(inps.reference)
    else:
        print('write_file or the reference_file should be provided at least one.')
        sys.exit(1)

    WIDTH_geo = int(meta_geo['WIDTH'])
    LENGTH_geo = int(meta_geo['LENGTH'])

    x = np.arange(0,WIDTH_geo)
    y = np.arange(0,LENGTH_geo)
    xv, yv = np.meshgrid(x, y)

    LAT = float(Corner_LAT) + yv*float(post_Lat)
    LON = float(Corner_LON) + xv*float(post_Lon)
    LAT = LAT.flatten()
    LON = LON.flatten()

    WIDTH  = int(meta['WIDTH'])
    LENGTH  = int(meta['LENGTH'])

    xx0 = rangeCoord.flatten()
    yy0 = azimuthCoord.flatten()

    zz01 = LAT.flatten()
    zz02 = LON.flatten()

    xx = xx0[xx0!=0]
    yy = yy0[xx0!=0]
    zz1 = zz01[xx0!=0] #lat
    zz2 = zz02[xx0!=0] # lon

    #points = (xx,yy)
    #points = np.zeros((len(xx),2))
    #points[:,0] = xx
    #points[:,1] = yy

    x = np.arange(0,WIDTH)
    y = np.arange(0,LENGTH)
    grid_x, grid_y = np.meshgrid(x, y)

    row_sample = 10
    col_sample = 10

    list_row = split_range(LENGTH, row_sample)
    list_col = split_range(WIDTH, col_sample)

    split_grid_y = split_box(grid_y,row_sample,col_sample)
    split_grid_x = split_box(grid_x,row_sample,col_sample)

    data_parallel = []
    for i, (ay, ax) in enumerate(zip(split_grid_y, split_grid_x)):
        # extend the search area by 5 pixels
        max_ax = max(ax.flatten()) + 5
        min_ax = min(ax.flatten()) - 5
        max_ay = max(ay.flatten()) + 5
        min_ay = min(ay.flatten()) - 5

        f0 = np.where((min_ax < xx) & (xx < max_ax) & (min_ay < yy) & (yy < max_ay))
        xx0 = xx[f0]
        yy0 = yy[f0]
        zz10 = zz1[f0]
        zz20 = zz2[f0]

        points0 = np.zeros((len(xx0),2))
        points0[:,0] = xx0
        points0[:,1] = yy0
        #print(split_grid_x[i].shape)

        data0 = (points0, zz10, zz20, ax, ay)
        data_parallel.append(data0)

    #grid_lat_all = []
    #grid_lon_all = []

    grid_lat = np.zeros((LENGTH,WIDTH), dtype=np.float32)
    grid_lon = np.zeros((LENGTH,WIDTH), dtype=np.float32)

    proNumb = inps.parallelNumb
    future = np.zeros((len(data_parallel),))
    future = list(future)
    future = parallel_process(data_parallel, function, n_jobs= proNumb, use_kwargs=False, front_num=1)

    for i in range(row_sample):
        for j in range(col_sample):
            k0 = i*col_sample + j
            kk = future[k0]
            y0 = min(list_row[i])
            y1 = max(list_row[i])
            x0 = min(list_col[j])
            x1 = max(list_col[j])
            #print(kk)
            try:
                lat0 = kk[0]
                lon0 = kk[1]
                grid_lat[y0:y1+1,x0:x1+1] = lat0
                grid_lon[y0:y1+1,x0:x1+1] = lon0
            except Exception as e:
                del e

    #grid_lat = griddata(points, zz1, (grid_x, grid_y), method='nearest')
    #grid_lon = griddata(points, zz2, (grid_x, grid_y), method='nearest')

    dataNames = get_dataNames(inps.write)
    datasetDict = dict()
    meta = readfile.read_attribute(inps.write)
    for k0 in dataNames:
        datasetDict[k0] = readfile.read(inps.write,datasetName = k0)[0]

    DEM = readfile.read(inps.write,datasetName = 'height')[0]
    grid_lat[DEM==0] = 0
    grid_lon[DEM==0] = 0
    grid_lat[grid_lat==0] = 'nan'
    grid_lon[grid_lon==0] = 'nan'
    datasetDict['latitude'] = grid_lat
    datasetDict['longitude'] = grid_lon
    write_h5(datasetDict, inps.write, metadata=meta, ref_file=None, compression=None)
    print('done.')

    return
