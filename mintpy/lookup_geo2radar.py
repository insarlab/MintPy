############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Yunmeng Cao, Jul 2019                            #
############################################################


import os
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


################################################################################
def write_h5(datasetDict, out_file, metadata=None, ref_file=None, compression=None):

    if os.path.isfile(out_file):
        print('delete exsited file: {}'.format(out_file))
        os.remove(out_file)

    with h5py.File(out_file, 'w') as f:
        for dsName in datasetDict.keys():
            data = datasetDict[dsName]
            f.create_dataset(dsName,data=data,compression=compression)

        for key, value in metadata.items():
            f.attrs[key] = str(value)
            #print(key + ': ' +  value)
    print('finished writing to {}'.format(out_file))

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
