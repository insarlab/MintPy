############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import numpy as np

from mintpy.utils import readfile
from mintpy.multilook import multilook_data
from mintpy.mask import mask_matrix


############################################################################
def design_matrix(dem, poly_order=1):
    """Design matrix for phase/elevation ratio estimation
    Parameters: dem : 1D array in size of (length*width, ), or
                      2D array in size of (length, width)
                poly_order : int
    Returns:    A : 2D array in size of (length*width, poly_order+1)
    """
    dem = np.reshape(dem, (-1, 1))
    A = np.ones((dem.size, 1), np.float64)
    for i in range(poly_order):
        Ai = np.array(dem**(i+1), np.float64)
        A = np.hstack((A, Ai))
    return A


def read_topographic_data(geom_file, metadata):
    print('read DEM from file: '+geom_file)
    dem = readfile.read(geom_file,
                        datasetName='height',
                        print_msg=False)[0]

    print('considering the incidence angle of each pixel ...')
    inc_angle = readfile.read(geom_file,
                              datasetName='incidenceAngle',
                              print_msg=False)[0]
    dem *= 1.0/np.cos(inc_angle*np.pi/180.0)

    ref_y = int(metadata['REF_Y'])
    ref_x = int(metadata['REF_X'])
    dem -= dem[ref_y, ref_x]

    # Design matrix for elevation v.s. phase
    # dem = dem.flatten()
    return dem


def estimate_phase_elevation_ratio(dem, ts_data, inps):
    """Estimate phase/elevation ratio for each acquisition of timeseries
    Parameters: dem     : 2D array in size of (          length, width)
                ts_data : 3D array in size of (num_date, length, width)
                inps    : Namespace
    Returns:    X       : 2D array in size of (poly_num+1, num_date)
    """
    num_date = ts_data.shape[0]

    # prepare phase and elevation data
    print('reading mask from file: '+inps.mask_file)
    mask = readfile.read(inps.mask_file, datasetName='mask')[0]
    dem = mask_matrix(np.array(dem), mask)
    ts_data = mask_matrix(np.array(ts_data), mask)

    # display
    # 1. effect of multilooking --> narrow phase range --> better ratio estimation
    debug_mode = False
    if debug_mode:
        import matplotlib.pyplot as plt
        d_index = 47   # np.argmax(topo_trop_corr)
        data = ts_data[d_index, :, :]
        title = inps.date_list[d_index]
        plt.figure()
        plt.plot(dem[~np.isnan(dem)],
                 data[~np.isnan(dem)],
                 '.', label='Number of Looks = 1')
        mli_dem = multilook_data(dem, 8, 8)
        mli_data = multilook_data(data, 8, 8)
        plt.plot(mli_dem[~np.isnan(mli_dem)],
                 mli_data[~np.isnan(mli_dem)],
                 '.', label='Number of Looks = 8')
        plt.legend()
        plt.xlabel('Elevation (m)')
        plt.ylabel('Range Change (m)')
        plt.title(title)
        out_file = 'phase_elevation_ratio_{}.png'.format(title)
        plt.savefig(out_file, bbox_inches='tight', transparent=True, dpi=300)
        print('save to {}'.format(out_file))
        plt.show()

    print('----------------------------------------------------------')
    print('Empirical tropospheric delay correction based on phase/elevation ratio (Doin et al., 2009)')
    print('polynomial order: {}'.format(inps.poly_order))

    if inps.num_multilook > 1:
        print('number of multilook: {} (multilook data for estimation only)'.format(inps.num_multilook))
        mask = multilook_data(mask, inps.num_multilook, inps.num_multilook)
        dem = multilook_data(dem, inps.num_multilook, inps.num_multilook)
        ts_data = multilook_data(ts_data, inps.num_multilook, inps.num_multilook)

    if inps.threshold > 0.:
        print('correlation threshold: {}'.format(inps.threshold))

    mask_nan = ~np.isnan(dem)
    dem = dem[mask_nan]
    ts_data = ts_data[:, mask_nan]

    # calculate correlation coefficient
    print('----------------------------------------------------------')
    print('calculate correlation of DEM with each acquisition')
    topo_trop_corr = np.zeros(num_date, np.float32)
    for i in range(num_date):
        phase = ts_data[i, :]
        cc = 0.
        if np.count_nonzero(phase) > 0:
            comp_data = np.vstack((dem, phase))
            cc = np.corrcoef(comp_data)[0, 1]
            topo_trop_corr[i] = cc
        print('{}: {:>5.2f}'.format(inps.date_list[i], cc))
    topo_trop_corr = np.abs(topo_trop_corr)
    print('average correlation magnitude: {:>5.2f}'.format(np.nanmean(topo_trop_corr)))

    # estimate ratio parameter
    print('----------------------------------------------------------')
    print('estimate phase/elevation ratio')
    A = design_matrix(dem=dem, poly_order=inps.poly_order)
    X = np.dot(np.linalg.pinv(A), ts_data.T)
    X = np.array(X, dtype=np.float32)
    X[:, topo_trop_corr < inps.threshold] = 0.
    return X


def estimate_tropospheric_delay(dem, X, metadata):
    poly_order = X.shape[0]-1
    num_date = X.shape[1]
    length, width = dem.shape

    print('estimate the stratified tropospheric delay')
    B = design_matrix(dem=dem, poly_order=poly_order)
    trop_data = np.array(np.dot(B, X).T, dtype=np.float32)

    ref_index = int(metadata['REF_Y']) * width + int(metadata['REF_X'])
    ref_value = trop_data[:, ref_index].reshape(-1, 1)
    trop_data -= np.tile(ref_value, (1, length*width))

    trop_data = np.reshape(trop_data, (num_date, length, width))
    return trop_data
