#! /usr/bin/env python
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Bhuvan Varugu, Zhang Yunjun, 2018                #
############################################################


import os
import sys
import argparse
import h5py
import numpy as np
from scipy.interpolate import griddata, RegularGridInterpolator as RGI
from mintpy.utils import ptime, readfile, utils as ut, network as pnet
import multiprocessing as mp
##########################################################

def subset_data(latmin, latmax, lonmin, lonmax, work_dir, gacos_file):

    command = 'subset.py {} -l {} {} -L {} {} -o {}/{}'.format(gacos_file, latmin,
                                                               latmax, lonmin, lonmax,
                                                               work_dir, os.path.basename(gacos_file))
    print(command)
    os.system(command)
    return


def read_params(filename):
    fil = open(filename + '.rsc', 'r')
    line = fil.readline()
    rscdict = {}
    while line:
        llist = line.split()
        if len(llist) > 0:
            rscdict[llist[0]] = llist[1]
        line = fil.readline()
    fil.close()
    nx = np.int(rscdict['WIDTH'])
    ny = np.int(rscdict['FILE_LENGTH'])
    lat = np.zeros((4, 1))
    lon = np.zeros((4, 1))
    lat[0] = np.float(rscdict['Y_FIRST'])
    lat[1] = np.float(rscdict['Y_FIRST'])
    lat[2] = np.float(rscdict['Y_FIRST']) + (ny - 1) * np.float(rscdict['Y_STEP'])
    lat[3] = np.float(rscdict['Y_FIRST']) + (ny - 1) * np.float(rscdict['Y_STEP'])
    lon[0] = np.float(rscdict['X_FIRST'])
    lon[1] = np.float(rscdict['X_FIRST']) + (nx - 1) * np.float(rscdict['X_STEP'])
    lon[2] = np.float(rscdict['X_FIRST'])
    lon[3] = np.float(rscdict['X_FIRST']) + (nx - 1) * np.float(rscdict['X_STEP'])
    return lat, lon, nx, ny


def get_delay(input_args):
    delay_file, atr, lookup_table, cinc, date, work_dir = input_args
    pts_new = lookup_table[0]
    length = lookup_table[1]
    width = lookup_table[2]

    [lat, lon, nx, ny] = read_params(delay_file)
    data1 = np.fromfile(delay_file, dtype=np.float32, sep=(""))
    data = data1.reshape(ny, nx)
    data = np.flipud(data)
    lats, step = np.linspace(lat[2], lat[0], num=ny, endpoint=True, retstep=True)
    lats = np.asarray(lats)
    lons, step = np.linspace(lon[0], lon[1], num=nx, endpoint=True, retstep=True)
    lons = np.asarray(lons)
    pts_old = tuple((lats.flatten(), lons.flatten()))

    rx = np.int(atr['REF_X'])
    ry = np.int(atr['REF_Y'])
    if 'X_FIRST' in list(atr.keys()):
        RGI_func = RGI(pts_old, data, method='linear', bounds_error=False)
        delay_geo = RGI_func(pts_new)
        delay_geo = np.flipud(delay_geo.reshape(length, width))
        #delay_geo = delay_geo - delay_geo[ry, rx]
        delay_geo = delay_geo / cinc
        delay = delay_geo
        del delay_geo, pts_new

    else:
        if atr['processor'] == 'isce':
            print('interferogram processor is isce')
            RGI_func = RGI(pts_old, data, method='linear', bounds_error=False)
            delay_rdr = RGI_func(pts_new)
            delay_rdr = delay_rdr.reshape(length, width)
            #delay_rdr = delay_rdr - delay_rdr[ry, rx]
            delay_rdr = delay_rdr / cinc
            delay = delay_rdr
            del delay_rdr, pts_new

        else:
            print('interferogram processor is roi_pac')
            RGI_func = RGI(pts_old, data, method='linear', bounds_error=False)
            pts_new2 = lookup_table[3]
            idx = lookup_table[4]
            delay_geo = RGI_func(pts_new)
            delay_geo = delay_geo.reshape(length, width)
            yy, xx = np.mgrid[0:length:length * 1j, 0:width:width * 1j]
            yx_rdr = np.hstack((yy.reshape(-1, 1), xx.reshape(-1, 1)))
            delay_geo = delay_geo[idx]
            delay_rdr = griddata(pts_new2, delay_geo, yx_rdr, method='linear')
            delay_rdr = delay_rdr.reshape(length, width)
            #delay_rdr = delay_rdr - delay_rdr[ry, rx]
            delay_rdr = delay_rdr / cinc
            delay = delay_rdr
            del delay_rdr, pts_new, pts_new2

    np.save(os.path.join(work_dir, date + '.npy'), delay)
    print(os.path.basename(delay_file) + ' Done')
    return


def lookup_table_read(atr, lookup_file):
    if 'X_FIRST' in list(atr.keys()):
        latd = np.zeros((4, 1))
        lond = np.zeros((4, 1))
        width = np.int(atr['WIDTH'])
        length = np.int(atr['FILE_LENGTH'])
        latd[0] = np.float(atr['Y_FIRST'])
        latd[2] = np.float(atr['Y_FIRST']) + (length - 1) * np.float(atr['Y_STEP'])
        lond[0] = np.float(atr['X_FIRST'])
        lond[1] = np.float(atr['X_FIRST']) + (width - 1) * np.float(atr['X_STEP'])
        xarr, step = np.linspace(lond[0], lond[1], num=width, endpoint=True, retstep=True)
        yarr, step = np.linspace(latd[2], latd[0], num=length, endpoint=True, retstep=True)
        lons3 = np.tile(xarr, length)
        lats3 = np.repeat(yarr, width)
        pts_new = np.hstack((lats3.reshape(-1, 1), lons3.reshape(-1, 1)))
        return pts_new, length, width

    else:
        if atr['processor'] == 'isce':
            width = np.int(atr['WIDTH'])
            length = np.int(atr['FILE_LENGTH'])
            lon_lut = readfile.read(lookup_file, datasetName='longitude')[0]
            lat_lut = readfile.read(lookup_file, datasetName='latitude')[0]
            pts_new = np.hstack((lat_lut.reshape(-1, 1), lon_lut.reshape(-1, 1)))
            return pts_new, length, width

        else:
            print('interferogram processor is roi_pac')
            [latd, lond, width, length] = read_params(lookup_file)
            xarr, step = np.linspace(lond[0], lond[1], num=nxd, endpoint=True, retstep=True)
            yarr, step = np.linspace(latd[2], latd[0], num=nyd, endpoint=True, retstep=True)
            lons3 = np.tile(xarr, length)
            lats3 = np.repeat(yarr, width)
            pts_new = np.hstack((lats3.reshape(-1, 1), lons3.reshape(-1, 1)))
            rg = readfile.read(lookup_file, datasetName='range')[0]
            az = readfile.read(lookup_file, datasetName='azimuth')[0]
            idx = (az > 0.0) * (az <= length) * (rg > 0.0) * (rg <= width)
            pts_new2 = np.hstack((az[idx].reshape(-1, 1), rg[idx].reshape(-1, 1)))
            return pts_new, length, width, pts_new2, idx


###############################################################
EXAMPLE = '''example:
  tropo_gacos.py timeseries.h5 -l geomap_*rlks.trans -i incidenceAngle.h5
'''
TEMPLATE = '''
mintpy.troposphericDelay.method        = GACOS   #[pyaps, height_correlation,GACOS] 
'''


def cmdLineParse():
    parser = argparse.ArgumentParser(description='Tropospheric correction using GACOS delays\n', \
                                     formatter_class=argparse.RawTextHelpFormatter, \
                                     epilog=EXAMPLE)

    parser.add_argument(dest='timeseries_file', nargs='?', help='timeseries HDF5 file, i.e. timeseries.h5')
    parser.add_argument('-i', dest='inc_angle', \
                        help='a file containing all incidence angles, or a number representing for the whole image.')
    parser.add_argument('-l', dest='lookup_file', \
                        help='a file containing all information to tranfer from radar to geo coordinates.')
    parser.add_argument('--GACOS-dir', dest='GACOS_dir', \
                        help='directory to downloaded GACOS delays data, i.e. ./../WEATHER/GACOS\n' + \
                             'use directory of input timeseries_file if not specified.')
    parser.add_argument('--date-list', dest='date_list_file', \
                        help='Read the first column of text file as list of date to download data\n' + \
                             'in YYYYMMDD or YYMMDD format')
    parser.add_argument('--ref-yx', dest='ref_yx', type=int, nargs=2, help='reference pixel in y/x')
    parser.add_argument('--template', dest='template_file', \
                        help='template file with input options below:\n' + TEMPLATE)
    parser.add_argument('-o', dest='out_file', help='Output file name for trospheric corrected timeseries.')

    inps = parser.parse_args()

    # Correcting TIMESERIES or DOWNLOAD DATA ONLY, required one of them
    if not inps.timeseries_file and not inps.download:
        parser.print_help()
        sys.exit(1)
    return inps


def main(argv):
    inps = cmdLineParse()

    if inps.timeseries_file:
        inps.timeseries_file = ut.get_file_list([inps.timeseries_file])[0]
        atr = readfile.read_attribute(inps.timeseries_file)
        k = atr['FILE_TYPE']
        if 'REF_Y' not in list(atr.keys()) and inps.ref_yx:
            print('No reference info found in input file, use input ref_yx: ' + str(inps.ref_yx))
            atr['REF_Y'] = inps.ref_yx[0]
            atr['REF_X'] = inps.ref_yx[1]
    
    # ****reading incidence angle file***/
    if os.path.isfile(inps.inc_angle):
        inps.inc_angle = readfile.read(inps.inc_angle, datasetName='incidenceAngle')[0]
        inps.inc_angle = np.nan_to_num(inps.inc_angle)
    else:
        inps.inps.inc_angle = float(inps.inc_angle)
        print('incidence angle: ' + str(inps.inc_angle))
    cinc = np.cos(inps.inc_angle * np.pi / 180.0)
    
    # ****look up file****/
    if inps.lookup_file:
        inps.lookup_file = ut.get_file_list([inps.lookup_file])[0]  # 'geomap_32rlks_tight.trans'
    ''' 
    altitude = float(atr['altitude'])
    slantRangeDistance = readfile.read(inps.lookup_file, datasetName='slantRangeDistance')[0]
    cinc = altitude/slantRangeDistance
    '''
    # ****GACOS****/
    delay_source = 'GACOS'
    # Get weather directory
    if not inps.GACOS_dir:
        if inps.timeseries_file:
            inps.GACOS_dir = os.path.dirname(os.path.abspath(inps.timeseries_file)) + '/../WEATHER/GACOS'
        elif inps.lookup_file:
            inps.GACOS_dir = os.path.dirname(os.path.abspath(inps.lookup_file)) + '/../WEATHER/GACOS'
        else:
            inps.GACOS_dir = os.path.abspath(os.getcwd())

    print('Store weather data into directory: ' + inps.GACOS_dir)

    # source_dir=os.path.dirname(os.path.abspath('timeseries_file'))+'/Agung/GACOS/data';print source_dir
    # os.makedirs(GACOS_dir)  -----------------------------------------------add part to copy/download weather data------#
    # ----get date list-----#
    if not inps.date_list_file:
        print('read date list info from: ' + inps.timeseries_file)
        h5 = h5py.File(inps.timeseries_file, 'r')
        if 'timeseries' in list(h5.keys()):
            date_list = sorted([i.decode('utf8') for i in h5['date'][:]])
        elif k in ['interferograms', 'coherence', 'wrapped']:
            ifgram_list = sorted([i.decode('utf8') for i in h5['date'][:]])
            date12_list = pnet.get_date12_list(inps.timeseries_file)
            m_dates = [i.split('-')[0] for i in date12_list]
            s_dates = [i.split('-')[1] for i in date12_list]
            date_list = ptime.yyyymmdd(sorted(list(set(m_dates + s_dates))))
        else:
            raise ValueError('Un-support input file type:' + k)
        h5.close()
    else:
        date_list = ptime.yyyymmdd(np.loadtxt(inps.date_list_file, dtype=str, usecols=(0,)).tolist())
        print('read date list info from: ' + inps.date_list_file)

    # ****cheacking availability of delays****/
    print('checking availability of delays')
    delay_file_list = []
    for d in date_list:
        if delay_source == 'GACOS':  delay_file = inps.GACOS_dir + '/' + d + '.ztd';
        delay_file_list.append(delay_file)
    delay_file_existed = ut.get_file_list(delay_file_list)

    if len(delay_file_existed) == len(date_list):
        print('no missing files')
    else:
        print('no. of date files found:', len(delay_file_existed));
        print('no. of dates:', len(date_list))

    # *****Calculating delays***/
    print('calculating delays')
    length = int(atr['FILE_LENGTH'])
    width = int(atr['WIDTH'])

    date_num = len(date_list)
    lookup_table = lookup_table_read(atr, inps.lookup_file)

    work_dir = os.path.dirname(inps.timeseries_file) + '/gacos_ztd'

    os.makedirs(work_dir, exist_ok=True)

    latmin = -4.5762 #np.min(lookup_table[0][0])

    latmax = 1.6733  #np.max(lookup_table[0][0])
    lonmin = -80.0468 #np.min(lookup_table[0][1])
    lonmax = -76.4248 #np.max(lookup_table[0][1])

    input_arg = []
    for i in range(date_num):
        delay_file = delay_file_existed[i]
        subset_data(latmin, latmax, lonmin, lonmax, work_dir, delay_file)
        delay_file = os.path.join(work_dir, os.path.basename(delay_file))
        date = date_list[i]
        #delay, date0 = get_delay((delay_file, atr, lookup_table, cinc, date, work_dir))
        input_arg.append((delay_file, atr, lookup_table, cinc, date, work_dir))

    pool = mp.Pool(40)
    pool.map(get_delay, input_arg)
    pool.close()

    # Write tropospheric delay to HDF5
    if inps.out_file:
        tropFile = inps.out_file
    else:
        tropFile = 'GACOSdelays' + '.h5'

    print('writing >>> %s' % (tropFile))
    h5trop = h5py.File(tropFile, 'w')
    trop_ts = h5trop.create_dataset('timeseries',
                                    shape=(date_num, length, width),
                                    maxshape=(None, length, width),
                                    chunks=True,
                                    dtype=np.float32)

    # reading wrf files for each epoch and getting delay
    prog_bar = ptime.progressBar(maxValue=date_num)

    ii = 0
    for date in date_list:
        delay = np.load(os.path.join(work_dir, date + '.npy'))
        i = date_list.index(date)
        trop_ts[i, :, :] = delay
        prog_bar.update(ii + 1, suffix=date)
        ii += 1

    prog_bar.close()

    print('Delays Calculated')
    # Convert relative phase delay on reference date
    try:
        ref_date = atr['REF_DATE']
    except:
        ref_date = date_list[0]

    print('convert to relative phase delay with reference date: ' + ref_date)
    ref_idx = date_list.index(ref_date)
    #trop_ts[:, :, :] -= np.tile(trop_ts[ref_idx, :, :], (date_num, 1, 1))

    dates = np.array(date_list, dtype=np.string_)
    h5trop.create_dataset('date', data=dates)

    with h5py.File(inps.timeseries_file, 'r') as f:
        bperp = f['bperp'][:]
    bperp = np.array(bperp, dtype=np.float32)
    h5trop.create_dataset('bperp', data=bperp)

    # remove metadata related with double reference
    # because absolute delay is calculated and saved
    for key in ['REF_DATE','REF_X','REF_Y','REF_LAT','REF_LON']:
        if key in atr.keys():
            atr.pop(key)

    # Write Attributes
    for key, value in atr.items():
        h5trop.attrs[key] = value

    h5trop.close()

    print('finished')
    return


###############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
