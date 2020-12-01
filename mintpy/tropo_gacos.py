#! /usr/bin/env python
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Bhuvan Varugu, Zhang Yunjun, 2018                #
# Update: Sara Mirzaee, 2020                               #
############################################################


import os
import sys
import argparse
import h5py
import numpy as np
from skimage.transform import resize
from scipy.interpolate import griddata, RegularGridInterpolator as RGI
from mintpy.utils import ptime, readfile, utils as ut, network as pnet
import multiprocessing as mp


##########################################################


def get_dataset_size(fname):
    atr = readfile.read_attribute(fname)
    return (atr['LENGTH'], atr['WIDTH'])


def correct_timeseries(dis_file, tropo_file, cor_dis_file):
    # diff.py can handle different reference in space and time
    # between the absolute tropospheric delay and the double referenced time-series
    print('\n------------------------------------------------------------------------------')
    print('correcting relative delay for input time-series using diff.py')
    from mintpy import diff

    iargs = [dis_file, tropo_file, '-o', cor_dis_file, '--force']
    print('diff.py', ' '.join(iargs))
    diff.main(iargs)
    return cor_dis_file


def correct_single_ifgram(dis_file, tropo_file, cor_dis_file):
    print('\n------------------------------------------------------------------------------')
    print('correcting relative delay for input interferogram')

    print('read data from {}'.format(dis_file))
    data, atr = readfile.read(dis_file, datasetName='phase')
    date1, date2 = ptime.yyyymmdd(atr['DATE12'].split('-'))

    print('calc tropospheric delay for {}-{} from {}'.format(date1, date2, tropo_file))
    tropo = readfile.read(tropo_file, datasetName=date2)[0]
    tropo -= readfile.read(tropo_file, datasetName=date1)[0]
    tropo *= -4. * np.pi / float(atr['WAVELENGTH'])

    print('write corrected data to {}'.format(cor_dis_file))
    writefile.write(data - tropo, cor_dis_file, atr)
    return cor_dis_file


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


def get_delay(delay_file, atr, cinc, pts_new=None):
    length = int(atr['FILE_LENGTH'])
    width = int(atr['WIDTH'])

    if 'X_FIRST' in list(atr.keys()):
        pixel_box = (0, 0, int(atr['WIDTH']), int(atr['LENGTH']))
        geo_box = ut.coordinate(atr).box_pixel2geo(pixel_box)
        atr_ztd = readfile.read_attribute(delay_file)
        coord = ut.coordinate(atr_ztd)
        box = coord.box_geo2pixel(geo_box)
        phs = readfile.read(delay_file, box=box)[0]
        out_shape = (length, width)
        delay_geo = resize(phs, out_shape, order=1, mode='constant', anti_aliasing=True, preserve_range=True)
        delay_geo = delay_geo / cinc
        return delay_geo

    else:
        [lat, lon, nx, ny] = read_params(delay_file)
        data1 = np.fromfile(delay_file, dtype=np.float32, sep=(""))
        data = data1.reshape(ny, nx)
        data = np.flipud(data)
        lats, step = np.linspace(lat[2], lat[0], num=ny, endpoint=True, retstep=True)
        lats = np.asarray(lats)
        lons, step = np.linspace(lon[0], lon[1], num=nx, endpoint=True, retstep=True)
        lons = np.asarray(lons)
        pts_old = tuple((lats.flatten(), lons.flatten()))

        RGI_func = RGI(pts_old, data, method='linear', bounds_error=False)
        delay_rdr = RGI_func(pts_new)
        delay_rdr = delay_rdr.reshape(length, width)
        delay_rdr = delay_rdr / cinc
        del pts_new
        return delay_rdr


def lookup_table_read(lookup_file):
    lon_lut = readfile.read(lookup_file, datasetName='longitude')[0]
    lat_lut = readfile.read(lookup_file, datasetName='latitude')[0]
    pts_new = np.hstack((lat_lut.reshape(-1, 1), lon_lut.reshape(-1, 1)))
    return pts_new


###############################################################
EXAMPLE = '''example:
  tropo_gacos.py -f timeseries.h5 -l inputs/geometryRadar.h5 --GACOS_dir ./GACOS
  tropo_gacos.py -f geo/geo_timeseries.h5 -l geo/geo_geometryRadar.h5 --GACOS_dir ./GACOS

'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Tropospheric correction using GACOS delays\n',
                                     formatter_class=argparse.RawTextHelpFormatter, epilog=EXAMPLE)

    parser.add_argument('-f', '--file', dest='dis_file', help='timeseries HDF5 file, i.e. timeseries.h5')
    parser.add_argument('-l', '--lookup', dest='lookup_file',
                        help='a file containing all information to tranfer from radar to geo coordinates.')
    parser.add_argument('--GACOS-dir', dest='GACOS_dir',
                        help='directory to downloaded GACOS delays data, i.e. ./GACOS\n' +
                             'use directory of input timeseries_file if not specified.')
    parser.add_argument('--date-list', dest='date_list_file',
                        help='Read the first column of text file as list of date to download data\n' +
                             'in YYYYMMDD or YYMMDD format')
    parser.add_argument('--ref-yx', dest='ref_yx', type=int, nargs=2, help='reference pixel in y/x')
    parser.add_argument('-o', dest='out_file', help='Output file name for trospheric corrected timeseries.')

    inps = parser.parse_args()

    # Correcting TIMESERIES or DOWNLOAD DATA ONLY, required one of them
    if not inps.dis_file and not inps.download:
        parser.print_help()
        sys.exit(1)
    return inps


def main(argv):
    inps = cmdLineParse()

    if inps.dis_file:
        inps.dis_file = ut.get_file_list([inps.dis_file])[0]
        atr = readfile.read_attribute(inps.dis_file)
        k = atr['FILE_TYPE']
        if 'REF_Y' not in list(atr.keys()) and inps.ref_yx:
            print('No reference info found in input file, use input ref_yx: ' + str(inps.ref_yx))
            atr['REF_Y'] = inps.ref_yx[0]
            atr['REF_X'] = inps.ref_yx[1]

    # ****reading incidence angle file***/
    inps.inc_angle = readfile.read(inps.lookup_file, datasetName='incidenceAngle')[0]
    cinc = np.cos(inps.inc_angle * np.pi / 180.0)

    # ****look up file****/
    if inps.lookup_file:
        inps.lookup_file = ut.get_file_list([inps.lookup_file])[0]

    # ****GACOS****/
    # Get weather directory
    if not inps.GACOS_dir:
        if inps.dis_file:
            inps.GACOS_dir = os.path.dirname(os.path.abspath(inps.dis_file)) + '/GACOS'
        elif inps.lookup_file:
            inps.GACOS_dir = os.path.dirname(os.path.abspath(inps.lookup_file)) + '/GACOS'
        else:
            inps.GACOS_dir = os.path.abspath(os.getcwd())

    print('Store weather data into directory: ' + inps.GACOS_dir)

    # source_dir=os.path.dirname(os.path.abspath('inps.dis_file'))+'/Agung/GACOS/data';print source_dir
    # os.makedirs(GACOS_dir)  -----------------------------------------------add part to copy/download weather data------#
    # ----get date list-----#
    if not inps.date_list_file:
        print('read date list info from: ' + inps.dis_file)
        h5 = h5py.File(inps.dis_file, 'r')
        if 'timeseries' in list(h5.keys()):
            date_list = sorted([i.decode('utf8') for i in h5['date'][:]])
        elif k in ['interferograms', 'coherence', 'wrapped']:
            ifgram_list = sorted([i.decode('utf8') for i in h5['date'][:]])
            date12_list = pnet.get_date12_list(inps.dis_file)
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
    tropFile = os.path.dirname(os.path.abspath(inps.dis_file)) + '/GACOS.h5'
    if os.path.isfile(tropFile) and get_dataset_size(tropFile) == get_dataset_size(inps.dis_file):
        print(tropFile + ' exists ...')
    else:
        delay_file_list = []
        for d in date_list:
            delay_file = inps.GACOS_dir + '/' + d + '.ztd'
            delay_file_list.append(delay_file)
        delay_file_existed = ut.get_file_list(delay_file_list)

        if len(delay_file_existed) == len(date_list):
            print('no missing files')
        else:
            print('no. of date files found:', len(delay_file_existed))
            print('no. of dates:', len(date_list))

        # *****Calculating delays***/
        print('calculating delays')
        length = int(atr['FILE_LENGTH'])
        width = int(atr['WIDTH'])

        date_num = len(date_list)

        # Write tropospheric delay to HDF5
        print('writing >>> %s' % (tropFile))

        h5trop = h5py.File(tropFile, 'w')
        trop_ts = h5trop.create_dataset('timeseries',
                                        shape=(date_num, length, width),
                                        maxshape=(None, length, width),
                                        chunks=True,
                                        dtype=np.float32)

        prog_bar = ptime.progressBar(maxValue=date_num)

        if 'X_FIRST' in list(atr.keys()):
            pts_new = None
        else:
            pts_new = lookup_table_read(inps.lookup_file)

        for i in range(date_num):
            delay_file = delay_file_existed[i]
            date = date_list[i]
            delay = get_delay(delay_file, atr, cinc, pts_new)
            trop_ts[i, :, :] = delay
            prog_bar.update(i + 1, suffix=date)

        dates = np.array(date_list, dtype=np.string_)
        h5trop.create_dataset('date', data=dates)

        with h5py.File(inps.dis_file, 'r') as f:
            bperp = f['bperp'][:]
        bperp = np.array(bperp, dtype=np.float32)
        h5trop.create_dataset('bperp', data=bperp)

        # remove metadata related with double reference
        # because absolute delay is calculated and saved
        for key in ['REF_DATE', 'REF_X', 'REF_Y', 'REF_LAT', 'REF_LON']:
            if key in atr.keys():
                atr.pop(key)

        # Write Attributes
        for key, value in atr.items():
            h5trop.attrs[key] = value

        h5trop.close()

    # correct tropo delay from displacement time-series
    if not inps.out_file:
        inps.out_file = inps.dis_file.split('.')[0] + '_GACOS.h5'

    if inps.dis_file:
        ftype = atr['FILE_TYPE']
        if ftype == 'timeseries':
            correct_timeseries(dis_file=inps.dis_file,
                               tropo_file=tropFile,
                               cor_dis_file=inps.out_file)

        elif ftype == '.unw':
            correct_single_ifgram(dis_file=inps.dis_file,
                                  tropo_file=tropFile,
                                  cor_dis_file=inps.out_file)
        else:
            print('input file {} is not timeseries nor .unw, correction is not supported yet.'.format(ftype))

    else:
        print('No input displacement file, skip correcting tropospheric delays.')

    return


###############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
