#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2015, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os
import sys
import re

try:
    import pyaps as pa
except:
    sys.exit('Cannot import pyaps into Python!')

import argparse
import h5py
import numpy as np
from pysar.objects import timeseries, geometry
from pysar.utils import readfile, writefile, ptime, utils as ut

standardWeatherModelNames = {'ERAI': 'ECMWF', 'ERAINT': 'ECMWF', 'ERAINTERIM': 'ECMWF',
                             'MERRA2': 'MERRA'
                             }


###############################################################
EXAMPLE = """example:
  tropcor_pyaps.py -d 20151002 20151003 --hour 12 -m ECMWF
  tropcor_pyaps.py -d date_list.txt     --hour 12 -m MERRA
  tropcor_pyaps.py -d 20151002 20151003 --hour 12 -m ECMWF --dem geometryRadar.h5 --ref-yx 30 40 -i geometryRadar.h5
  tropcor_pyaps.py -m ECMWF --dem geometryRadar.h5 -i geometryRadar.h5 -f timeseries.h5
"""

REFERENCE = """reference:
  Jolivet, R., R. Grandin, C. Lasserre, M.-P. Doin and G. Peltzer (2011), Systematic InSAR tropospheric
  phase delay corrections from global meteorological reanalysis data, Geophys. Res. Lett., 38, L17311,
  doi:10.1029/2011GL048757
"""

TEMPLATE = """
## 7. Tropospheric Delay Correction (optional and recommended)
## correct tropospheric delay using the following methods:
## a. pyaps - use weather re-analysis data (Jolivet et al., 2011, GRL, need to install PyAPS; Dee et al., 2011)
## b. height_correlation - correct stratified tropospheric delay (Doin et al., 2009, J Applied Geop)
## c. base_trop_cor - (not recommend) baseline error and stratified tropo simultaneously (Jo et al., 2010, Geo J)
pysar.troposphericDelay.method       = auto  #[pyaps / height_correlation / base_trop_cor / no], auto for pyaps
pysar.troposphericDelay.weatherModel = auto  #[ECMWF / MERRA / NARR], auto for ECMWF, for pyaps method
pysar.troposphericDelay.polyOrder    = auto  #[1 / 2 / 3], auto for 1, for height_correlation method
pysar.troposphericDelay.looks        = auto  #[1-inf], auto for 8, Number of looks to be applied to interferogram 
"""

DATA_INFO = """
  re-analysis_dataset        coverage   temporal_resolution    spatial_resolution      latency     analysis
------------------------------------------------------------------------------------------------------------
ERA-Interim (by ECMWF)        Global      00/06/12/18 UTC      0.75 deg (~83 km)       2-month      4D-var
MERRA(2) (by NASA Goddard)    Global      00/06/12/18 UTC      0.5*0.625 (~50 km)     2-3 weeks     3D-var

To download MERRA2, you need an Earthdata account, and pre-authorize the "NASA GESDISC DATA ARCHIVE" application, following https://disc.gsfc.nasa.gov/earthdata-login.
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Tropospheric correction using weather models\n' +
                                     '  PyAPS is used to download and calculate the delay for each time-series epoch.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=REFERENCE+'\n'+DATA_INFO+'\n'+EXAMPLE)
    # For data download
    parser.add_argument('-m', '--model', '-s', dest='trop_model', default='ECMWF',
                        choices={'ECMWF', 'MERRA', 'NARR', 'ERA', 'MERRA1'},
                        help='source of the atmospheric data.\nNARR is working for 1979-Jan to 2014-Oct.')
    parser.add_argument('-d', '--date-list', dest='date_list', nargs='*',
                        help='Read the first column of text file as list of date to download data\n' +
                             'in YYYYMMDD or YYMMDD format')
    parser.add_argument('--hour', help='time of data in HH, e.g. 12, 06')
    parser.add_argument('-w', '--dir', '--weather-dir', dest='weather_dir',
                        help='directory to put downloaded weather data, i.e. ./../WEATHER\n' +
                             'use directory of input timeseries_file if not specified.')

    # For delay calculation
    parser.add_argument('--dem', dest='dem_file',
                        help='DEM file, i.e. radar_4rlks.hgt, srtm1.dem')
    parser.add_argument('-i', dest='inc_angle', default='30.0',
                        help='a file containing all incidence angles, or\n'+
                             'a number representing for the whole image.')
    parser.add_argument('--ref-yx', dest='ref_yx', type=int,
                        nargs=2, help='reference pixel in y/x')
    parser.add_argument('--delay', dest='delay_type', default='comb', choices={'comb', 'dry', 'wet'},
                        help='Delay type to calculate, comb contains both wet and dry delays')

    # For delay correction
    parser.add_argument('-f', '--file', dest='timeseries_file',
                        help='timeseries HDF5 file, i.e. timeseries.h5')
    parser.add_argument('-o', dest='outfile',
                        help='Output file name for trospheric corrected timeseries.')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    if all(not i for i in [inps.date_list, inps.timeseries_file, inps.dem_file]):
        parser.print_help()
        sys.exit(1)
    return inps


###############################################################
def check_inputs(inps):
    parser = create_parser()
    atr = dict()
    if inps.timeseries_file:
        atr = readfile.read_attribute(inps.timeseries_file)
    elif inps.dem_file:
        atr = readfile.read_attribute(inps.dem_file)

    # Get Grib Source
    inps.trop_model = ut.standardize_trop_model(inps.trop_model,
                                                standardWeatherModelNames)
    print('weather model: '+inps.trop_model)

    # output file name
    if not inps.outfile:
        fbase = os.path.splitext(inps.timeseries_file)[0]
        inps.outfile = '{}_{}.h5'.format(fbase, inps.trop_model)

    # hour
    if not inps.hour:
        if 'CENTER_LINE_UTC' in atr.keys():
            inps.hour = ptime.closest_weather_product_time(atr['CENTER_LINE_UTC'],
                                                           inps.trop_model)
        else:
            parser.print_usage()
            print('ERROR: no input for hour')
            sys.exit(1)
    print('time of cloest available product: {}:00 UTC'.format(inps.hour))

    # date list
    if inps.timeseries_file:
        print('read date list from timeseries file: {}'.format(inps.timeseries_file))
        ts_obj = timeseries(inps.timeseries_file)
        ts_obj.open(print_msg=False)
        inps.date_list = ts_obj.dateList
    elif len(inps.date_list) == 1:
        if os.path.isfile(inps.date_list[0]):
            print('read date list from text file: {}'.format(inps.date_list[0]))
            inps.date_list = ptime.yyyymmdd(np.loadtxt(inps.date_list[0],
                                                       dtype=bytes,
                                                       usecols=(0,)).astype(str).tolist())
        else:
            parser.print_usage()
            raise Exception('ERROR: input date list < 2')

    # weather directory
    if not inps.weather_dir:
        if inps.timeseries_file:
            inps.weather_dir = os.path.join(os.path.dirname(os.path.abspath(inps.timeseries_file)),
                                            '../WEATHER')
        elif inps.dem_file:
            inps.weather_dir = os.path.join(os.path.dirname(os.path.abspath(inps.dem_file)),
                                            '../WEATHER')
        else:
            inps.weather_dir = os.path.abspath(os.getcwd())
    print('weather data directory: '+inps.weather_dir)

    if 'REF_Y' in atr.keys():
        inps.ref_yx = [int(atr['REF_Y']), int(atr['REF_X'])]
        print('reference pixel: {}'.format(inps.ref_yx))

    # Incidence angle: to map the zenith delay to the slant delay
    if os.path.isfile(inps.inc_angle):
        print('incidence angle from file: {}'.format(inps.inc_angle))
        inps.inc_angle = readfile.read(inps.inc_angle,
                                       datasetName='incidenceAngle',
                                       print_msg=False)[0]
    else:
        print('incidence angle from input: {}'.format(inps.inc_angle))
        inps.inc_angle = float(inps.inc_angle)
    inps.inc_angle = inps.inc_angle*np.pi/180.0

    # Coordinate system: geocoded or not
    inps.geocoded = False
    if 'Y_FIRST' in atr.keys():
        inps.geocoded = True
    print('geocoded: {}'.format(inps.geocoded))

    # Prepare DEM file in ROI_PAC format for PyAPS to read
    if inps.dem_file:
        inps.dem_file = prepare_roipac_dem(inps.dem_file, inps.geocoded)

    return inps, atr


###############################################################
def date_list2grib_file(date_list, hour, trop_model, grib_dir):
    grib_file_list = []
    for d in date_list:
        grib_file = grib_dir+'/'
        if   trop_model == 'ECMWF' :  grib_file += 'ERA-Int_%s_%s.grb' % (d, hour)
        elif trop_model == 'MERRA' :  grib_file += 'merra-%s-%s.nc4' % (d, hour)
        elif trop_model == 'NARR'  :  grib_file += 'narr-a_221_%s_%s00_000.grb' % (d, hour)
        elif trop_model == 'ERA'   :  grib_file += 'ERA_%s_%s.grb' % (d, hour)
        elif trop_model == 'MERRA1':  grib_file += 'merra-%s-%s.hdf' % (d, hour)
        grib_file_list.append(grib_file)
    return grib_file_list


def dload_grib_pyaps(date_list, hour, trop_model='ECMWF', weather_dir='./'):
    """Download weather re-analysis grib files using PyAPS
    Inputs:
        date_list   : list of string in YYYYMMDD format
        hour        : string in HH:MM or HH format
        trop_model : string, 
        weather_dir : string,
    Output:
        grib_file_list : list of string
    """
    print('*'*50+'\nDownloading weather model data using PyAPS (Jolivet et al., 2011, GRL) ...')
    # Grib data directory
    grib_dir = weather_dir+'/'+trop_model
    if not os.path.isdir(grib_dir):
        os.makedirs(grib_dir)
        print('making directory: '+grib_dir)

    # Date list to grib file list
    grib_file_list = date_list2grib_file(date_list, hour, trop_model, grib_dir)

    # Get date list to download (skip already downloaded files)
    grib_file_existed = ut.get_file_list(grib_file_list)
    if grib_file_existed:
        grib_filesize_digit = ut.most_common([len(str(os.path.getsize(i))) for i in grib_file_existed])
        grib_filesize_max2 = ut.most_common([str(os.path.getsize(i))[0:2] for i in grib_file_existed])
        grib_file_corrupted = [i for i in grib_file_existed
                               if (len(str(os.path.getsize(i))) != grib_filesize_digit
                                   or str(os.path.getsize(i))[0:2] != grib_filesize_max2)]
        print('file size mode: %se%d bytes' % (grib_filesize_max2, grib_filesize_digit-2))
        print('number of grib files existed    : %d' % len(grib_file_existed))
        if grib_file_corrupted:
            print('------------------------------------------------------------------------------')
            print('corrupted grib files detected! Delete them and re-download...')
            print('number of grib files corrupted  : %d' % len(grib_file_corrupted))
            for i in grib_file_corrupted:
                rmCmd = 'rm '+i
                print(rmCmd)
                os.system(rmCmd)
                grib_file_existed.remove(i)
            print('------------------------------------------------------------------------------')
    grib_file2download = sorted(list(set(grib_file_list) - set(grib_file_existed)))
    date_list2download = [str(re.findall('\d{8}', i)[0]) for i in grib_file2download]
    print('number of grib files to download: %d' % len(date_list2download))
    print('------------------------------------------------------------------------------\n')

    # Download grib file using PyAPS
    if   trop_model == 'ECMWF' :  pa.ECMWFdload( date_list2download, hour, grib_dir)
    elif trop_model == 'MERRA' :  pa.MERRAdload( date_list2download, hour, grib_dir)
    elif trop_model == 'NARR'  :  pa.NARRdload(  date_list2download, hour, grib_dir)
    elif trop_model == 'ERA'   :  pa.ERAdload(   date_list2download, hour, grib_dir)
    elif trop_model == 'MERRA1':  pa.MERRA1dload(date_list2download, hour, grib_dir)
    return grib_file_list


def prepare_roipac_dem(demFile, geocoded=False):
    print('convert input DEM to ROIPAC format')
    dem, atr = readfile.read(demFile, datasetName='height')
    if geocoded:
        ext = '.dem'
    else:
        ext = '.hgt'
    demFileOut = '{}4pyaps{}'.format(os.path.splitext(demFile)[0], ext)
    demFileOut = writefile.write(dem, out_file=demFileOut, metadata=atr)
    return demFileOut


def get_delay(grib_file, inps):
    """Get delay matrix using PyAPS for one acquisition
    Inputs:
        grib_file - strng, grib file path
        atr       - dict, including the following attributes:
                    dem_file    - string, DEM file path
                    trop_model - string, Weather re-analysis data source
                    delay_type  - string, comb/dry/wet
                    ref_y/x     - string, reference pixel row/col number
                    inc_angle   - np.array, 0/1/2 D
    Output:
        phs - 2D np.array, absolute tropospheric phase delay relative to ref_y/x
    """
    if inps.geocoded:
        aps = pa.PyAPS_geo(grib_file, inps.dem_file, grib=inps.trop_model,
                           verb=True, Del=inps.delay_type)
    else:
        aps = pa.PyAPS_rdr(grib_file, inps.dem_file, grib=inps.trop_model,
                           verb=True, Del=inps.delay_type)
    phs = np.zeros((aps.ny, aps.nx), dtype=np.float32)
    aps.getdelay(phs, inc=0.0)

    # Get relative phase delay in space
    phs -= phs[inps.ref_yx[0], inps.ref_yx[1]]
    phs /= np.cos(inps.inc_angle)  # Project into LOS direction
    phs *= -1    # reverse the sign for consistency between different phase correction steps/methods
    return phs


def get_delay_timeseries(inps, atr):
    """Calculate delay time-series and write it to HDF5 file."""
    print('*'*50+'\nCalcualting delay for each epoch using PyAPS ...')
    if any(i is None for i in [inps.dem_file, inps.inc_angle, inps.ref_yx]):
        print('No DEM / incidenceAngle / ref_yx found, exit.')
        return

    length = int(atr['LENGTH'])
    width = int(atr['WIDTH'])
    date_num = len(inps.date_list)
    drop_data = np.zeros((date_num, length, width), np.float32)
    for i in range(date_num):
        grib_file = inps.grib_file_list[i]
        date = inps.date_list[i]
        print('calculate phase delay on %s from file %s' %
              (date, os.path.basename(grib_file)))
        drop_data[i] = get_delay(grib_file, inps)

    # Convert relative phase delay on reference date
    try:
        inps.ref_date = atr['REF_DATE']
    except:
        inps.ref_date = inps.date_list[0]
    print('convert to relative phase delay with reference date: '+inps.ref_date)
    inps.ref_idx = inps.date_list.index(inps.ref_date)
    drop_data -= np.tile(drop_data[inps.ref_idx, :, :], (date_num, 1, 1))

    # Write tropospheric delay to HDF5
    tropFile = os.path.join(os.path.dirname(inps.dem_file), inps.trop_model+'.h5')
    ts_obj = timeseries(tropFile)
    ts_obj.write2hdf5(data=drop_data,
                     dates=inps.date_list,
                     metadata=atr,
                     refFile=inps.timeseries_file)

    # Delete temporary DEM file in ROI_PAC format
    if '4pyaps' in inps.dem_file:
        rmCmd = 'rm {f} {f}.rsc'.format(f=inps.dem_file)
        print(rmCmd)
        os.system(rmCmd)

    return drop_data


def correct_delay(timeseries_file, trop_data, out_file=None):
    print('*'*50+'\nCorrecting delay for input time-series')
    ts_obj = timeseries(timeseries_file)
    ts_data = ts_obj.read()
    mask = ts_data == 0.
    ts_data -= trop_data
    ts_data[mask] = 0.

    ts_obj = timeseries(out_file)
    ts_obj.write2hdf5(data=ts_data, refFile=timeseries_file)
    return out_file


###############################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    inps, atr = check_inputs(inps)

    inps.grib_file_list = dload_grib_pyaps(inps.date_list,
                                           inps.hour,
                                           inps.trop_model,
                                           inps.weather_dir)

    drop_data = get_delay_timeseries(inps, atr)

    if atr['FILE_TYPE'] == 'timeseries':
        inps.outfile = correct_delay(inps.timeseries_file, drop_data, inps.outfile)

    return inps.outfile


###############################################################
if __name__ == '__main__':
    main()
