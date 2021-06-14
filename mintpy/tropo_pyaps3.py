#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2015               #
############################################################


import os
import sys
import re
import subprocess
import argparse
import h5py
import numpy as np
from mintpy.objects import timeseries, geometry
from mintpy.utils import ptime, readfile, writefile, utils as ut

try:
    import pyaps3 as pa
except ImportError:
    raise ImportError('Cannot import pyaps3!')


WEATHER_MODEL_HOURS = {
    'ERA5'   : [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23],
    'ERAINT' : [0, 6, 12, 18],
    'MERRA'  : [0, 6, 12, 18],
}


###############################################################
EXAMPLE = """example:
  # download datasets, calculate tropospheric delays and correct time-series file.
  tropo_pyaps3.py -f timeseries.h5 -g inputs/geometryRadar.h5
  tropo_pyaps3.py -f filt_fine.unw -g ../../../mintpy/inputs/geometryRadar.h5

  # download datasets, calculate tropospheric delays
  tropo_pyaps3.py -d date.list         --hour 12 -m ERA5  -g inputs/geometryGeo.h5
  tropo_pyaps3.py -d 20151002 20151003 --hour 12 -m MERRA -g inputs/geometryRadar.h5

  # download datasets (covering the whole world)
  tropo_pyaps3.py -d date.list --hour 12
  tropo_pyaps3.py -d SAFE_files.txt
  # download datasets (covering the area of interest)
  tropo_pyaps3.py -d SAFE_files.txt -g inputs/geometryRadar.h5
"""

SAFE_FILE = """SAFE_files.txt:
    /data/SanAndreasSenDT42/SLC/S1B_IW_SLC__1SDV_20191117T140737_20191117T140804_018968_023C8C_82DC.zip
    /data/SanAndreasSenDT42/SLC/S1A_IW_SLC__1SDV_20191111T140819_20191111T140846_029864_036803_69CA.zip
    ...
"""

REFERENCE = """reference:
  Jolivet, R., R. Grandin, C. Lasserre, M.-P. Doin and G. Peltzer (2011), Systematic InSAR tropospheric
  phase delay corrections from global meteorological reanalysis data, Geophys. Res. Lett., 38, L17311,
  doi:10.1029/2011GL048757

  Jolivet, R., P. S. Agram, N. Y. Lin, M. Simons, M. P. Doin, G. Peltzer, and Z. Li (2014), Improving
  InSAR geodesy using global atmospheric models, Journal of Geophysical Research: Solid Earth, 119(3),
  2324-2341, doi:10.1002/2013JB010588.

  # ERA-5
  Hersbach, H., Bell, B., Berrisford, P., Hirahara, S., Horányi, A., Muñoz-Sabater, J., et al. (2020). 
  The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society, 146(730), 1999–2049.
  https://doi.org/10.1002/qj.3803
"""

DATA_INFO = """Global Atmospheric Models:
  re-analysis_dataset      coverage  temp_resolution  spatial_resolution       latency       assimilation
  --------------------------------------------------------------------------------------------------------
  ERA-5(T) (ECMWF)          global       hourly        0.25 deg (~31 km)   3 months (5 days)    4D-Var
  ERA-Int  (ECMWF)          global       6-hourly      0.75 deg (~79 km)        2 months        4D-Var
  MERRA(2) (NASA Goddard)   global       6-hourly     0.5*0.625 (~50 km)       2-3 weeks        3D-Var
  NARR     (NOAA, working from Jan 1979 to Oct 2014)

Notes for data access:
  For MERRA2, you need an Earthdata account, and pre-authorize the "NASA GESDISC DATA ARCHIVE" application
      following https://disc.gsfc.nasa.gov/earthdata-login.
  For ERA-5 from CDS, you need to agree to the Terms of Use of every datasets that you intend to download.
"""

WEATHER_DIR_DEMO = """--weather-dir ~/data/aux
atmosphere/
    /ERA5
        ERA5_N20_N40_E120_E140_20060624_14.grb
        ERA5_N20_N40_E120_E140_20060924_14.grb
        ...
    /MERRA
        merra-20110126-06.nc4
        merra-20110313-06.nc4
        ...
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Tropospheric correction using weather models\n' +
                                     '  PyAPS is used to download and calculate the delay for each acquisition.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=REFERENCE+'\n'+DATA_INFO+'\n'+EXAMPLE)

    parser.add_argument('-f', '--file', dest='dis_file',
                        help='timeseries HDF5 file, i.e. timeseries.h5')
    parser.add_argument('-d', '--date-list', dest='date_list', type=str, nargs='*',
                        help='List of dates in YYYYMMDD or YYMMDD format. It can be:\n'
                             'a) list of strings in YYYYMMDD or YYMMDD format OR\n'
                             'b) a text file with the first column as list of date in YYYYMMDD or YYMMDD format OR\n'
                             'c) a text file with Sentinel-1 SAFE filenames\ne.g.: '+SAFE_FILE)
    parser.add_argument('--hour', type=str, help='time of data in HH, e.g. 12, 06')
    parser.add_argument('-o', dest='cor_dis_file',
                        help='Output file name for trospheric corrected timeseries.')

    # delay calculation
    delay = parser.add_argument_group('delay calculation')
    delay.add_argument('-m', '--model', '-s', dest='tropo_model', default='ERA5',
                       choices={'ERA5'},
                       #choices={'ERA5', 'MERRA', 'NARR'},
                       help='source of the atmospheric model (default: %(default)s).')
    delay.add_argument('--delay', dest='delay_type', default='comb', choices={'comb', 'dry', 'wet'},
                       help='Delay type to calculate, comb contains both wet and dry delays (default: %(default)s).')

    delay.add_argument('-w', '--dir', '--weather-dir', dest='weather_dir', default='${WEATHER_DIR}',
                       help='parent directory of downloaded weather data file (default: %(default)s).\n' +
                            'e.g.: '+WEATHER_DIR_DEMO)
    delay.add_argument('-g','--geomtry', dest='geom_file', type=str,
                       help='geometry file including height, incidenceAngle and/or latitude and longitude')
    delay.add_argument('--custom-height', dest='custom_height', type=float,
                       help='[for testing] specify a custom height value for delay calculation.')

    delay.add_argument('--tropo-file', dest='tropo_file', type=str,
                       help='tropospheric delay file name')
    delay.add_argument('--verbose', dest='verbose', action='store_true', help='Verbose message.')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    ## print model info
    msg = 'weather model: {}'.format(inps.tropo_model)
    if inps.delay_type == 'dry':
        msg += ' - dry (hydrostatic) delay'
    elif inps.delay_type == 'wet':
        msg += ' - wet delay'
    else:
        msg += ' - dry (hydrostatic) and wet delay'
    print(msg)

    ## weather_dir
    # expand path for ~ and environment variables in the path
    inps.weather_dir = os.path.expanduser(inps.weather_dir)
    inps.weather_dir = os.path.expandvars(inps.weather_dir)
    # fallback value if WEATHER_DIR is not defined as environmental variable
    if inps.weather_dir == '${WEATHER_DIR}':
        inps.weather_dir = './'
    inps.weather_dir = os.path.abspath(inps.weather_dir)
    print('weather directory: {}'.format(inps.weather_dir))

    ## ignore invalid filename inputs
    for key in ['dis_file', 'geom_file']:
        fname = vars(inps)[key]
        if fname and not os.path.isfile(fname):
            raise FileNotFoundError('input file not exist: {}'.format(fname))

    ## required options (for date/time): --file OR --date-list
    if (not inps.dis_file 
            and any(vars(inps)[key] is None for key in ['date_list'])):
        raise SystemExit('ERROR: --file OR --date-list is required.\n\n'+EXAMPLE)

    ## output filename - tropo delay file
    if inps.geom_file and not inps.tropo_file:
        inps.tropo_file = os.path.join(os.path.dirname(inps.geom_file), '{}.h5'.format(inps.tropo_model))
    if inps.tropo_file:
        print('output tropospheric delay file: {}'.format(inps.tropo_file))

    ## output filename - corrected displacement file
    if inps.dis_file and not inps.cor_dis_file:
        fbase, fext = os.path.splitext(inps.dis_file)
        inps.cor_dis_file = '{fbase}_{suffix}{fext}'.format(fbase=fbase, suffix=inps.tropo_model, fext=fext)
    if inps.cor_dis_file:
        print('output corrected time-series file: {}'.format(inps.cor_dis_file))

    return inps


###############################################################
def read_inps2date_time(inps):
    """Read dates and time info from input arguments.
    Related options: --file OR --date-list, --hour

    Parameters: inps      - Namespace for the input arguments
    Returns:    date_list - list of str, dates in YYYYMMDD format
                hour      - str, hour in 2-digit with zero padding
    """

    # if --file is specified
    if inps.dis_file:
        # 1) ignore --date-list and --hour
        for key in ['date_list', 'hour']:
            if vars(inps)[key] is not None:
                vars(inps)[key] = None
                msg = 'input "{:<10}" is ignored'.format(key)
                msg += ', use info from file {} instead'.format(inps.dis_file)
                print(msg)

        # 2) read dates/time from time-series file
        print('read dates/time info from file: {}'.format(inps.dis_file))
        atr = readfile.read_attribute(inps.dis_file)
        if atr['FILE_TYPE'] == 'timeseries':
            ts_obj = timeseries(inps.dis_file)
            ts_obj.open(print_msg=False)
            inps.date_list = ts_obj.dateList
        else:
            inps.date_list = ptime.yyyymmdd(atr['DATE12'].split('-'))
        inps.hour = closest_weather_model_hour(atr['CENTER_LINE_UTC'], grib_source=inps.tropo_model)

    # read dates if --date-list is text file
    if len(inps.date_list) == 1 and os.path.isfile(inps.date_list[0]):
        date_file = inps.date_list[0]
        if date_file.startswith('SAFE_'):
            print('read date list and hour info from Sentinel-1 SAFE filenames: {}'.format(date_file))
            inps.date_list, inps.hour = safe2date_time(date_file, inps.tropo_model)
        else:
            print('read date list from text file: {}'.format(date_file))
            inps.date_list = np.loadtxt(date_file, dtype=bytes, usecols=(0,)).astype(str).tolist()
            inps.date_list = ptime.yyyymmdd(inps.date_list)

    # at east 2 dates are required (for meaningful calculation)
    if len(inps.date_list) < 2:
        raise AttributeError('input number of dates < 2!')

    # print time info
    if inps.hour is None:
        raise AttributeError('time info (--hour) not found!')
    print('time of cloest available product: {}:00 UTC'.format(inps.hour))

    return inps.date_list, inps.hour


def get_grib_info(inps):
    """Read the following info from inps
        inps.grib_dir
        inps.atr
        inps.snwe
        inps.grib_files
    """
    # grib data directory, under weather_dir
    inps.grib_dir = os.path.join(inps.weather_dir, inps.tropo_model)
    if not os.path.isdir(inps.grib_dir):
        os.makedirs(inps.grib_dir)
        print('make directory: {}'.format(inps.grib_dir))

    # read metadata
    if inps.dis_file:
        inps.atr = readfile.read_attribute(inps.dis_file)
    elif inps.geom_file:
        inps.atr = readfile.read_attribute(inps.geom_file)
    else:
        inps.atr = dict()

    # area extent for ERA5 grib data download
    if inps.atr:
        inps.snwe = get_snwe(inps.atr, geom_file=inps.geom_file)
    else:
        inps.snwe = None

    # grib file list
    inps.grib_files = get_grib_filenames(date_list=inps.date_list,
                                         hour=inps.hour,
                                         model=inps.tropo_model,
                                         grib_dir=inps.grib_dir,
                                         snwe=inps.snwe)
    return inps


def get_grib_filenames(date_list, hour, model, grib_dir, snwe=None):
    """Get default grib file names based on input info.
    Parameters: date_list  - list of str, date in YYYYMMDD format
                hour       - str, hour in 2-digit with zero padding
                model      - str, global atmospheric model name
                grib_dir   - str, local directory to save grib files
                snwe       - tuple of 4 int, for ERA5 only.
    Returns:    grib_files - list of str, local grib file path
    """
    # area extent
    area = snwe2str(snwe)

    grib_files = []
    for d in date_list:
        if model == 'ERA5':
            if area:
                grib_file = 'ERA5{}_{}_{}.grb'.format(area, d, hour)
            else:
                grib_file = 'ERA5_{}_{}.grb'.format(d, hour)

        elif model == 'ERAINT': grib_file = 'ERA-Int_{}_{}.grb'.format(d, hour)
        elif model == 'MERRA' : grib_file = 'merra-{}-{}.nc4'.format(d, hour)
        elif model == 'NARR'  : grib_file = 'narr-a_221_{}_{}00_000.grb'.format(d, hour)
        elif model == 'ERA'   : grib_file = 'ERA_{}_{}.grb'.format(d, hour)
        elif model == 'MERRA1': grib_file = 'merra-{}-{}.hdf'.format(d, hour)
        grib_files.append(os.path.join(grib_dir, grib_file))
    return grib_files


###############################################################
def closest_weather_model_hour(sar_acquisition_time, grib_source='ERA5'):
    """Find closest available time of weather product from SAR acquisition time
    Inputs:
        sar_acquisition_time - string, SAR data acquisition time in seconds
        grib_source          - string, Grib Source of weather reanalysis product
    Output:
        grib_hr              - string, time of closest available weather product
    Example:
        '06' = closest_weather_model_hour(atr['CENTER_LINE_UTC'])
        '12' = closest_weather_model_hour(atr['CENTER_LINE_UTC'], 'NARR')
    """
    # get hour/min of SAR acquisition time
    sar_time = float(sar_acquisition_time)

    # find closest time in available weather products
    grib_hr_list = WEATHER_MODEL_HOURS[grib_source]
    grib_hr = int(min(grib_hr_list, key=lambda x: abs(x-sar_time/3600.)))

    # add zero padding
    grib_hr = "{:02d}".format(grib_hr)
    return grib_hr


def safe2date_time(safe_file, tropo_model):
    """generate date_list and hour from safe_list"""

    def seconds_UTC(seconds):
        """generate second list"""
        if isinstance(seconds, list):
            secondsOut = []
            for second in seconds:
                secondsOut.append(second)
        else:
            print('\nUn-recognized CENTER_LINE_UTC input!')
            return None

        return secondsOut

    date_list = ptime.yyyymmdd(np.loadtxt(safe_file, dtype=bytes, converters={0:define_date}).astype(str).tolist())
    second_list = seconds_UTC(np.loadtxt(safe_file, dtype=bytes, converters={0:define_second}).astype(str).tolist())

    # change second into hour
    hour_list = [closest_weather_model_hour(float(second), tropo_model) for second in second_list]
    hour = ut.most_common(hour_list)

    return date_list, hour


def define_date(string):
    """extract date from *.SAFE"""
    filename = string.split(str.encode('.'))[0].split(str.encode('/'))[-1]
    date = filename.split(str.encode('_'))[5][0:8]

    return date


def define_second(string):
    """extract CENTER_LINE_UTC from *.SAFE"""
    filename = string.split(str.encode('.'))[0].split(str.encode('/'))[-1]
    time1 = filename.split(str.encode('_'))[5][9:15]
    time2 = filename.split(str.encode('_'))[6][9:15]
    time1_second = int(time1[0:2]) * 3600 + int(time1[2:4]) * 60 + int(time1[4:6])
    time2_second = int(time2[0:2]) * 3600 + int(time2[2:4]) * 60 + int(time2[4:6])
    CENTER_LINE_UTC = (time1_second + time2_second) / 2

    return CENTER_LINE_UTC


def ceil2multiple(x, step=10):
    """Given a number x, find the smallest number in multiple of step >= x."""
    assert isinstance(x, (int, np.int16, np.int32, np.int64)), 'input number is not int: {}'.format(type(x))
    if x % step == 0:
        return x
    return x + (step - x % step)


def floor2multiple(x, step=10):
    """Given a number x, find the largest number in multiple of step <= x."""
    assert isinstance(x, (int, np.int16, np.int32, np.int64)), 'input number is not int: {}'.format(type(x))
    return x - x % step


def get_snwe(meta, geom_file=None, min_buffer=2, step=10):
    # get bounding box
    lat0, lat1, lon0, lon1 = get_bounding_box(meta, geom_file=geom_file)

    # lat/lon0/1 --> SNWE
    S = np.floor(min(lat0, lat1) - min_buffer).astype(int)
    N = np.ceil( max(lat0, lat1) + min_buffer).astype(int)
    W = np.floor(min(lon0, lon1) - min_buffer).astype(int)
    E = np.ceil( max(lon0, lon1) + min_buffer).astype(int)

    # SNWE in multiple of 10
    if step > 1:
        S = floor2multiple(S, step=step)
        W = floor2multiple(W, step=step)
        N = ceil2multiple(N, step=step)
        E = ceil2multiple(E, step=step)
    return (S, N, W, E)


def snwe2str(snwe):
    """Get area extent in string"""
    if not snwe:
        return None

    area = ''
    s, n, w, e = snwe

    if s < 0:
        area += '_S{}'.format(abs(s))
    else:
        area += '_N{}'.format(abs(s))

    if n < 0:
        area += '_S{}'.format(abs(n))
    else:
        area += '_N{}'.format(abs(n))

    if w < 0:
        area += '_W{}'.format(abs(w))
    else:
        area += '_E{}'.format(abs(w))

    if e < 0:
        area += '_W{}'.format(abs(e))
    else:
        area += '_E{}'.format(abs(e))
    return area


def get_bounding_box(meta, geom_file=None):
    """Get lat/lon range (roughly), in the same order of data file
    lat0/lon0 - starting latitude/longitude (first row/column)
    lat1/lon1 - ending latitude/longitude (last row/column)
    """
    length, width = int(meta['LENGTH']), int(meta['WIDTH'])
    if 'Y_FIRST' in meta.keys():
        # geo coordinates
        lat0 = float(meta['Y_FIRST'])
        lon0 = float(meta['X_FIRST'])
        lat_step = float(meta['Y_STEP'])
        lon_step = float(meta['X_STEP'])
        lat1 = lat0 + lat_step * (length - 1)
        lon1 = lon0 + lon_step * (width - 1)

        # 'Y_FIRST' not in 'degree'
        # e.g. meters for UTM projection from ASF HyP3
        if not meta['Y_UNIT'].lower().startswith('deg'):
            lat0, lon0 = ut.to_latlon(meta['OG_FILE_PATH'], lon0, lat0)
            lat1, lon1 = ut.to_latlon(meta['OG_FILE_PATH'], lon1, lat1)

    else:
        # radar coordinates
        if geom_file and os.path.isfile(geom_file):
            geom_dset_list = readfile.get_dataset_list(geom_file)
        else:
            geom_dset_list = []

        if 'latitude' in geom_dset_list:
            lats = readfile.read(geom_file, datasetName='latitude')[0]
            lons = readfile.read(geom_file, datasetName='longitude')[0]
            lats[lats == 0] = np.nan
            lons[lons == 0] = np.nan
            lat0 = np.nanmin(lats)
            lat1 = np.nanmax(lats)
            lon0 = np.nanmin(lons)
            lon1 = np.nanmax(lons)

        else:
            lats = [float(meta['LAT_REF{}'.format(i)]) for i in [1,2,3,4]]
            lons = [float(meta['LON_REF{}'.format(i)]) for i in [1,2,3,4]]
            lat0 = np.mean(lats[0:2])
            lat1 = np.mean(lats[2:4])
            lon0 = np.mean(lons[0:3:2])
            lon1 = np.mean(lons[1:4:2])

    return lat0, lat1, lon0, lon1


###############################################################
def check_exist_grib_file(gfile_list, print_msg=True):
    """Check input list of grib files, and return the existing ones with right size."""
    gfile_exist = ut.get_file_list(gfile_list)
    if gfile_exist:
        file_sizes = [os.path.getsize(i) for i in gfile_exist] # if os.path.getsize(i) > 10e6]
        if file_sizes:
            comm_size = ut.most_common([i for i in file_sizes])
            if print_msg:
                print('common file size: {} bytes'.format(comm_size))
                print('number of grib files existed    : {}'.format(len(gfile_exist)))

            gfile_corrupt = []
            for gfile in gfile_exist:
                if os.path.getsize(gfile) < comm_size * 0.9:
                    gfile_corrupt.append(gfile)
        else:
            gfile_corrupt = gfile_exist

        if gfile_corrupt:
            if print_msg:
                print('------------------------------------------------------------------------------')
                print('corrupted grib files detected! Delete them and re-download...')
                print('number of grib files corrupted  : {}'.format(len(gfile_corrupt)))

            for gfile in gfile_corrupt:
                print('remove {}'.format(gfile))
                os.remove(gfile)
                gfile_exist.remove(gfile)

            if print_msg:
                print('------------------------------------------------------------------------------')
    return gfile_exist


def dload_grib_files(grib_files, tropo_model='ERA5', snwe=None):
    """Download weather re-analysis grib files using PyAPS
    Parameters: grib_files : list of string of grib files
    Returns:    grib_files : list of string
    """
    print('\n------------------------------------------------------------------------------')
    print('downloading weather model data using PyAPS ...')

    # Get date list to download (skip already downloaded files)
    grib_files_exist = check_exist_grib_file(grib_files, print_msg=True)
    grib_files2dload = sorted(list(set(grib_files) - set(grib_files_exist)))
    date_list2dload = [str(re.findall('\d{8}', os.path.basename(i))[0]) for i in grib_files2dload]
    print('number of grib files to download: %d' % len(date_list2dload))
    print('------------------------------------------------------------------------------\n')

    # Download grib file using PyAPS
    if len(date_list2dload) > 0:
        hour = re.findall('\d{8}[-_]\d{2}', os.path.basename(grib_files2dload[0]))[0].replace('-', '_').split('_')[1]
        grib_dir = os.path.dirname(grib_files2dload[0])

        # try 3 times to download, then use whatever downloaded to calculate delay
        i = 0
        while i < 3:
            i += 1
            try:
                if tropo_model in ['ERA5', 'ERAINT']:
                    pa.ECMWFdload(date_list2dload, hour, grib_dir,
                                  model=tropo_model,
                                  snwe=snwe,
                                  flist=grib_files2dload)

                elif tropo_model == 'MERRA':
                    pa.MERRAdload(date_list2dload, hour, grib_dir)

                elif tropo_model == 'NARR':
                    pa.NARRdload(date_list2dload, hour, grib_dir)
            except:
                if i < 3:
                    print('WARNING: the {} attampt to download failed, retry it.\n'.format(i))
                else:
                    print('\n\n'+'*'*50)
                    print('WARNING: downloading failed for 3 times, stop trying and continue.')
                    print('*'*50+'\n\n')
                pass

    # check potentially corrupted files
    grib_files = check_exist_grib_file(grib_files, print_msg=False)
    return grib_files


def get_delay(grib_file, tropo_model, delay_type, dem, inc, lat, lon, mask=None, verbose=False):
    """Get delay matrix using PyAPS for one acquisition
    Parameters: grib_file       - str, grib file path
                tropo_model     - str, GAM model
                delay_type      - str, dry/wet/comb
                dem/inc/lat/lon - 2D np.ndarray in float32 for DEM, incidence angle, latitude/longitude
                verbose         - bool, verbose message
    Returns:    pha             - 2D np.ndarray in float32, single path tropospheric delay
                                  temporally absolute, spatially referenced to ref_y/x
    """
    if verbose:
        print('GRIB FILE: {}'.format(grib_file))

    # initiate pyaps object
    aps_obj = pa.PyAPS(grib_file,
                       grib=tropo_model,
                       Del=delay_type,
                       dem=dem,
                       inc=inc,
                       lat=lat,
                       lon=lon,
                       mask=mask,
                       verb=verbose)

    # estimate delay
    pha = np.zeros((aps_obj.ny, aps_obj.nx), dtype=np.float32)
    aps_obj.getdelay(pha)

    # reverse the sign for consistency between different phase correction steps/methods
    pha *= -1
    return pha


def calc_delay_timeseries(inps):
    """Calculate delay time-series and write it to HDF5 file.
    Parameters: inps : namespace, all input parameters
    Returns:    tropo_file : str, file name of ECMWF.h5
    """
    def get_dataset_size(fname):
        atr = readfile.read_attribute(fname)
        shape = (int(atr['LENGTH']), int(atr['WIDTH']))
        return shape

    def run_or_skip(grib_files, tropo_file, geom_file):
        print('update mode: ON')
        print('output file: {}'.format(tropo_file))
        flag = 'skip'

        # check existance and modification time
        if ut.run_or_skip(out_file=tropo_file, in_file=grib_files, print_msg=False) == 'run':
            flag = 'run'
            print('1) output file either do NOT exist or is NOT newer than all GRIB files.')

        else:
            print('1) output file exists and is newer than all GRIB files.')

            # check dataset size in space / time
            date_list = [str(re.findall('\d{8}', os.path.basename(i))[0]) for i in grib_files]
            if (get_dataset_size(tropo_file) != get_dataset_size(geom_file) 
                    or any(i not in timeseries(tropo_file).get_date_list() for i in date_list)):
                flag = 'run'
                print('2) output file does NOT have the same len/wid as the geometry file {} or does NOT contain all dates'.format(geom_file))
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
        print('run or skip: {}'.format(flag))
        return flag

    if run_or_skip(inps.grib_files, inps.tropo_file, inps.geom_file) == 'skip':
        return


    ## 1. prepare geometry data
    geom_obj = geometry(inps.geom_file)
    geom_obj.open()
    inps.inc = geom_obj.read(datasetName='incidenceAngle')
    inps.dem = geom_obj.read(datasetName='height')

    # for testing
    if inps.custom_height:
        print('use input custom height of {} m for vertical integration'.format(inps.custom_height))
        inps.dem[:] = inps.custom_height

    if 'latitude' in geom_obj.datasetNames:
        # for lookup table in radar-coord (isce, doris)
        inps.lat = geom_obj.read(datasetName='latitude')
        inps.lon = geom_obj.read(datasetName='longitude')

    elif 'Y_FIRST' in geom_obj.metadata:
        # for lookup table in geo-coded (gamma, roipac) and obs. in geo-coord
        inps.lat, inps.lon = ut.get_lat_lon(geom_obj.metadata)

        # convert coordinates to lat/lon, e.g. from UTM for ASF HyPP3
        if not geom_obj.metadata['Y_UNIT'].startswith('deg'):
            inps.lat, inps.lon = ut.to_latlon(inps.atr['OG_FILE_PATH'], inps.lon, inps.lat)

    else:
        # for lookup table in geo-coded (gamma, roipac) and obs. in radar-coord
        inps.lat, inps.lon = ut.get_lat_lon_rdc(inps.atr)

    # mask of valid pixels
    mask = np.multiply(inps.inc != 0, ~np.isnan(inps.inc))


    ## 2. prepare output file
    # metadata
    atr = inps.atr.copy()
    atr['FILE_TYPE'] = 'timeseries'
    atr['UNIT'] = 'm'

    # remove metadata related with double reference
    # because absolute delay is calculated and saved
    for key in ['REF_DATE','REF_X','REF_Y','REF_LAT','REF_LON']:
        if key in atr.keys():
            atr.pop(key)

    # instantiate time-series
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    num_date = len(inps.grib_files)
    date_list = [str(re.findall('\d{8}', os.path.basename(i))[0]) for i in inps.grib_files]
    dates = np.array(date_list, dtype=np.string_)
    ds_name_dict = {
        "date"       : [dates.dtype, (num_date,), dates],
        "timeseries" : [np.float32,  (num_date, length, width), None],
    }
    writefile.layout_hdf5(inps.tropo_file, ds_name_dict, metadata=atr)


    ## 3. calculate phase delay
    print('\n------------------------------------------------------------------------------')
    print('calculating absolute delay for each date using PyAPS (Jolivet et al., 2011; 2014) ...')
    print('number of grib files used: {}'.format(num_date))

    prog_bar = ptime.progressBar(maxValue=num_date, print_msg=~inps.verbose)
    for i in range(num_date):
        grib_file = inps.grib_files[i]

        # calc tropo delay
        tropo_data = get_delay(grib_file,
                               tropo_model=inps.tropo_model,
                               delay_type=inps.delay_type,
                               dem=inps.dem,
                               inc=inps.inc,
                               lat=inps.lat,
                               lon=inps.lon,
                               mask=mask,
                               verbose=inps.verbose)

        # write tropo delay to file
        block = [i, i+1, 0, length, 0, width]
        writefile.write_hdf5_block(inps.tropo_file,
                                   data=tropo_data,
                                   datasetName='timeseries',
                                   block=block,
                                   print_msg=False)

        prog_bar.update(i+1, suffix=os.path.basename(grib_file))
    prog_bar.close()

    return inps.tropo_file


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
    writefile.write(data-tropo, cor_dis_file, atr)
    return cor_dis_file


###############################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    # read dates / time info
    read_inps2date_time(inps)

    # get corresponding grib files info
    get_grib_info(inps)

    # download
    inps.grib_files = dload_grib_files(inps.grib_files, 
                                       tropo_model=inps.tropo_model,
                                       snwe=inps.snwe)

    # calculate tropo delay and save to h5 file
    if inps.geom_file:
        calc_delay_timeseries(inps)
    else:
        print('No input geometry file, skip calculating and correcting tropospheric delays.')
        return

    # correct tropo delay from displacement time-series
    if inps.dis_file:
        ftype = inps.atr['FILE_TYPE']
        if ftype == 'timeseries':
            correct_timeseries(dis_file=inps.dis_file,
                               tropo_file=inps.tropo_file,
                               cor_dis_file=inps.cor_dis_file)

        elif ftype == '.unw':
            correct_single_ifgram(dis_file=inps.dis_file,
                                  tropo_file=inps.tropo_file,
                                  cor_dis_file=inps.cor_dis_file)
        else:
            print('input file {} is not timeseries nor .unw, correction is not supported yet.'.format(ftype))

    else:
        print('No input displacement file, skip correcting tropospheric delays.')

    
    return

###############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
