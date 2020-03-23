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
import numpy as np
from mintpy.objects import timeseries, geometry
from mintpy.utils import readfile, ptime, utils as ut

try:
    import pyaps3 as pa
except ImportError:
    raise ImportError('Cannot import pyaps3!')


WEATHER_MODEL_HOURS = {
    'ERA5'   : [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23],
    'ERAINT' : [0, 6, 12, 18],
    'MERRA'  : [0, 6, 12, 18],
}

verbose = False


###############################################################
EXAMPLE = """example:
  # download reanalysys dataset, calculate tropospheric delays and correct time-series file.
  tropo_pyaps3.py -f timeseries.h5 -g inputs/geometryRadar.h5

  # download reanalysys dataset, calculate tropospheric delays
  tropo_pyaps3.py -d date.list         --hour 12 -m ERA5  -g inputs/geometryRadar.h5
  tropo_pyaps3.py -d 20151002 20151003 --hour 12 -m MERRA -g inputs/geometryRadar.h5

  # download reanalysys dataset
  tropo_pyaps3.py -d date.list --hour 12
"""

REFERENCE = """reference:
  Jolivet, R., R. Grandin, C. Lasserre, M.-P. Doin and G. Peltzer (2011), Systematic InSAR tropospheric
  phase delay corrections from global meteorological reanalysis data, Geophys. Res. Lett., 38, L17311,
  doi:10.1029/2011GL048757

  Jolivet, R., P. S. Agram, N. Y. Lin, M. Simons, M. P. Doin, G. Peltzer, and Z. Li (2014), Improving
  InSAR geodesy using global atmospheric models, Journal of Geophysical Research: Solid Earth, 119(3),
  2324-2341.
"""

DATA_INFO = """Global Atmospheric Models:
  re-analysis_dataset         coverage   temporal_resolution  spatial_resolution      latency     analysis
  --------------------------------------------------------------------------------------------------------
  ERA-5    (ECMWF)             Global      Hourly              0.25 deg (~31 km)       3-month      4D-var
  ERA-Int  (ECMWF)             Global      6-hourly            0.75 deg (~79 km)       2-month      4D-var
  MERRA(2) (NASA Goddard)      Global      6-hourly           0.5*0.625 (~50 km)      2-3 weeks     3D-var
  NARR     (NOAA, working from Jan 1979 to Oct 2014)

Notes for data access:
  For MERRA2, you need an Earthdata account, and pre-authorize the "NASA GESDISC DATA ARCHIVE" application
      following https://disc.gsfc.nasa.gov/earthdata-login.
  For ERA-5 from CDS, you need to agree to the Terms of Use of every datasets that you intend to download.
"""

WEATHER_DIR_DEMO = """--weather-dir ~/atmosphere
atmosphere/
    /ERA5
        ERA5_N20_N40_E120_E140_20060624_14.grb
        ERA5_N20_N40_E120_E140_20060924_14.grb
        ...
    /ECMWF
        ERA-Int_20030329_06.grb
        ERA-Int_20030503_06.grb
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

    parser.add_argument('-f', '--file', dest='timeseries_file',
                        help='timeseries HDF5 file, i.e. timeseries.h5')
    parser.add_argument('-d', '--date-list', dest='date_list', type=str, nargs='*',
                        help='List of dates in YYYYMMDD or YYMMDD format. It can be:\n'
                             'a) list of strings OR\n'
                             'b) a text file with the first column as list of date')
    parser.add_argument('--hour', type=str, help='time of data in HH, e.g. 12, 06')
    parser.add_argument('-o', dest='cor_timeseries_file',
                        help='Output file name for trospheric corrected timeseries.')

    # delay calculation
    delay = parser.add_argument_group('delay calculation')
    delay.add_argument('-m', '--model', '-s', dest='tropo_model', default='ERA5',
                       choices={'ERA5'},
                       #choices={'ERA5', 'ERAINT', 'MERRA', 'NARR'},
                       help='source of the atmospheric model (default: %(default)s).')
    delay.add_argument('--delay', dest='delay_type', default='comb', choices={'comb', 'dry', 'wet'},
                       help='Delay type to calculate, comb contains both wet and dry delays (default: %(default)s).')
    delay.add_argument('-w', '--dir', '--weather-dir', dest='weather_dir', default='${WEATHER_DIR}',
                       help='parent directory of downloaded weather data file (default: %(default)s).\n' +
                            'e.g.: '+WEATHER_DIR_DEMO)
    delay.add_argument('-g','--geomtry', dest='geom_file', type=str,
                       help='geometry file including height, incidenceAngle and/or latitude and longitude')
    delay.add_argument('--tropo-file', dest='tropo_file', type=str,
                       help='tropospheric delay file name')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    ## print model info
    print('weather model: {}'.format(inps.tropo_model))

    ## weather_dir
    # expand path for ~ and environmental variables in the path
    inps.weather_dir = os.path.expanduser(inps.weather_dir)
    inps.weather_dir = os.path.expandvars(inps.weather_dir)
    # fallback value if WEATHER_DIR is not defined as environmental variable
    if inps.weather_dir == '${WEATHER_DIR}':
        inps.weather_dir = './'
    inps.weather_dir = os.path.abspath(inps.weather_dir)
    print('weather directory: {}'.format(inps.weather_dir))

    ## ignore invalid filename inputs
    for key in ['timeseries_file', 'geom_file']:
        if vars(inps)[key] and not os.path.isfile(vars(inps)[key]):
            vars(inps)[key] = None

    ## required options (for date/time): --file OR --date-list, --hour
    if (not inps.timeseries_file 
            and any(vars(inps)[key] is None for key in ['date_list', 'hour'])):
        msg = 'ERROR: --file OR (--date-list, --hour) is required.'
        msg += '\n\n'+EXAMPLE
        raise SystemExit(msg)

    ## output filename - tropo delay file
    if inps.geom_file and not inps.tropo_file:
        inps.tropo_file = os.path.join(os.path.dirname(inps.geom_file), '{}.h5'.format(inps.tropo_model))
    if inps.tropo_file:
        print('output tropospheric delay file: {}'.format(inps.tropo_file))

    ## output filename - corrected displacement file
    if inps.timeseries_file and not inps.cor_timeseries_file:
        fbase = os.path.splitext(inps.timeseries_file)[0]
        inps.cor_timeseries_file = '{}_{}.h5'.format(fbase, inps.tropo_model)
    if inps.cor_timeseries_file:
        print('output corrected time-series file: {}'.format(inps.cor_timeseries_file))

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
    if inps.timeseries_file:
        # 1) ignore --date-list and --hour
        for key in ['date_list', 'hour']:
            if vars(inps)[key] is not None:
                vars(inps)[key] = None
                msg = 'input "{:<10}" is ignored'.format(key)
                msg += ', use info from file {} instead'.format(inps.timeseries_file)
                print(msg)

        # 2) read dates/time from time-series file
        print('read dates/time info from file: {}'.format(inps.timeseries_file))
        ts_obj = timeseries(inps.timeseries_file)
        ts_obj.open(print_msg=False)
        inps.date_list = ts_obj.dateList
        inps.hour = closest_weather_model_hour(ts_obj.metadata['CENTER_LINE_UTC'],
                                               grib_source=inps.tropo_model)

    # read dates if --date-list is text file
    if len(inps.date_list) == 1 and os.path.isfile(inps.date_list[0]):
        print('read date list from text file: {}'.format(inps.date_list[0]))
        inps.date_list = np.loadtxt(inps.date_list[0],
                                    dtype=bytes,
                                    usecols=(0,)).astype(str).tolist()
        inps.date_list = ptime.yyyymmdd(inps.date_list)

    # at east 2 dates are required (for meaningful calculation)
    if len(inps.date_list) < 2:
        raise SystemExit('ERROR: input date list < 2')      

    # print time info
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
    if inps.timeseries_file:
        inps.atr = readfile.read_attribute(inps.timeseries_file)
    elif inps.geom_file:
        inps.atr = readfile.read_attribute(inps.geom_file)
    else:
        inps.atr = dict()

    # area extent for ERA5 grib data download
    if inps.atr:
        inps.snwe = get_snwe(inps.atr)
    else:
        inps.snwe = None

    # grib file list
    inps.grib_files = get_grib_filenames(date_list=inps.date_list,
                                         hour=inps.hour,
                                         model=inps.tropo_model,
                                         grib_dir=inps.grib_dir,
                                         snwe=inps.snwe)
    return inps


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


def ceil2multiple(x, step=10):
    """Return the closest number in multiple of step in the larger direction"""
    assert isinstance(x, (int, np.int16, np.int32, np.int64)), 'input number is not int: {}'.format(type(x))
    if x % step == 0:
        return x
    return x + (step - x % step)


def floor2multiple(x, step=10):
    """Return the closest number in multiple of step in the lesser direction"""
    assert isinstance(x, (int, np.int16, np.int32, np.int64)), 'input number is not int: {}'.format(type(x))
    return x - x % step


def get_snwe(meta, min_buffer=2, step=10):
    # get bounding box
    lat0, lat1, lon0, lon1 = get_bounding_box(meta)

    # lat/lon0/1 --> SNWE
    S = np.floor(min(lat0, lat1) - min_buffer).astype(int)
    N = np.ceil( max(lat0, lat1) + min_buffer).astype(int)
    W = np.floor(min(lon0, lon1) - min_buffer).astype(int)
    E = np.ceil( max(lon0, lon1) + min_buffer).astype(int)

    # SNWE in multiple of 5
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
            for i in gfile_corrupt:
                rmCmd = 'rm '+i
                print(rmCmd)
                os.system(rmCmd)
                gfile_exist.remove(i)
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
    date_list2dload = [str(re.findall('\d{8}', i)[0]) for i in grib_files2dload]
    print('number of grib files to download: %d' % len(date_list2dload))
    print('------------------------------------------------------------------------------\n')

    # Download grib file using PyAPS
    if len(date_list2dload) > 0:
        hour = re.findall('\d{8}[-_]\d{2}', grib_files2dload[0])[0].replace('-', '_').split('_')[1]
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
                pass

    # check potentially corrupted files
    grib_files = check_exist_grib_file(grib_files, print_msg=False)
    return grib_files


def get_delay(grib_file, inps):
    """Get delay matrix using PyAPS for one acquisition
    Inputs:
        grib_file - strng, grib file path
        atr       - dict, including the following attributes:
                    dem_file    - string, DEM file path
                    tropo_model - string, Weather re-analysis data source
                    delay_type  - string, comb/dry/wet
                    ref_y/x     - string, reference pixel row/col number
                    inc_angle   - np.array, 0/1/2 D
    Output:
        pha - 2D np.array, absolute tropospheric phase delay relative to ref_y/x
    """
    # initiate pyaps object
    aps_obj = pa.PyAPS(grib_file,
                       grib=inps.tropo_model,
                       Del=inps.delay_type,
                       dem=inps.dem,
                       inc=inps.inc,
                       lat=inps.lat,
                       lon=inps.lon,
                       verb=verbose)

    # estimate delay
    pha = np.zeros((aps_obj.ny, aps_obj.nx), dtype=np.float32)
    aps_obj.getdelay(pha)

    # reverse the sign for consistency between different phase correction steps/methods
    pha *= -1
    return pha


def get_bounding_box(meta):
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
    else:
        # radar coordinates
        lats = [float(meta['LAT_REF{}'.format(i)]) for i in [1,2,3,4]]
        lons = [float(meta['LON_REF{}'.format(i)]) for i in [1,2,3,4]]
        lat0 = np.mean(lats[0:2])
        lat1 = np.mean(lats[2:4])
        lon0 = np.mean(lons[0:3:2])
        lon1 = np.mean(lons[1:4:2])
    return lat0, lat1, lon0, lon1


def calculate_delay_timeseries(inps):
    """Calculate delay time-series and write it to HDF5 file.
    Parameters: inps : namespace, all input parameters
    Returns:    tropo_file : str, file name of ECMWF.h5
    """
    def get_dataset_size(fname):
        atr = readfile.read_attribute(fname)
        shape = (int(atr['LENGTH']), int(atr['WIDTH']))
        return shape

    # check existing tropo delay file
    if (ut.run_or_skip(out_file=inps.tropo_file, in_file=inps.grib_files, print_msg=False) == 'skip'
            and get_dataset_size(inps.tropo_file) == get_dataset_size(inps.geom_file)):
        print('{} file exists and is newer than all GRIB files, skip updating.'.format(inps.tropo_file))
        return

    # prepare geometry data
    geom_obj = geometry(inps.geom_file)
    geom_obj.open()
    inps.dem = geom_obj.read(datasetName='height')
    inps.inc = geom_obj.read(datasetName='incidenceAngle')

    if 'latitude' in geom_obj.datasetNames:
        # for dataset in geo OR radar coord with lookup table in radar-coord (isce, doris)
        inps.lat = geom_obj.read(datasetName='latitude')
        inps.lon = geom_obj.read(datasetName='longitude')
    elif 'Y_FIRST' in geom_obj.metadata:
        # for geo-coded dataset (gamma, roipac)
        inps.lat, inps.lon = ut.get_lat_lon(geom_obj.metadata)
    else:
        # for radar-coded dataset (gamma, roipac)
        inps.lat, inps.lon = ut.get_lat_lon_rdc(geom_obj.metadata)

    # calculate phase delay
    length, width = int(inps.atr['LENGTH']), int(inps.atr['WIDTH'])
    num_date = len(inps.grib_files)
    date_list = [str(re.findall('\d{8}', i)[0]) for i in inps.grib_files]
    tropo_data = np.zeros((num_date, length, width), np.float32)
    print('\n------------------------------------------------------------------------------')
    print('calcualting absolute delay for each date using PyAPS (Jolivet et al., 2011; 2014) ...')
    print('number of grib files used: {}'.format(num_date))
    prog_bar = ptime.progressBar(maxValue=num_date)
    for i in range(num_date):
        grib_file = inps.grib_files[i]
        tropo_data[i] = get_delay(grib_file, inps)
        prog_bar.update(i+1, suffix=os.path.basename(grib_file))
    prog_bar.close()

    # remove metadata related with double reference
    # because absolute delay is calculated and saved
    for key in ['REF_DATE','REF_X','REF_Y','REF_LAT','REF_LON']:
        if key in inps.atr.keys():
            inps.atr.pop(key)

    # Write tropospheric delay to HDF5
    ts_obj = timeseries(inps.tropo_file)
    ts_obj.write2hdf5(data=tropo_data,
                      dates=date_list,
                      metadata=inps.atr,
                      refFile=inps.timeseries_file)
    return inps.tropo_file


def correct_timeseries(timeseries_file, tropo_file, out_file):
    # diff.py can handle different reference in space and time
    # between the absolute tropospheric delay and the double referenced time-series
    print('\n------------------------------------------------------------------------------')
    print('correcting relative delay for input time-series using diff.py')
    cmd = 'diff.py {} {} -o {} --force'.format(timeseries_file,
                                               tropo_file,
                                               out_file)
    print(cmd)
    status = subprocess.Popen(cmd, shell=True).wait()
    if status is not 0:
        raise Exception(('Error while correcting timeseries file '
                         'using diff.py with tropospheric delay file.'))
    return out_file


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
        calculate_delay_timeseries(inps)
    else:
        print('No input geometry file, skip calculating and correcting tropospheric delays.')
        return

    # correct tropo delay from displacement time-series
    if inps.timeseries_file and inps.atr['FILE_TYPE'] == 'timeseries':
        correct_timeseries(inps.timeseries_file,
                           tropo_file=inps.tropo_file,
                           out_file=inps.cor_timeseries_file)
    else:
        print('No input timeseries file, skip correcting tropospheric delays.')

    return


###############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
