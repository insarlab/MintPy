#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2015-2019, Zhang Yunjun, Heresh Fattahi     #
# Author:  Zhang Yunjun, Heresh Fattahi                    #
############################################################


import os
import re
import subprocess
import argparse
import numpy as np
from pysar.objects import timeseries, geometry
from pysar.utils import readfile, writefile, ptime, utils as ut

try:
    import pyaps3 as pa
except ImportError:
    raise ImportError('Cannot import pyaps3!')


weatherModelHours = {
    'ERA5'   : [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23],
    'ERAINT' : [0, 6, 12, 18],
    'MERRA'  : [0, 6, 12, 18],
}

verbose = False


###############################################################
EXAMPLE = """example:
  # download reanalysys dataset, calculate tropospheric delays and correct time-series file.
  tropo_pyaps3.py -f timeseries.h5 -g INPUTS/geometryRadar.h5

  # download reanalysys dataset, calculate tropospheric delays
  tropo_pyaps3.py -d date_list.txt     --hour 12 -m ERA5  -g INPUTS/geometryRadar.h5 --ref-yx 30 40
  tropo_pyaps3.py -d 20151002 20151003 --hour 12 -m MERRA -g INPUTS/geometryRadar.h5 --ref-yx 30 40

  # download reanalysys dataset
  tropo_pyaps3.py -d date_list.txt     --hour 12 -m ECMWF
"""

REFERENCE = """reference:
  Jolivet, R., R. Grandin, C. Lasserre, M.-P. Doin and G. Peltzer (2011), Systematic InSAR tropospheric
  phase delay corrections from global meteorological reanalysis data, Geophys. Res. Lett., 38, L17311,
  doi:10.1029/2011GL048757

  Jolivet, R., P. S. Agram, N. Y. Lin, M. Simons, M. P. Doin, G. Peltzer, and Z. Li (2014), Improving
  InSAR geodesy using global atmospheric models, Journal of Geophysical Research: Solid Earth, 119(3),
  2324-2341.
"""

TEMPLATE = """
## 7. Tropospheric Delay Correction (optional and recommended)
## For pyaps method, correction is applied to dates with data available, and skipped for dates (usually recent) without it.
pysar.troposphericDelay.method       = auto  #[pyaps / height_correlation / base_trop_cor / no], auto for pyaps
pysar.troposphericDelay.weatherModel = auto  #[ERA5 / ERAINT / MERRA / NARR], auto for ERA5, for pyaps method
pysar.troposphericDelay.weatherDir   = auto  #[path2directory], auto for "./../WEATHER"
"""

DATA_INFO = """
  re-analysis_dataset       coverage   temporal_resolution  spatial_resolution      latency     analysis
------------------------------------------------------------------------------------------------------------
ERA-5    (by ECMWF)          Global      Hourly              0.25 deg (~31 km)       3-month      4D-var
ERA-Int  (by ECMWF)          Global      6-hourly            0.75 deg (~79 km)       2-month      4D-var
MERRA(2) (by NASA Goddard)   Global      6-hourly           0.5*0.625 (~50 km)      2-3 weeks     3D-var

To download MERRA2, you need an Earthdata account, and pre-authorize the "NASA GESDISC DATA ARCHIVE" application, following https://disc.gsfc.nasa.gov/earthdata-login.
"""

WEATHER_DIR_DEMO = """--weather-dir ~/WEATHER
WEATHER/
    /ECMWF
        ERA-Int_20030329_06.grb
        ERA-Int_20030503_06.grb
    /MERRA
        merra-20110126-06.nc4
        merra-20110313-06.nc4
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Tropospheric correction using weather models\n' +
                                     '  PyAPS is used to download and calculate the delay for each time-series epoch.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=REFERENCE+'\n'+DATA_INFO+'\n'+EXAMPLE)
    # For data download
    parser.add_argument('-m', '--model', '-s', dest='trop_model', default='ERA5',
                        choices={'ERA5', 'ERAINT', 'MERRA', 'NARR'},
                        help='source of the atmospheric data.\nNARR is working for 1979-Jan to 2014-Oct.')
    parser.add_argument('-d', '--date-list', dest='date_list', nargs='*',
                        help='Read the first column of text file as list of date to download data\n' +
                             'in YYYYMMDD or YYMMDD format')
    parser.add_argument('--hour', help='time of data in HH, e.g. 12, 06')
    parser.add_argument('-w', '--dir', '--weather-dir', dest='weather_dir', default='${WEATHER_DIR}',
                        help='parent directory of downloaded weather data file. Default: ${WEATHER_DIR}\n' +
                             'e.g.: '+WEATHER_DIR_DEMO)

    # For delay calculation
    parser.add_argument('-g','--geomtry', dest='geom_file', type=str,
                        help='geometry file including height, incidenceAngle and/or latitude and longitude')
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

    # check the input requirements
    key_list = ['date_list', 'hour']
    # with timeseries file
    if inps.timeseries_file:
        for key in key_list+['ref_yx']:
            if vars(inps)[key]:
                print(('input "{:<10}" is ignored because it will be extracted from '
                       'timeseries file {}').format(key, inps.timeseries_file))

    # without timeseries file
    elif any(not vars(inps)[key] for key in key_list):
        msg = 'ERROR: No input timeseries file, all the following options are required: \n{}'.format(key_list)
        msg += '\n\n'+EXAMPLE
        raise SystemExit(msg)

    ## default values
    print('weather model: '+inps.trop_model)

    # weather_dir
    inps.weather_dir = os.path.expanduser(inps.weather_dir)
    inps.weather_dir = os.path.expandvars(inps.weather_dir)
    # Fallback value if WEATHER_DIR is not defined as environmental variable
    if inps.weather_dir == '${WEATHER_DIR}':
        inps.weather_dir = './'
    print('weather data directory: '+inps.weather_dir)

    return inps


###############################################################
def check_inputs(inps):
    parser = create_parser()

    # output directories/files
    atr = dict()
    pysar_dir = None
    if inps.timeseries_file:
        atr = readfile.read_attribute(inps.timeseries_file)
        pysar_dir = os.path.dirname(inps.timeseries_file)
        if not inps.outfile:
            fbase = os.path.splitext(inps.timeseries_file)[0]
            inps.outfile = '{}_{}.h5'.format(fbase, inps.trop_model)
    elif inps.geom_file:
        atr = readfile.read_attribute(inps.geom_file)
        pysar_dir = os.path.join(os.path.dirname(inps.geom_file), '..')
    else:
        pysar_dir = os.path.abspath(os.getcwd())

    # trop_file
    inps.trop_file = os.path.join(pysar_dir, 'INPUTS/{}.h5'.format(inps.trop_model))
    print('output tropospheric delay file: {}'.format(inps.trop_file))

    # hour
    if not inps.hour:
        if 'CENTER_LINE_UTC' in atr.keys():
            inps.hour = closest_weather_model_hour(atr['CENTER_LINE_UTC'], inps.trop_model)
        else:
            parser.print_usage()
            raise Exception('no input for hour')
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

    # Grib data directory
    inps.grib_dir = os.path.join(inps.weather_dir, inps.trop_model)
    if not os.path.isdir(inps.grib_dir):
        os.makedirs(inps.grib_dir)
        print('making directory: '+inps.grib_dir)

    # area extent for ERA5 grib data download
    inps.snwe = get_snwe(atr)

    # Date list to grib file list
    inps.grib_file_list = date_list2grib_file(date_list=inps.date_list,
                                              hour=inps.hour,
                                              model=inps.trop_model,
                                              grib_dir=inps.grib_dir,
                                              snwe=inps.snwe)

    if 'REF_Y' in atr.keys():
        inps.ref_yx = [int(atr['REF_Y']), int(atr['REF_X'])]
        print('reference pixel: {}'.format(inps.ref_yx))
    return inps, atr


###############################################################
def closest_weather_model_hour(sar_acquisition_time, grib_source='ERA5'):
    """Find closest available time of weather product from SAR acquisition time
    Inputs:
        sar_acquisition_time - string, SAR data acquisition time in seconds
        grib_source - string, Grib Source of weather reanalysis product
    Output:
        grib_hr - string, time of closest available weather product 
    Example:
        '06' = closest_weather_model_hour(atr['CENTER_LINE_UTC'])
        '12' = closest_weather_model_hour(atr['CENTER_LINE_UTC'], 'NARR')
    """
    # Get hour/min of SAR acquisition time
    sar_time = float(sar_acquisition_time)

    # Find closest time in available weather products
    grib_hr_list = weatherModelHours[grib_source]
    grib_hr = int(min(grib_hr_list, key=lambda x: abs(x-sar_time/3600.)))

    # Adjust time output format
    grib_hr = "%02d" % grib_hr
    return grib_hr


def date_list2grib_file(date_list, hour, model, grib_dir, snwe=None):
    # area extent
    area = snwe2str(snwe)

    grib_files = []
    for d in date_list:
        if model == 'ERA5':
            if area:
                grib_file = 'ERA-5{}_{}_{}.grb'.format(area, d, hour)
            else:
                grib_file = 'ERA-5_{}_{}.grb'.format(d, hour)

        elif model == 'ERAINT': grib_file = 'ERA-Int_{}_{}.grb'.format(d, hour)
        elif model == 'MERRA' : grib_file = 'merra-{}-{}.nc4'.format(d, hour)
        elif model == 'NARR'  : grib_file = 'narr-a_221_{}_{}00_000.grb'.format(d, hour)
        elif model == 'ERA'   : grib_file = 'ERA_{}_{}.grb'.format(d, hour)
        elif model == 'MERRA1': grib_file = 'merra-{}-{}.hdf'.format(d, hour)
        grib_files.append(os.path.join(grib_dir, grib_file))
    return grib_files


def ceil_to_5(x):
    """Return the closest number in multiple of 5 in the larger direction"""
    assert isinstance(x, (int, np.int16, np.int32, np.int64)), 'input number is not int: {}'.format(type(x))
    if x % 5 == 0:
        return x
    return x + (5 - x % 5)


def floor_to_5(x):
    """Return the closest number in multiple of 5 in the lesser direction"""
    assert isinstance(x, (int, np.int16, np.int32, np.int64)), 'input number is not int: {}'.format(type(x))
    return x - x % 5


def get_snwe(meta, min_buffer=2, multi_5=True):
    # get bounding box
    lat0, lat1, lon0, lon1 = get_bounding_box(meta)

    # lat/lon0/1 --> SNWE
    S = np.floor(min(lat0, lat1) - min_buffer).astype(int)
    N = np.ceil( max(lat0, lat1) + min_buffer).astype(int)
    W = np.floor(min(lon0, lon1) - min_buffer).astype(int)
    E = np.ceil( max(lon0, lon1) + min_buffer).astype(int)

    # SNWE in multiple of 5
    if multi_5:
        S = floor_to_5(S)
        W = floor_to_5(W)
        N = ceil_to_5(N)
        E = ceil_to_5(E)
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


def dload_grib_files(grib_file_list, trop_model='ERA5', snwe=None):
    """Download weather re-analysis grib files using PyAPS
    Parameters: grib_file_list : list of string of grib files
    Returns:    grib_file_list : list of string
    """
    print('\n------------------------------------------------------------------------------')
    print('downloading weather model data using PyAPS ...')

    # Get date list to download (skip already downloaded files)
    grib_file_exist = check_exist_grib_file(grib_file_list, print_msg=True)
    grib_file2dload = sorted(list(set(grib_file_list) - set(grib_file_exist)))
    date_list2dload = [str(re.findall('\d{8}', i)[0]) for i in grib_file2dload]
    print('number of grib files to download: %d' % len(date_list2dload))
    print('------------------------------------------------------------------------------\n')

    # Download grib file using PyAPS
    if len(date_list2dload) > 0:
        hour = re.findall('\d{8}[-_]\d{2}', grib_file2dload[0])[0].replace('-', '_').split('_')[1]
        grib_dir = os.path.dirname(grib_file2dload[0])

        # try 3 times to download, then use whatever downloaded to calculate delay
        i = 0
        while i < 3:
            i += 1
            try:
                if trop_model in ['ERA5', 'ERAINT']:
                    pa.ECMWFdload(date_list2dload, hour, grib_dir,
                                  model=trop_model,
                                  snwe=snwe,
                                  flist=grib_file2dload)

                elif trop_model == 'MERRA':
                    pa.MERRAdload(date_list2dload, hour, grib_dir)

                elif trop_model == 'NARR':
                    pa.NARRdload(date_list2dload, hour, grib_dir)
            except:
                pass

    grib_file_list = check_exist_grib_file(grib_file_list, print_msg=False)
    return grib_file_list


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
    # initiate pyaps object
    aps_obj = pa.PyAPS(grib_file,
                       grib=inps.trop_model,
                       Del=inps.delay_type,
                       dem=inps.dem,
                       inc=inps.inc,
                       lat=inps.lat,
                       lon=inps.lon,
                       verb=verbose)

    # estimate delay
    phs = np.zeros((aps_obj.ny, aps_obj.nx), dtype=np.float32)
    aps_obj.getdelay(phs)

    # Get relative phase delay in space
    phs -= phs[inps.ref_yx[0], inps.ref_yx[1]]
    phs *= -1    # reverse the sign for consistency between different phase correction steps/methods
    return phs

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


def get_lat_lon(meta):
    """Get 2D array of lat and lon from metadata"""
    length, width = int(meta['LENGTH']), int(meta['WIDTH'])
    lat0, lat1, lon0, lon1 = get_bounding_box(meta)
    lat, lon = np.mgrid[lat0:lat1:length*1j, lon0:lon1:width*1j]
    return lat, lon


def get_delay_timeseries(inps, atr):
    """Calculate delay time-series and write it to HDF5 file.
    Parameters: inps : namespace, all input parameters
                atr  : dict, metadata to be saved in trop_file
    Returns:    trop_file : str, file name of ECMWF.h5
    """
    def get_dataset_size(fname):
        atr = readfile.read_attribute(fname)
        return (atr['LENGTH'], atr['WIDTH'])

    # check 1 - existing tropo delay file
    if (ut.run_or_skip(out_file=inps.trop_file, in_file=inps.grib_file_list, print_msg=False) == 'skip' 
            and get_dataset_size(inps.trop_file) == get_dataset_size(inps.geom_file)):
        print('{} file exists and is newer than all GRIB files, skip updating.'.format(inps.trop_file))
        return

    # check 2 - geometry file
    if any(i is None for i in [inps.geom_file, inps.ref_yx]):
        print('No DEM / incidenceAngle / ref_yx found, skip calculating tropospheric delays.')
        if not os.path.isfile(inps.trop_file):
            inps.trop_file = None
        return

    # prepare geometry data
    geom_obj = geometry(inps.geom_file)
    geom_obj.open()
    inps.dem = geom_obj.read(datasetName='height')
    inps.inc = geom_obj.read(datasetName='incidenceAngle')
    if 'latitude' in geom_obj.datasetNames:
        inps.lat = geom_obj.read(datasetName='latitude')
        inps.lon = geom_obj.read(datasetName='longitude')
    else:
        inps.lat, inps.lon = get_lat_lon(geom_obj.metadata)

    # calculate phase delay
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    num_date = len(inps.grib_file_list)
    date_list = [str(re.findall('\d{8}', i)[0]) for i in inps.grib_file_list]
    trop_data = np.zeros((num_date, length, width), np.float32)

    print('calcualting delay for each date using PyAPS (Jolivet et al., 2011; 2014) ...')
    print('number of grib files used: {}'.format(num_date))
    if verbose:
        prog_bar = ptime.progressBar(maxValue=num_date)
    for i in range(num_date):
        grib_file = inps.grib_file_list[i]
        trop_data[i] = get_delay(grib_file, inps)
        if verbose:
            prog_bar.update(i+1, suffix=os.path.basename(grib_file))
    if verbose:
        prog_bar.close()

    # Convert relative phase delay on reference date
    inps.ref_date = atr.get('REF_DATE', date_list[0])
    print('convert to relative phase delay with reference date: '+inps.ref_date)
    inps.ref_idx = date_list.index(inps.ref_date)
    trop_data -= np.tile(trop_data[inps.ref_idx, :, :], (num_date, 1, 1))

    # Write tropospheric delay to HDF5
    atr['REF_Y'] = inps.ref_yx[0]
    atr['REF_X'] = inps.ref_yx[1]
    ts_obj = timeseries(inps.trop_file)
    ts_obj.write2hdf5(data=trop_data,
                      dates=date_list,
                      metadata=atr,
                      refFile=inps.timeseries_file)
    return


def correct_timeseries(timeseries_file, trop_file, out_file):
    print('\n------------------------------------------------------------------------------')
    print('correcting delay for input time-series by calling diff.py')
    cmd = 'diff.py {} {} -o {} --force'.format(timeseries_file,
                                               trop_file,
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
    inps, atr = check_inputs(inps)

    inps.grib_file_list = dload_grib_files(inps.grib_file_list,
                                           inps.trop_model,
                                           snwe=inps.snwe)

    get_delay_timeseries(inps, atr)

    if atr and atr['FILE_TYPE'] == 'timeseries':
        inps.outfile = correct_timeseries(inps.timeseries_file,
                                          inps.trop_file,
                                          out_file=inps.outfile)
    else:
        print('No input timeseries file, skip correcting tropospheric delays.')

    return inps.outfile


###############################################################
if __name__ == '__main__':
    main()
