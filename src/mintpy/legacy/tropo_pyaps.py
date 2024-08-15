#!/usr/bin/env python
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, Zhang Yunjun, 2015               #
############################################################


import argparse
import os
import re
import sys

import numpy as np

try:
    import pyaps as pa
except ImportError:
    raise ImportError('Cannot import pyaps!')

import mintpy.cli.diff
from mintpy.objects import geometry, timeseries
from mintpy.utils import ptime, readfile, utils as ut, writefile

standardWeatherModelNames = {
    'ERAI': 'ECMWF', 'ERAINT': 'ECMWF', 'ERAINTERIM': 'ECMWF',
    'MERRA2': 'MERRA',
}


###############################################################
EXAMPLE = """example:
  # download reanalysys dataset, calculate tropospheric delays and correct time-series file.
  tropo_pyaps.py -f timeseries.h5 -m ECMWF -g inputs/geometryRadar.h5 -w ${WEATHER_DIR}

  # download reanalysys dataset, calculate tropospheric delays
  tropo_pyaps.py -d date_list.txt     --hour 12 -m ECMWF -g inputs/geometryRadar.h5 --ref-yx 30 40
  tropo_pyaps.py -d 20151002 20151003 --hour 12 -m MERRA -g inputs/geometryRadar.h5 --ref-yx 30 40

  # download reanalysys dataset
  tropo_pyaps.py -d date_list.txt     --hour 12 -m ECMWF
"""

REFERENCE = """reference:
  Jolivet, R., R. Grandin, C. Lasserre, M.-P. Doin and G. Peltzer (2011), Systematic InSAR tropospheric
  phase delay corrections from global meteorological reanalysis data, Geophys. Res. Lett., 38, L17311,
  doi:10.1029/2011GL048757

  Jolivet, R., P. S. Agram, N. Y. Lin, M. Simons, M. P. Doin, G. Peltzer, and Z. Li (2014), Improving
  InSAR geodesy using global atmospheric models, Journal of Geophysical Research: Solid Earth, 119(3),
  2324-2341.
"""

DATA_INFO = """
  re-analysis_dataset        coverage   temporal_resolution    spatial_resolution      latency     analysis
------------------------------------------------------------------------------------------------------------
ERA-Interim (by ECMWF)        Global      00/06/12/18 UTC      0.75 deg (~83 km)       2-month      4D-var
MERRA(2) (by NASA Goddard)    Global      00/06/12/18 UTC      0.5*0.625 (~50 km)     2-3 weeks     3D-var

To download MERRA2, you need an Earthdata account, and pre-authorize the "NASA GESDISC DATA ARCHIVE" application
    following https://disc.gsfc.nasa.gov/earthdata-login.
"""

WEATHER_DIR_DEMO = """--weather-dir ~/data/aux
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
    parser.add_argument('-m', '--model', '-s', dest='trop_model', default='ECMWF',
                        choices={'ECMWF', 'MERRA', 'NARR', 'ERA', 'MERRA1'},
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
        msg = f'No input timeseries file, all the following options are required: \n{key_list}'
        msg += '\n\n'+EXAMPLE
        raise ValueError(msg)


    ## default values
    # Get Grib Source
    inps.trop_model = standardize_trop_model(inps.trop_model, standardWeatherModelNames)
    print('weather model: '+inps.trop_model)

    # weather_dir
    inps.weather_dir = os.path.expanduser(inps.weather_dir)
    inps.weather_dir = os.path.expandvars(inps.weather_dir)
    # Fallback value if WEATHER_DIR is not defined as environment variable
    if inps.weather_dir == '${WEATHER_DIR}':
        inps.weather_dir = './'
    print('weather data directory: '+inps.weather_dir)

    return inps


###############################################################
def check_inputs(inps):
    parser = create_parser()

    # output directories/files
    atr = dict()
    mintpy_dir = None
    if inps.timeseries_file:
        atr = readfile.read_attribute(inps.timeseries_file)
        mintpy_dir = os.path.dirname(inps.timeseries_file)
        if not inps.outfile:
            fbase = os.path.splitext(inps.timeseries_file)[0]
            inps.outfile = f'{fbase}_{inps.trop_model}.h5'
    elif inps.geom_file:
        atr = readfile.read_attribute(inps.geom_file)
        mintpy_dir = os.path.join(os.path.dirname(inps.geom_file), '..')
    else:
        mintpy_dir = os.path.abspath(os.getcwd())

    # trop_file
    inps.trop_file = os.path.join(mintpy_dir, f'inputs/{inps.trop_model}.h5')
    print(f'output tropospheric delay file: {inps.trop_file}')

    # hour
    if not inps.hour:
        if 'CENTER_LINE_UTC' in atr.keys():
            inps.hour = closest_weather_product_time(atr['CENTER_LINE_UTC'], inps.trop_model)
        else:
            parser.print_usage()
            raise Exception('no input for hour')
    print(f'time of cloest available product: {inps.hour}:00 UTC')

    # date list
    if inps.timeseries_file:
        print(f'read date list from timeseries file: {inps.timeseries_file}')
        ts_obj = timeseries(inps.timeseries_file)
        ts_obj.open(print_msg=False)
        inps.date_list = ts_obj.dateList
    elif len(inps.date_list) == 1:
        if os.path.isfile(inps.date_list[0]):
            print(f'read date list from text file: {inps.date_list[0]}')
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

    # Date list to grib file list
    inps.grib_file_list = date_list2grib_file(inps.date_list,
                                              inps.hour,
                                              inps.trop_model,
                                              inps.grib_dir)

    if 'REF_Y' in atr.keys():
        inps.ref_yx = [int(atr['REF_Y']), int(atr['REF_X'])]
        print(f'reference pixel: {inps.ref_yx}')

    # Coordinate system: geocoded or not
    inps.geocoded = False
    if 'Y_FIRST' in atr.keys():
        inps.geocoded = True
    print(f'geocoded: {inps.geocoded}')

    # Prepare DEM, inc_angle, lat/lon file for PyAPS to read
    if inps.geom_file:
        geom_obj = geometry(inps.geom_file)
        geom_obj.open()

        print('converting DEM/incAngle for PyAPS to read')
        # DEM
        dem = readfile.read(inps.geom_file, datasetName='height', print_msg=False)[0]
        inps.dem_file = 'pyapsDem.hgt'
        writefile.write(dem, inps.dem_file, metadata=atr)

        # inc_angle
        if 'incidenceAngle' in geom_obj.datasetNames:
            inps.inc_angle = readfile.read(inps.geom_file, datasetName='incidenceAngle', print_msg=False)[0]
        else:
            atr = readfile.read_attribute(inps.timeseries_file)
            inps.inc_angle = ut.incidence_angle(atr, dem=dem, dimension=0)
            inps.inc_angle = np.ones(dem.shape, dtype=np.float32) * inps.inc_angle
        inps.inc_angle_file = 'pyapsIncAngle.flt'
        writefile.write(inps.inc_angle, inps.inc_angle_file, metadata=atr)

        # latitude
        if 'latitude' in geom_obj.datasetNames:
            data = readfile.read(inps.geom_file, datasetName='latitude', print_msg=False)[0]
            print('converting lat for PyAPS to read')
            inps.lat_file = 'pyapsLat.flt'
            writefile.write(data, inps.lat_file, metadata=atr)
        else:
            inps.lat_file = None

        # longitude
        if 'longitude' in geom_obj.datasetNames:
            data = readfile.read(inps.geom_file, datasetName='longitude', print_msg=False)[0]
            print('converting lon for PyAPS to read')
            inps.lon_file = 'pyapsLon.flt'
            writefile.write(data, inps.lon_file, metadata=atr)
        else:
            inps.lon_file = None
    return inps, atr


###############################################################
def closest_weather_product_time(sar_acquisition_time, grib_source='ECMWF'):
    """Find closest available time of weather product from SAR acquisition time
    Inputs:
        sar_acquisition_time - string, SAR data acquisition time in seconds
        grib_source - string, Grib Source of weather reanalysis product
    Output:
        grib_hr - string, time of closest available weather product
    Example:
        '06' = closest_weather_product_time(atr['CENTER_LINE_UTC'])
        '12' = closest_weather_product_time(atr['CENTER_LINE_UTC'], 'NARR')
    """
    # Get hour/min of SAR acquisition time
    sar_time = float(sar_acquisition_time)

    # Find closest time in available weather products
    grib_hr_list = [0, 6, 12, 18]
    grib_hr = int(min(grib_hr_list, key=lambda x: abs(x-sar_time/3600.)))

    # Adjust time output format
    grib_hr = "%02d" % grib_hr
    return grib_hr


def standardize_trop_model(tropModel, standardWeatherModelNames):
    tropModel = tropModel.replace('-', '').upper()
    if tropModel in standardWeatherModelNames.keys():
        tropModel = standardWeatherModelNames[tropModel]
    return tropModel


def date_list2grib_file(date_list, hour, trop_model, grib_dir):
    grib_file_list = []
    for d in date_list:
        grib_file = grib_dir+'/'
        if   trop_model == 'ECMWF' :  grib_file += f'ERA-Int_{d}_{hour}.grb'
        elif trop_model == 'MERRA' :  grib_file += f'merra-{d}-{hour}.nc4'
        elif trop_model == 'NARR'  :  grib_file += f'narr-a_221_{d}_{hour}00_000.grb'
        elif trop_model == 'ERA'   :  grib_file += f'ERA_{d}_{hour}.grb'
        elif trop_model == 'MERRA1':  grib_file += f'merra-{d}-{hour}.hdf'
        grib_file_list.append(grib_file)
    return grib_file_list


def grib_file_name2trop_model_name(grib_file):
    grib_file = os.path.basename(grib_file)
    if grib_file.startswith('ERA-Int'):  trop_model = 'ECMWF'
    elif grib_file.startswith('merra'):  trop_model = 'MERRA'
    elif grib_file.startswith('narr'):   trop_model = 'NARR'
    elif grib_file.startswith('ERA_'):   trop_model = 'ERA'
    return trop_model


def check_exist_grib_file(gfile_list, print_msg=True):
    """Check input list of grib files, and return the existing ones with right size."""
    gfile_exist = ut.get_file_list(gfile_list)
    if gfile_exist:
        file_sizes = [os.path.getsize(i) for i in gfile_exist
                      if os.path.getsize(i) > 10e6]
        if file_sizes:
            comm_size = ut.most_common([i for i in file_sizes])
            if print_msg:
                print(f'common file size: {comm_size} bytes')
                print(f'number of grib files existed    : {len(gfile_exist)}')

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
                print(f'number of grib files corrupted  : {len(gfile_corrupt)}')

            for gfile in gfile_corrupt:
                os.remove(gfile)
                gfile_exist.remove(gfile)

            if print_msg:
                print('------------------------------------------------------------------------------')
    return gfile_exist


def dload_grib_pyaps(grib_file_list):
    """Download weather re-analysis grib files using PyAPS
    Parameters: grib_file_list : list of string of grib files
    Returns:    grib_file_list : list of string
    """
    print('\n------------------------------------------------------------------------------')
    print('downloading weather model data using PyAPS ...')

    # Get date list to download (skip already downloaded files)
    grib_file_exist = check_exist_grib_file(grib_file_list, print_msg=True)
    grib_file2dload = sorted(list(set(grib_file_list) - set(grib_file_exist)))
    date_list2dload = [str(re.findall(r'\d{8}', i)[0]) for i in grib_file2dload]
    print('number of grib files to download: %d' % len(date_list2dload))
    print('------------------------------------------------------------------------------\n')

    # Download grib file using PyAPS
    if len(date_list2dload) > 0:
        hour = re.findall(r'\d{8}[-_]\d{2}', grib_file2dload[0])[0].replace('-', '_').split('_')[1]
        grib_dir = os.path.dirname(grib_file2dload[0])

        # try 3 times to download, then use whatever downloaded to calculate delay
        trop_model = grib_file_name2trop_model_name(grib_file2dload[0])
        i = 0
        while i < 3:
            i += 1
            try:
                if   trop_model == 'ECMWF' :  pa.ECMWFdload( date_list2dload, hour, grib_dir)
                elif trop_model == 'MERRA' :  pa.MERRAdload( date_list2dload, hour, grib_dir)
                elif trop_model == 'NARR'  :  pa.NARRdload(  date_list2dload, hour, grib_dir)
                elif trop_model == 'ERA'   :  pa.ERAdload(   date_list2dload, hour, grib_dir)
                elif trop_model == 'MERRA1':  pa.MERRA1dload(date_list2dload, hour, grib_dir)
            except:
                pass

    grib_file_list = check_exist_grib_file(grib_file_list, print_msg=False)
    return grib_file_list


def get_delay(grib_file, inps):
    """Get delay matrix using PyAPS for one acquisition
    Inputs:
        grib_file - string, grib file path
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
    if inps.geocoded:
        aps = pa.PyAPS_geo(grib_file, inps.dem_file, grib=inps.trop_model,
                           demtype=np.float32, demfmt='RMG',
                           verb=False, Del=inps.delay_type)
    else:
        aps = pa.PyAPS_rdr(grib_file, inps.dem_file, grib=inps.trop_model,
                           demtype=np.float32, demfmt='RMG',
                           verb=False, Del=inps.delay_type)

    # estimate delay
    phs = np.zeros((aps.ny, aps.nx), dtype=np.float32)
    if not inps.geocoded and inps.lat_file is not None:
        aps.getgeodelay(phs,
                        lat=inps.lat_file,
                        lon=inps.lon_file,
                        inc=inps.inc_angle_file)
    else:
        aps.getdelay(phs, inc=0.)
        phs /= np.cos(inps.inc_angle*np.pi/180.)

    # Get relative phase delay in space
    phs -= phs[inps.ref_yx[0], inps.ref_yx[1]]
    phs *= -1    # reverse the sign for consistency between different phase correction steps/methods
    return phs


def get_delay_timeseries(inps, atr):
    """Calculate delay time-series and write it to HDF5 file.
    Parameters: inps : namespace, all input parameters
                atr  : dict, metadata to be saved in trop_file
    Returns:    trop_file : str, file name of ECMWF.h5
    """
    def get_dataset_size(fname):
        atr = readfile.read_attribute(fname)
        return (atr['LENGTH'], atr['WIDTH'])

    if (ut.run_or_skip(out_file=inps.trop_file, in_file=inps.grib_file_list, print_msg=False) == 'skip'
            and get_dataset_size(inps.trop_file) == get_dataset_size(inps.geom_file)):
        print(f'{inps.trop_file} file exists and is newer than all GRIB files, skip updating.')
    else:
        if any(i is None for i in [inps.geom_file, inps.ref_yx]):
            print('No DEM / incidenceAngle / ref_yx found, skip calculating tropospheric delays.')
            if not os.path.isfile(inps.trop_file):
                inps.trop_file = None
            return

        # calculate phase delay
        length, width = int(atr['LENGTH']), int(atr['WIDTH'])
        num_date = len(inps.grib_file_list)
        date_list = [str(re.findall(r'\d{8}', i)[0]) for i in inps.grib_file_list]
        trop_data = np.zeros((num_date, length, width), np.float32)

        print('calculating delay for each date using PyAPS (Jolivet et al., 2011; 2014) ...')
        print(f'number of grib files used: {num_date}')
        prog_bar = ptime.progressBar(maxValue=num_date)
        for i in range(num_date):
            grib_file = inps.grib_file_list[i]
            trop_data[i] = get_delay(grib_file, inps)
            prog_bar.update(i+1, suffix=os.path.basename(grib_file))
        prog_bar.close()

        # Convert relative phase delay on reference date
        try:
            inps.ref_date = atr['REF_DATE']
        except:
            inps.ref_date = date_list[0]
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

    # Delete temporary DEM file in ROI_PAC format
    if inps.geom_file:
        temp_files =[fname for fname in [inps.dem_file,
                                         inps.inc_angle_file,
                                         inps.lat_file,
                                         inps.lon_file]
                     if (fname is not None and 'pyaps' in fname)]
        if temp_files:
            print(f'delete temporary geometry files: {temp_files}')
            for temp_file in temp_files:
                os.remove(temp_file)
                os.remove(temp_file+'.rsc')
    return


###############################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    inps, atr = check_inputs(inps)

    inps.grib_file_list = dload_grib_pyaps(inps.grib_file_list)

    if inps.trop_file:
        get_delay_timeseries(inps, atr)

    if inps.dis_file:
        print('\n'+'-'*80)
        print('Applying tropospheric correction to displacement file...')
        if ut.run_or_skip(inps.cor_dis_file, [inps.dis_file, inps.tropo_file]) == 'run':
            # diff.py can handle different reference in space and time
            # e.g. the absolute delay and the double referenced time-series
            print('correcting delay for using diff.py')
            iargs = [inps.timeseries_file, inps.tropo_file, '-o', inps.outfile, '--force']
            print('diff.py', ' '.join(iargs))
            mintpy.cli.diff.main(iargs)

        else:
            print(f'Skip re-applying and use existed corrected displacement file: {inps.cor_dis_file}.')
    else:
        print('No input displacement file, skip correcting tropospheric delays.')

    return inps.outfile


###############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
