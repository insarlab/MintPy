#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2015, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os
import sys
import argparse
import re

try:
    import pyaps as pa
except:
    sys.exit('Cannot import pyaps into Python!')

import h5py
import numpy as np

import pysar.utils.datetime as ptime
import pysar.utils.readfile as readfile
import pysar.utils.writefile as writefile
import pysar.utils.utils as ut


###############################################################
def get_delay(grib_file, atr, inps_dict):
    '''Get delay matrix using PyAPS for one acquisition
    Inputs:
        grib_file - strng, grib file path
        atr       - dict, including the following attributes:
                    dem_file    - string, DEM file path
                    grib_source - string, Weather re-analysis data source
                    delay_type  - string, comb/dry/wet
                    ref_y/x     - string, reference pixel row/col number
                    inc_angle   - np.array, 0/1/2 D
    Output:
        phs - 2D np.array, absolute tropospheric phase delay relative to ref_y/x
    '''
    if 'X_FIRST' in atr.keys():
        aps = pa.PyAPS_geo(grib_file, inps_dict['dem_file'], grib=inps_dict['grib_source'],\
                           verb=True, Del=inps_dict['delay_type'])
    else:
        aps = pa.PyAPS_rdr(grib_file, inps_dict['dem_file'], grib=inps_dict['grib_source'],\
                           verb=True, Del=inps_dict['delay_type'])
    phs = np.zeros((aps.ny, aps.nx), dtype=np.float32)
    aps.getdelay(phs, inc=0.0)

    # Get relative phase delay in space
    yref = int(atr['ref_y'])
    xref = int(atr['ref_x'])
    phs -= phs[yref, xref]

    # project into LOS direction
    phs /= np.cos(inps_dict['inc_angle'])
    
    # reverse the sign for consistency between different phase correction steps/methods
    phs *= -1
    
    return phs


def date_list2grib_file(date_list, hour, grib_source, grib_dir):
    grib_file_list = []
    for d in date_list:
        grib_file = grib_dir+'/'
        if   grib_source == 'ECMWF' :  grib_file += 'ERA-Int_%s_%s.grb' % (d, hour)
        elif grib_source == 'ERA'   :  grib_file += 'ERA_%s_%s.grb' % (d, hour)
        elif grib_source == 'NARR'  :  grib_file += 'narr-a_221_%s_%s00_000.grb' % (d, hour)
        elif grib_source == 'MERRA' :  grib_file += 'merra-%s-%s.nc4' % (d, hour)
        elif grib_source == 'MERRA1':  grib_file += 'merra-%s-%s.hdf' % (d, hour)
        grib_file_list.append(grib_file)
    return grib_file_list


def dload_grib(date_list, hour, grib_source='ECMWF', weather_dir='./'):
    '''Download weather re-analysis grib files using PyAPS
    Inputs:
        date_list   : list of string in YYYYMMDD format
        hour        : string in HH:MM or HH format
        grib_source : string, 
        weather_dir : string,
    Output:
        grib_file_list : list of string
    '''
    ## Grib data directory
    weather_dir = os.path.abspath(weather_dir)
    grib_dir = weather_dir+'/'+grib_source
    if not os.path.isdir(grib_dir):
        print('making directory: '+grib_dir)
        os.makedirs(grib_dir)

    ## Date list to grib file list
    grib_file_list = date_list2grib_file(date_list, hour, grib_source, grib_dir)

    ## Get date list to download (skip already downloaded files)
    grib_file_existed = ut.get_file_list(grib_file_list)
    if grib_file_existed:
        grib_filesize_digit = ut.mode([len(str(os.path.getsize(i))) for i in grib_file_existed])
        grib_filesize_max2 = ut.mode([str(os.path.getsize(i))[0:2] for i in grib_file_existed])
        grib_file_corrupted = [i for i in grib_file_existed if (len(str(os.path.getsize(i))) != grib_filesize_digit or\
                                                                str(os.path.getsize(i))[0:2] != grib_filesize_max2)]
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

    ## Download grib file using PyAPS
    if   grib_source == 'ECMWF' :  pa.ECMWFdload( date_list2download, hour, grib_dir)
    elif grib_source == 'ERA'   :  pa.ERAdload(   date_list2download, hour, grib_dir)
    elif grib_source == 'NARR'  :  pa.NARRdload(  date_list2download, hour, grib_dir)
    elif grib_source == 'MERRA' :  pa.MERRAdload( date_list2download, hour, grib_dir)
    elif grib_source == 'MERRA1':  pa.MERRA1dload(date_list2download, hour, grib_dir)

    return grib_file_existed


###############################################################
EXAMPLE='''example:
  tropcor_pyaps.py timeseries.h5 -d geometryRadar.h5 -i geometryRadar.h5
  tropcor_pyaps.py timeseries.h5 -d geometryGeo.h5   -i geometryGeo.h5   --weather-dir /famelung/data/WEATHER
  tropcor_pyaps.py -d srtm1.dem -i 30 --hour 00 --ref-yx 2000 2500 --date-list date_list.txt

  tropcor_pyaps.py timeseries.h5 -d demRadar.h5 -s NARR
  tropcor_pyaps.py timeseries.h5 -d demRadar.h5 -s MERRA --delay dry -i 23
  tropcor_pyaps.py timeseries_LODcor.h5 -d demRadar.h5

  tropcor_pyaps.py -s ECMWF --hour 18 --date-list date_list.txt --download
  tropcor_pyaps.py -s ECMWF --hour 18 --date-list bl_list.txt   --download
'''

REFERENCE='''reference:
  Jolivet, R., R. Grandin, C. Lasserre, M.-P. Doin and G. Peltzer (2011), Systematic InSAR tropospheric
  phase delay corrections from global meteorological reanalysis data, Geophys. Res. Lett., 38, L17311,
  doi:10.1029/2011GL048757
'''

TEMPLATE='''
## 7. Tropospheric Delay Correction (optional and recommended)
## correct tropospheric delay using the following methods:
## a. pyaps - use weather re-analysis data (Jolivet et al., 2011, GRL, need to install PyAPS; Dee et al., 2011)
## b. height_correlation - correct stratified tropospheric delay (Doin et al., 2009, J Applied Geop)
## c. base_trop_cor - (not recommend) baseline error and stratified tropo simultaneously (Jo et al., 2010, Geo J)
pysar.troposphericDelay.method       = auto  #[pyaps / height_correlation / base_trop_cor / no], auto for pyaps
pysar.troposphericDelay.weatherModel = auto  #[ECMWF / MERRA / NARR], auto for ECMWF, for pyaps method
pysar.troposphericDelay.polyOrder    = auto  #[1 / 2 / 3], auto for 1, for height_correlation method
pysar.troposphericDelay.looks        = auto  #[1-inf], auto for 8, Number of looks to be applied to interferogram 
'''

DATA_INFO='''
  re-analysis_dataset      coverage   temporal_resolution    spatial_resolution      latency     analysis
------------------------------------------------------------------------------------------------------------
ERA-Interim (by ECMWF)      Global      00/06/12/18 UTC      0.75 deg (~83 km)       2-month       4D-var
MERRA2 (by NASA Goddard)    Global      00/06/12/18 UTC      0.5 * 0.625 (~50 km)   2-3 weeks      3D-var

To download MERRA2, you need an Earthdata account, and pre-authorize the "NASA GESDISC DATA ARCHIVE" application, following https://disc.gsfc.nasa.gov/earthdata-login.
'''


def cmdLineParse():
    parser = argparse.ArgumentParser(description='Tropospheric correction using weather models\n'+\
                                     '  PyAPS is used to download and calculate the delay for each time-series epoch.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=REFERENCE+'\n'+DATA_INFO+'\n'+EXAMPLE)

    parser.add_argument(dest='timeseries_file', nargs='?', help='timeseries HDF5 file, i.e. timeseries.h5')
    parser.add_argument('-d','--dem', dest='dem_file',\
                        help='DEM file, i.e. radar_4rlks.hgt, srtm1.dem')
    parser.add_argument('-i', dest='inc_angle', default='30',\
                        help='a file containing all incidence angles, or a number representing for the whole image.')
    parser.add_argument('--weather-dir', dest='weather_dir', \
                        help='directory to put downloaded weather data, i.e. ./../WEATHER\n'+\
                             'use directory of input timeseries_file if not specified.')
    parser.add_argument('--delay', dest='delay_type', default='comb', choices={'comb','dry','wet'},\
                        help='Delay type to calculate, comb contains both wet and dry delays')
    parser.add_argument('--download', action='store_true', help='Download weather data only.')
    parser.add_argument('--date-list', dest='date_list_file',\
                        help='Read the first column of text file as list of date to download data\n'+\
                             'in YYYYMMDD or YYMMDD format')
    parser.add_argument('--ref-yx', dest='ref_yx', type=int, nargs=2, help='reference pixel in y/x')

    parser.add_argument('-s', dest='weather_model',\
                        default='ECMWF', choices={'ECMWF','ERA-Interim','ERA','MERRA','MERRA1','NARR'},\
                        help='source of the atmospheric data.\n'+\
                             'By the time of 2018-Mar-06, ERA and ECMWF data download link is working.\n'+\
                             'NARR is working for 1979-Jan to 2014-Oct.\n'+\
                             'MERRA(2) is not working.')
    parser.add_argument('--hour', help='time of data in HH, e.g. 12, 06')

    parser.add_argument('--template', dest='template_file',\
                        help='template file with input options below:\n'+TEMPLATE)
    parser.add_argument('-o', dest='out_file', help='Output file name for trospheric corrected timeseries.')

    inps = parser.parse_args()

    # Calculate DELAY or DOWNLOAD DATA ONLY, required one of them
    if not inps.download and not inps.dem_file and ( not inps.timeseries_file or not inps.date_list_file ):
        parser.print_help()
        sys.exit(1)
    return inps


###############################################################
def main(argv):
    inps = cmdLineParse()

    k = None
    atr = dict()
    if inps.timeseries_file:
        inps.timeseries_file = ut.get_file_list([inps.timeseries_file])[0]
        atr = readfile.read_attribute(inps.timeseries_file)
        k = atr['FILE_TYPE']
    elif inps.dem_file:
        inps.dem_file = ut.get_file_list([inps.dem_file])[0]
        atr = readfile.read_attribute(inps.dem_file)
    if 'ref_y' not in atr.keys() and inps.ref_yx:
        print('No reference info found in input file, use input ref_yx: '+str(inps.ref_yx))
        atr['ref_y'] = inps.ref_yx[0]
        atr['ref_x'] = inps.ref_yx[1]

    ##Read Incidence angle: to map the zenith delay to the slant delay
    if os.path.isfile(inps.inc_angle):
        inps.inc_angle = readfile.read(inps.inc_angle, epoch='incidenceAngle')[0]
    else:
        inps.inc_angle = float(inps.inc_angle)
        print('incidence angle: '+str(inps.inc_angle))
    inps.inc_angle = inps.inc_angle*np.pi/180.0

    ##Prepare DEM file in ROI_PAC format for PyAPS to read
    if inps.dem_file:
        inps.dem_file = ut.get_file_list([inps.dem_file])[0]
        if os.path.splitext(inps.dem_file)[1] in ['.h5']:
            print('convert DEM file to ROIPAC format')
            dem, atr_dem = readfile.read(inps.dem_file, epoch='height')
            if 'Y_FIRST' in atr.keys():
                atr_dem['FILE_TYPE'] = '.dem'
            else:
                atr_dem['FILE_TYPE'] = '.hgt'
            outname = os.path.splitext(inps.dem_file)[0]+'4pyaps'+atr_dem['FILE_TYPE']
            inps.dem_file = writefile.write(dem, atr_dem, outname)

    print('*******************************************************************************')
    print('Downloading weather model data ...')

    ## Get Grib Source
    if   inps.weather_model in ['ECMWF','ERA-Interim']:   inps.grib_source = 'ECMWF'
    elif inps.weather_model == 'ERA'  :                   inps.grib_source = 'ERA'
    elif inps.weather_model == 'MERRA':                   inps.grib_source = 'MERRA'
    elif inps.weather_model == 'NARR' :                   inps.grib_source = 'NARR'
    else: raise Reception('Unrecognized weather model: '+inps.weather_model)
    print('grib source: '+inps.grib_source)

    # Get weather directory
    if not inps.weather_dir:
        if inps.timeseries_file:
            inps.weather_dir = os.path.dirname(os.path.abspath(inps.timeseries_file))+'/../WEATHER'
        elif inps.dem_file:
            inps.weather_dir = os.path.dirname(os.path.abspath(inps.dem_file))+'/../WEATHER'
        else:
            inps.weather_dir = os.path.abspath(os.getcwd())
    print('Store weather data into directory: '+inps.weather_dir)

    # Get date list to download
    if not inps.date_list_file:
        print('read date list info from: '+inps.timeseries_file)
        h5 = h5py.File(inps.timeseries_file, 'r')
        if 'timeseries' in list(h5.keys()):
            date_list = sorted(h5[k].keys())
        elif k in ['interferograms','coherence','wrapped']:
            ifgram_list = sorted(h5[k].keys())
            date12_list = ptime.list_ifgram2date12(ifgram_list)
            m_dates = [i.split('-')[0] for i in date12_list]
            s_dates = [i.split('-')[1] for i in date12_list]
            date_list = ptime.yyyymmdd(sorted(list(set(m_dates + s_dates))))
        else:
            raise ValueError('Un-support input file type:'+k)
        h5.close()
    else:
        date_list = ptime.yyyymmdd(np.loadtxt(inps.date_list_file, dtype=bytes, usecols=(0,)).astype(str).tolist())
        print('read date list info from: '+inps.date_list_file)

    # Get Acquisition time - hour
    if not inps.hour:
        inps.hour = ptime.closest_weather_product_time(atr['CENTER_LINE_UTC'], inps.grib_source)
    print('Time of cloest available product: '+inps.hour)

    ## Download data using PyAPS
    inps.grib_file_list = dload_grib(date_list, inps.hour, inps.weather_model, inps.weather_dir)

    if inps.download:
        print('Download completed, exit as planned.')
        return

    print('*******************************************************************************')
    print('Calcualting delay for each epoch.')

    ## Calculate tropo delay using pyaps
    length = int(atr['LENGTH'])
    width = int(atr['WIDTH'])
    date_num = len(date_list)
    trop_ts = np.zeros((date_num, length, width), np.float32)
    for i in range(date_num):
        grib_file = inps.grib_file_list[i] 
        date = date_list[i]
        print('calculate phase delay on %s from file %s' % (date, os.path.basename(grib_file)))
        trop_ts[i] = get_delay(grib_file, atr, vars(inps))

    ## Convert relative phase delay on reference date
    try:    ref_date = atr['REF_DATE']
    except: ref_date = date_list[0]
    print('convert to relative phase delay with reference date: '+ref_date)
    ref_idx = date_list.index(ref_date)
    trop_ts -= np.tile(trop_ts[ref_idx,:,:], (date_num, 1, 1))

    ## Write tropospheric delay to HDF5
    tropFile = inps.grib_source+'.h5'
    print('writing >>> %s' % (tropFile))
    h5trop = h5py.File(tropFile, 'w')
    group_trop = h5trop.create_group('timeseries')
    print('number of acquisitions: '+str(date_num))
    prog_bar = ptime.progress_bar(maxValue=date_num)
    for i in range(date_num):
        date = date_list[i]
        group_trop.create_dataset(date, data=trop_ts[i], compression='gzip')
        prog_bar.update(i+1, suffix=date)
    prog_bar.close()
    # Write Attributes
    for key,value in iter(atr.items()):
        group_trop.attrs[key] = value
    h5trop.close()

    ## Write corrected Time series to HDF5
    if k == 'timeseries':
        if not inps.out_file:
            inps.out_file = os.path.splitext(inps.timeseries_file)[0]+'_'+inps.grib_source+'.h5'
        print('writing >>> %s' % (inps.out_file))
        h5ts = h5py.File(inps.timeseries_file, 'r')
        h5tsCor = h5py.File(inps.out_file, 'w')    
        group_tsCor = h5tsCor.create_group('timeseries')
        print('number of acquisitions: '+str(date_num))
        prog_bar = ptime.progress_bar(maxValue=date_num)
        for i in range(date_num):
            date = date_list[i]
            ts = h5ts['timeseries'].get(date)[:]
            group_tsCor.create_dataset(date, data=ts-trop_ts[i], compression='gzip')
            prog_bar.update(i+1, suffix=date)
        prog_bar.close()
        h5ts.close()
        # Write Attributes
        for key,value in iter(atr.items()):
            group_tsCor.attrs[key] = value
        h5tsCor.close()

    # Delete temporary DEM file in ROI_PAC format
    if '4pyaps' in inps.dem_file:
        rmCmd = 'rm %s %s.rsc' % (inps.dem_file, inps.dem_file)
        print(rmCmd)
        os.system(rmCmd)
    print('Done.')
    return inps.out_file


###############################################################
if __name__ == '__main__':
    main(sys.argv[1:])

