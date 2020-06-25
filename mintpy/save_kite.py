#!/usr/bin/env python3
#################################################################
# Program is used for preparing data for kite software          #
# Author: Lv Xiaoran                                            #
# Created: June 2020                                            #
#################################################################

import os
import argparse
import shutil
import numpy as np
import pyproj
import scipy.io as sio

import mintpy
from mintpy.utils import readfile, writefile, utils as ut
######################################################################################
EXAMPLE = """example:
  
  Note: roipac: theta and phi is just float, which will meet error when using scene.export_csv or scene.export_geojson
        matlab: theta and phi can be float or ndarray, float type will meet the same error as roi_pac; ndarray is OK.
        Now we suggest to use matlab format!
        The situation seems is caused by kite script, we will see whether this bug is modified in the following new version.
  
  save_kite.py geo_20171117_20200205.unw -g ./inputs/geometryRadar.h5 -t roi_pac -o geo_20171117_20200205_kite
  save_kite.py geo_20171117_20200205.unw -g ./inputs/geometryRadar.h5 -t matlab -u 38 -o geo_20171117_20200205_kite 
  
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Prepare data for Kite software',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs=1, type=str, help='geocoded unw or h5 files to be converted\n')

    parser.add_argument('-g', '--geometryRadar', dest='geometry', type=str, nargs=1,
                        help='geometry file')
    parser.add_argument('-t', '--type', dest='type', type=str, nargs=1,
                        help='processor type.\n1.roi_pac\n2.matlab')
    parser.add_argument('-u','--utm_zone', dest='utmzone', type=int, nargs='*',
                        help='utm zone')
    
    parser.add_argument('-o','--outfile',dest='outfile',nargs=1,
                        help='outfile name')
    
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

def process_roi_pac(inps):
    """prepare data for kite using roi_pac format and calculate look anlge for 4 corners from incidence angle"""
    # geometry data 
    geometryRadar = inps.geometry[0]
    print('processing geometry data {}'.format(geometryRadar))
    inc_angle = readfile.read(geometryRadar, datasetName='incidenceAngle')[0]
    lat = readfile.read(geometryRadar, datasetName='latitude')[0]
    lon = readfile.read(geometryRadar, datasetName='longitude')[0]
    
    # displacement data
    disp = inps.file[0]
    print('processing displacement data {}'.format(disp))
    disp_data, atr = readfile.read(disp)
    
    # calculate position of REF1/2/3/4
    
    lat_ref1 = float(atr['LAT_REF1'])
    lon_ref1 = float(atr['LON_REF1'])
    ref1_row, ref1_col = get_same_position(lat_ref1, lon_ref1, lat, lon)
    inc_ref1 = inc_angle[ref1_row, ref1_col]
    atr['LOOK_REF1'] = str(90.0 - inc_ref1)

    lat_ref2 = float(atr['LAT_REF2'])
    lon_ref2 = float(atr['LON_REF2'])
    ref2_row, ref2_col = get_same_position(lat_ref2, lon_ref2, lat, lon)
    inc_ref2 = inc_angle[ref2_row, ref2_col]    
    atr['LOOK_REF2'] = str(90.0 - inc_ref2)

    lat_ref3 = float(atr['LAT_REF3'])
    lon_ref3 = float(atr['LON_REF3'])
    ref3_row, ref3_col = get_same_position(lat_ref3, lon_ref3, lat, lon)
    inc_ref3 = inc_angle[ref3_row, ref3_col]
    atr['LOOK_REF3'] = str(90.0 - inc_ref3)    

    lat_ref4 = float(atr['LAT_REF4'])
    lon_ref4 = float(atr['LON_REF4'])
    ref4_row, ref4_col = get_same_position(lat_ref4, lon_ref4, lat, lon)
    inc_ref4 = inc_angle[ref4_row, ref4_col]
    atr['LOOK_REF4'] = str(90.0 - inc_ref4)   

    # write atr according to kite request
    atr['HEADING_DEG'] = atr['HEADING']
    atr['Y_UNIT'] = atr['Y_UNIT'][0:-1]
    del atr['HEADING']
    del atr['X_UNIT']
    del atr['LAT_REF1']
    del atr['LON_REF1']
   
    # output data name
    outname = inps.outfile[0] + '.unw'
    outfile = os.getcwd() + '/' + outname
    writefile.write(disp_data, out_file=outfile, metadata=atr) 
    
    return

def get_same_position(lat_ref, lon_ref, lat, lon):
    """get the [row, col] of REF1/2/3/4"""
    lat_pos = np.where(lat == lat_ref)
    lon_pos = np.where(lon == lon_ref)
    
    lat_row = lat_pos[0] 
    lat_col = lat_pos[1]
    lon_row = lon_pos[0]
    lon_col = lon_pos[1]

    pos_row = lat_row[np.in1d(lat_row, lon_row)][0]
    pos_col = lat_col[np.in1d(lat_col, lon_col)][0]
   
    return pos_row, pos_col 

def process_matlab(inps):
    """prepare data for kite using matlab format"""
    # geometry data 
    geometryRadar = inps.geometry[0]
    print('processing geometry data {}'.format(geometryRadar))
    inc_angle = readfile.read(geometryRadar, datasetName='incidenceAngle')[0]
    theta = inc_angle
    #lat = readfile.read(geometryRadar, datasetName='latitude')[0]
    #lon = readfile.read(geometryRadar, datasetName='longitude')[0]
    azi_angle = readfile.read(geometryRadar, datasetName='azimuthAngle')[0]

    # heading angle
    head_angle = ut.azimuth2heading_angle(azi_angle)
    phi = -1 * head_angle + 180
 
    # displacement data
    disp_file = inps.file[0]
    print('processing displacement data {}'.format(disp_file))
    disp_data, atr = readfile.read(disp_file)
    if atr['UNIT'] == 'radian':
        phase2range = float(atr['WAVELENGTH']) / (-4 * np.pi)
        disp_data *= phase2range 
   
    row, colm = np.shape(disp_data) 
    print('the size of data/phi(~heading)/theta(incidence) is {} {}'.format(row,colm))    
    
    utmzone = inps.utmzone[0]
    print('the utmzone is {}'.format(utmzone))
    
    # get the latitude and longitude
    lon_min = float(atr['X_FIRST'])
    lon_step = float(atr['X_STEP'])
    lon_len = int(atr['WIDTH'])

    lat_max = float(atr['Y_FIRST'])
    lat_step = float(atr['Y_STEP'])
    lat_len = int(atr['LENGTH'])

    if lon_len >= lat_len:
        lon_max = lon_min + lon_step * lon_len
        lat_min = lat_max - lat_step * lon_len
        lons = np.arange(lon_min, lon_max, lon_step)
        lats = np.arange(lat_min, lat_max, lat_step)
        lon_utm, lat_utm_tmp = latlon2UTM(lats, lons, utmzone)
        lat_utm = lat_utm_tmp[lon_len - lat_len:]
    else:
        lon_max = lon_min + lon_step * lat_len
        lat_min = lat_max - lat_step * lat_len
        lons = np.arange(lon_min, lon_max, lon_step)
        lats = np.arange(lat_min, lat_max, lat_step)
        print('lon parameters {} {} {}'.format(lon_min,lon_max,lon_step))
        print('lat parameters {} {} {}'.format(lat_min,lat_max,lat_step))
        print(np.shape(lons)[0])
        print(np.shape(lats)[0])
        lon_utm_tmp, lat_utm = latlon2UTM(lats, lons, utmzone)
        lon_utm = lon_utm_tmp[0:lon_len - lat_len]
        
    print('the length of lon_utm is {}'.format(np.shape(lon_utm)[0]))
    print('the length of lat_utm is {}'.format(np.shape(lat_utm)[0]))

    UTMZone = str(utmzone)

    #    write mat file
    #    ================== ==================== ===================== =====
    #    Property           Matlab ``.mat`` name type                  unit
    #    ================== ==================== ===================== =====
    #    Scene.displacement ``ig_``              n x m array           [m]
    #    Scene.phi          ``phi``              float or n x m array  [rad]
    #    Scene.theta        ``theta``            float or n x m array  [rad]
    #    Scene.frame.x      ``xx``               n x 1 vector          [m]
    #    Scene.frame.y      ``yy``               m x 1 vector          [m]
    #    Scene.utm_zone     ``utm_zone``         str ('33T')
    #    ================== ==================== ===================== =====  
    mdict = {}
    # required by GBIS
    mdict['ig_'] = disp_data
    mdict['phi'] = phi
    mdict['theta'] = theta
    mdict['xx'] = lon_utm
    mdict['yy'] = lat_utm
    mdict['utm_zone'] = UTMZone
    
    # save to mat file
    outfile = inps.outfile[0] + '.mat'
    sio.savemat(outfile, mdict, long_field_names=True)
    print('save to file: {}.mat'.format(os.path.abspath(inps.outfile)))
    return     

def latlon2UTM(lat, lon, UTMZone):
    """convert GCS (Geographic Coordinate System) to UMT project"""
    # GCS 
    p1 = pyproj.Proj(proj='latlong', datum='WGS84')
    
    # UTM
    p2 = pyproj.Proj(proj='utm', zone=UTMZone, datum='WGS84')    
    
    # convert
    x1, y1 = p1(lon,lat)
    lon_utm, lat_utm = pyproj.transform(p1, p2, x1, y1)

    return lon_utm, lat_utm
######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    
    # generate dem.jpeg
    if inps.type[0] == 'roi_pac':
        process_roi_pac(inps)
    elif inps.type[0] == 'matlab':
        process_matlab(inps)
######################################################################################
if __name__ == '__main__':
    main()
