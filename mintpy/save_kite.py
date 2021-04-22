#!/usr/bin/env python3
"""
Created on Mon Mar 22 11:13:48 2021

@author: Marin Govorcin
"""
import sys
import argparse
import numpy as np
from mintpy.utils import readfile
from datetime import datetime

try:
    from kite.scene import Scene, SceneConfig
except ImportError:
    raise Exception('Kite and pyrocko are missing.')

EXAMPLE = """example:
  with Velocity [use step to export coseismic displacement]:
  save_kite.py geo/geo_velocity.h5 -g geo/geo_geometry.h5
  save_kite.py geo/geo_velocity.h5 -g geo/geo_geometry.h5 -m geo/geo_maskTempCoh.h5
  save_kite.py geo/geo_velocity.h5 -d step20210104 -g geo/geo_geometry.h5 -m geo/geo_maskTempCoh.h5
  save_kite.py geo/geo_velocity.h5 -d step20210104 -g geo/geo_geometry.h5 -m geo/geo_maskTempCoh.h5 -o dsc
  
  with Timeseries:
  save_kite.py geo/geo_timeseries_ERA5_ramp_demErr.h5 -d 20101120_20110220 -g geo/geo_geometry.h5
  save_kite.py geo/geo_timeseries_ERA5_ramp_demErr.h5 -d 20101120_20110220 -g geo/geo_geometry.h5 -m geo/geo_maskTempCoh.h5
  save_kite.py geo/geo_timeseries_ERA5_ramp_demErr.h5 -d 20101120_20110220 -g geo/geo_geometry.h5 -m geo/geo_maskTempCoh.h5 -o dsc
  
  IMPORT TO kite
  spool outfile_name              % /do quadtree,covariance/aps and than File>Save Scene and it is ready for GROND or BEAT
  
"""

d2r = np.pi/180.
r2d = 180/np.pi


def create_parser():
    parser = argparse.ArgumentParser(description='Generate KITE npz and yaml from MintPy h5 file.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='file to be converted, in geo coordinate.')
    parser.add_argument('-d', '--dset', '--dataset', dest='dset',
                        help='velocity step, date of timeseries, or date12 of unw. interferograms to be converted')
    parser.add_argument('-g', '--geom', dest='geom',
                        help='insar geometry; incidence, azimuth, height')
    parser.add_argument('-m', '--mask', dest='mask',
                        help='mask file, use mask.py to create on')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output filename')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    
    return inps

#########################################################################################################

def mintpy2kite(ifg,attr,date1,date2,inc,azi):
    print('\n---------------PREPARING KITE CONTAINER-----------')
    # fill the Kite container
    config = SceneConfig()
    config.frame.llLat = float(attr['Y_FIRST']) + float(attr['Y_STEP'])*float(attr['LENGTH'])
    config.frame.llLon = float(attr['X_FIRST'])
    config.frame.dE = float(attr['X_STEP'])
    config.frame.dN = float(attr['Y_STEP']) * -1
    config.frame.spacing = 'degree'

    config.meta.scene_title = attr['PROJECT_NAME']
    config.meta.scene_id = attr['trackNumber']

    if date1 is not None:
        config.meta.time_master = datetime.strptime(date1 + ' ' +attr['stopUTC'].split(' ')[1],'%Y%m%d %H:%M:%S.%f')
        config.meta.time_slave =  datetime.strptime(date2 + ' ' + attr['stopUTC'].split(' ')[1],'%Y%m%d %H:%M:%S.%f')


    if attr['ORBIT_DIRECTION'] == 'ASCENDING':
        config.meta.orbital_node =  'Ascending'
    else:
        config.meta.orbital_node =  'Descending'

    config.meta.wavelength =  attr['WAVELENGTH']
    config.meta.satellite_name = attr['PLATFORM']

    scene = Scene(
    	theta = np.flipud(np.pi/2 - inc * d2r),
    	phi = np.flipud(azi*d2r)+np.pi/2,
    	displacement = np.flipud(ifg),
    	config = config)

    
    print('Kite Scene info:')
    print('Scene title: {}'.format(config.meta.scene_title))
    print('Scene id: {}'.format(config.meta.scene_id))
    print('Scene orbit: {}'.format(config.meta.orbital_node))
    print('Scene platform: {}'.format(config.meta.satellite_name))
    print('Scene wavelength: {}'.format(config.meta.wavelength))
    print('Scene frame cols: {0:g}, rows: {1:g}'.format(float(attr['WIDTH']),float(attr['length'])))
    print('Scene frame first LAT: {0:.2f} | last LAT: {1:.2f}'.format(config.frame.llLat, config.frame.llLat + config.frame.dN * float(attr['LENGTH'])))
    print('Scene frame first LON: {0:.2f} | last LON: {1:.2f}'.format(config.frame.llLon, config.frame.llLon + config.frame.dE * float(attr['WIDTH'])))
    print('Scene frame spacing: x {}'.format(config.frame.dE), ' y:{}'.format(config.frame.dN), 'units: {}'.format(config.frame.spacing))
    print('Scene min disp:   {0:.2f}m | mean disp:   {1:.2f}m | max disp:    {2:.2f}m'.format(round(np.nanmin(scene.displacement),2),round(np.nanmean(scene.displacement),2),round(np.nanmax(scene.displacement),2)))
    print('Scene min  Inc:   {0:.2f}° | mean  Inc:   {1:.2f}° | max  Inc:   {2:.2f}°'.format(np.nanmin(scene.theta)*r2d,np.nanmean(scene.theta)*r2d,np.nanmax(scene.theta)*r2d))
    print('Scene min  Azi:   {0:.2f}° | mean  Azi:   {1:.2f}° | max  Azi: {2:.2f}°'.format(np.nanmin(scene.phi)*r2d,np.nanmean(scene.phi)*r2d,np.nanmax(scene.phi)*r2d))

    return scene

#########################################################################################################

def main(iargs=None):
    inps = cmd_line_parse(iargs)

    print('\n-------------------READ INPUTS -------------------')
    # check reference date
    print('Read metadata from file: {}'.format(inps.file))
    attr = readfile.read_attribute(inps.file)
    if attr['FILE_TYPE'] == 'timeseries' and inps.dset:
        inps.ref_date, inps.dset = inps.dset.split('_')
    else:
        inps.ref_date, inps.dset = None, None

    # read data
    print('Read data from file: {}'.format(inps.file))
    array, attr = readfile.read(inps.file, datasetName=inps.dset)
    if attr['FILE_TYPE'] == 'timeseries' and inps.ref_date:
        array -= readfile.read(inps.file, datasetName=inps.ref_date)[0]
    
    
    #Mask data
    if inps.mask is not None:
        mask = readfile.read(inps.mask)[0] 
        print('Masking data')
        array[mask==0] = np.nan
        
    if inps.ref_date is not None: 
       print('\nFirst  InSAR date: {}'.format(inps.ref_date))
       print('Second InSAR date: {}'.format(inps.dset))

    # output filename
    if not inps.outfile:
        inps.outfile = attr['PROJECT_NAME']
    
    # read geometry Inc, Heading
    print('\nIncidence angle read data from file: {}'.format(inps.geom))
    inc = readfile.read(inps.geom, datasetName='incidenceAngle')[0]
    print('Mean satellite incidence angle; {0:.2f}°'.format(np.nanmean(inc)))
    print('Azimuth angle read data from file: {}'.format(inps.geom))
    azi = readfile.read(inps.geom, datasetName='azimuthAngle')[0]
    print('Mean satellite heading angle; {0:.2f}°'.format(90+np.nanmean(azi)*-1))

    # prepare kite container
    scene = mintpy2kite(array,attr,inps.ref_date,inps.dset,inc,azi)

    print('\n---------------SAVING KITE CONTAINER-----------')
    print('Save KITE data in file: {}'.format(inps.outfile))

    scene.save(inps.outfile)

    # Create kite container
    return
    
if __name__ == "__main__":
    main(sys.argv[1:])

    



