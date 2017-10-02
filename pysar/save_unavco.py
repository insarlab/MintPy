#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2016, Yunjun Zhang                          #
# Author:  Yunjun Zhang                                    #
############################################################


import os
import sys
import argparse
import re
import datetime as dt

import h5py
import numpy as np

import pysar._readfile as readfile
import pysar.info as info


################################################################
INT_ZERO = np.int16(0)
FLOAT_ZERO = np.float32(0.0)
CPX_ZERO = np.complex64(0.0)

################################################################
def get_mission_name(meta_dict):
    '''Get mission name in UNAVCO InSAR Archive format from attribute mission/PLATFORM
    Input:  meta_dict : dict, attributes
    Output: mission   : string, mission name in standard UNAVCO format.
    '''
    mission = None

    key_list = meta_dict.keys()
    if 'mission' in key_list:
        value = meta_dict['mission'].lower()
    elif 'PLATFORM' in key_list:
        value = meta_dict['PLATFORM'].lower()
    else:
        print 'No PLATFORM nor mission attribute found, can not identify mission name.'
        print 'return None'
        return mission

    ## Convert to UNAVCO Mission name
    ## ERS, ENV, S1, RS1, RS2, CSK, TSX, JERS, ALOS, ALOS2
    if value.startswith('ers'):
        mission = 'ERS'
    elif value.startswith(('env','asar')):
        mission = 'ENV'
    elif value.startswith(('s1','sen')):
        mission = 'S1'
    elif value.startswith(('rs','rsat','radarsat')):
        mission = 'RS'
        if value.endswith('1'):
            mission += '1'
        else:
            mission += '2'
    elif value.startswith(('csk','cos')):
        mission = 'CSK'
    elif value.startswith(('tsx','tdx','terra','tandem')):
        mission = 'TSX'
    elif value.startswith('jers'):
        mission = 'JERS'
    elif value.startswith(('alos','palsar')):
        if value.endswith('2'):
            mission = 'ALOS2'
        else:
            mission = 'ALOS'
    else:
        print 'Un-recognized PLATFORM attribute: '+value
        print 'return None'
    return mission


def metadata_pysar2unavco(pysar_meta_dict,dateList):
    ## Extract UNAVCO format metadata from PySAR attributes dictionary and dateList 

    for key in pysar_meta_dict.keys():
        if 'unavco.' in key:
            pysar_meta_dict[key.split('unavco.')[1]] = pysar_meta_dict[key]

    unavco_meta_dict = dict()
    
    #################################
    ##### Required metadata
    #################################
    ##### Given manually
    ## mission
    # ERS,ENV,S1,RS1,RS2,CSK,TSX,JERS,ALOS,ALOS2
    try: unavco_meta_dict['mission'] = get_mission_name(pysar_meta_dict)
    except ValueError:
        print 'Missing required attribute: mission'

    ## beam_mode/swath
    unavco_meta_dict['beam_mode']  = pysar_meta_dict['beam_mode']
    unavco_meta_dict['beam_swath'] = pysar_meta_dict['beam_swath']

    ## relative_orbit, or track number
    #atr_dict['relative_orbit'] = int(re.match(r'(\w+)T([0-9+])',atr['PROJECT_NAME']).groups()[1])
    unavco_meta_dict['relative_orbit'] = int(pysar_meta_dict['relative_orbit'])
    
    ## processing info
    unavco_meta_dict['processing_type']     = pysar_meta_dict['processing_type']
    unavco_meta_dict['processing_software'] = pysar_meta_dict['processing_software']

    ##### Grabbed by script
    ## date info
    unavco_meta_dict['first_date'] = dt.datetime.strptime(dateList[0], '%Y%m%d').isoformat()[0:10]
    unavco_meta_dict['last_date']  = dt.datetime.strptime(dateList[-1],'%Y%m%d').isoformat()[0:10]

    ## footprint
    lons = [pysar_meta_dict['LON_REF1'], pysar_meta_dict['LON_REF3'], pysar_meta_dict['LON_REF4'],\
            pysar_meta_dict['LON_REF2'], pysar_meta_dict['LON_REF1']]
    lats = [pysar_meta_dict['LAT_REF1'], pysar_meta_dict['LAT_REF3'], pysar_meta_dict['LAT_REF4'],\
            pysar_meta_dict['LAT_REF2'], pysar_meta_dict['LAT_REF1']]
    unavco_meta_dict['scene_footprint'] = "POLYGON((" + ",".join([lon+' '+lat for lon,lat in zip(lons,lats)]) + "))"

    unavco_meta_dict['history'] = dt.datetime.utcnow().isoformat()[0:10]


    #################################
    ##### Recommended metadata
    #################################
    ##### Given manually
    try:    unavco_meta_dict['frame'] = int(pysar_meta_dict['frame'])
    except: pass
    try:    unavco_meta_dict['atmos_correct_method'] = pysar_meta_dict['atmos_correct_method']
    except: pass
    try:    unavco_meta_dict['post_processing_method'] = pysar_meta_dict['post_processing_method']
    except: unavco_meta_dict['post_processing_method'] = 'PySAR'
    try:  unavco_meta_dict['processing_dem'] = pysar_meta_dict['processing_dem']
    except: pass
    try:  unavco_meta_dict['unwrap_method'] = pysar_meta_dict['unwrap_method']
    except: pass

    ##### Grabbed by script
    try: unavco_meta_dict['flight_direction'] = pysar_meta_dict['ORBIT_DIRECTION'][0].upper()
    except: pass
    if pysar_meta_dict['ANTENNA_SIDE'] == '-1':  unavco_meta_dict['look_direction'] = 'R'
    else:                                        unavco_meta_dict['look_direction'] = 'L'
    try: unavco_meta_dict['polarization'] = pysar_meta_dict['POLARIZATION']
    except: pass
    try: unavco_meta_dict['prf'] = float(pysar_meta_dict['PRF'])
    except: pass
    try: unavco_meta_dict['wavelength'] = float(pysar_meta_dict['WAVELENGTH'])
    except: pass

    #################################
    ##### insarmaps metadata
    #################################
    # footprint for data coverage
    if 'X_FIRST' in pysar_meta_dict.keys():
        lon0 = float(pysar_meta_dict['X_FIRST'])
        lat0 = float(pysar_meta_dict['Y_FIRST'])
        lon1 = lon0 + float(pysar_meta_dict['X_STEP'])*int(pysar_meta_dict['WIDTH'])
        lat1 = lat0 + float(pysar_meta_dict['Y_STEP'])*int(pysar_meta_dict['FILE_LENGTH'])
        lons = [str(lon0), str(lon1), str(lon1), str(lon0), str(lon0)]
        lats = [str(lat0), str(lat0), str(lat1), str(lat1), str(lat0)]
        unavco_meta_dict['data_footprint'] = "POLYGON((" + ",".join([lon+' '+lat for lon,lat in zip(lons,lats)]) + "))"
    else:
        print 'Input file is not geocoded, no data_footprint without X/Y_FIRST/STEP info.'

    return unavco_meta_dict


def get_unavco_filename(timeseriesFile):
    '''Get output file name of UNAVCO InSAR Archive'''
    ##### Prepare Metadata
    pysar_meta_dict = readfile.read_attribute(timeseriesFile)
    k = pysar_meta_dict['FILE_TYPE']
    h5_timeseries = h5py.File(timeseriesFile,'r')
    dateList = sorted(h5_timeseries[k].keys())
    unavco_meta_dict = metadata_pysar2unavco(pysar_meta_dict, dateList)
    h5_timeseries.close()

    meta_dict = pysar_meta_dict.copy()
    meta_dict.update(unavco_meta_dict)

    #### Open HDF5 File
    SAT = meta_dict['mission']
    SW  = meta_dict['beam_mode']    # should be like FB08 for ALOS, need to find out, Yunjun, 2016-12-26
    RELORB = "%03d"%(int(meta_dict['relative_orbit']))
    FRAME  = "%04d"%(int(meta_dict['frame']))
    DATE1 = dt.datetime.strptime(meta_dict['first_date'],'%Y-%m-%d').strftime('%Y%m%d')
    DATE2 = dt.datetime.strptime(meta_dict['last_date'], '%Y-%m-%d').strftime('%Y%m%d')
    TBASE = "%04d"%(0)
    BPERP = "%05d"%(0)
    outName = SAT+'_'+SW+'_'+RELORB+'_'+FRAME+'_'+DATE1+'-'+DATE2+'_'+TBASE+'_'+BPERP+'.he5'

    return outName


################################################################
TEMPALTE='''
pysar.save.unavco  = auto   #[yes / update / no], auto for no, save timeseries to UNAVCO InSAR Archive format
                            #update for enabling update mode, a.k.a. put XXXXXXXX as endDate in filename if endDate < 1 year
'''

EXAMPLE='''example:
  save_unavco.py geo_timeseries_ECMWF_demErr_refDate_plane.h5 -i geo_incidenceAngle
                 -d demGeo.h5 -c geo_temporalCoherence.h5 -m geo_maskTempCoh.h5
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Convert PySAR product into UNAVCO InSAR Archive format',\
                                     formatter_class=argparse.RawDescriptionHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('timeseries', default='timeseries.h5', help='Timeseries file')

    parser.add_argument('-i','--incidence_angle', default='incidence_angle.h5', help='Incidence angle file')
    parser.add_argument('-d','--dem', default='dem.h5', help='DEM file')
    parser.add_argument('-c','--coherence', default='temporal_coherence.h5',
                        help='Coherence/correlation file, i.e. spatial_coherence.h5, temporal_coherence.h5')
    parser.add_argument('-m','--mask', default='mask.h5',help='Mask file')
    parser.add_argument('--update', action='store_true',\
                        help='Enable update mode, a.k.a. put XXXXXXXX as endDate in filename if endDate < 1 year')

    inps = parser.parse_args()
    return inps


################################################################
def main(argv):
    inps = cmdLineParse()

    #print '\n**************** Output to UNAVCO **************'
    ##### Prepare Metadata
    pysar_meta_dict = readfile.read_attribute(inps.timeseries)
    k = pysar_meta_dict['FILE_TYPE']
    h5_timeseries = h5py.File(inps.timeseries,'r')
    dateList = sorted(h5_timeseries[k].keys())
    unavco_meta_dict = metadata_pysar2unavco(pysar_meta_dict, dateList)
    print '## UNAVCO Metadata:'
    print '-----------------------------------------'
    info.print_attributes(unavco_meta_dict)

    meta_dict = pysar_meta_dict.copy()
    meta_dict.update(unavco_meta_dict)

    #### Open HDF5 File
    SAT = meta_dict['mission']
    SW  = meta_dict['beam_mode']    # should be like FB08 for ALOS, need to find out, Yunjun, 2016-12-26
    RELORB = "%03d"%(int(meta_dict['relative_orbit']))
    FRAME  = "%04d"%(int(meta_dict['frame']))
    TBASE = "%04d"%(0)
    BPERP = "%05d"%(0)
    DATE1 = dt.datetime.strptime(meta_dict['first_date'],'%Y-%m-%d').strftime('%Y%m%d')
    end_date = dt.datetime.strptime(meta_dict['last_date'], '%Y-%m-%d')
    if inps.update and (dt.datetime.utcnow() - end_date) < dt.timedelta(days=365):
        DATE2 = 'XXXXXXXX'
    else:
        DATE2 = end_date.strftime('%Y%m%d')
    outName = SAT+'_'+SW+'_'+RELORB+'_'+FRAME+'_'+DATE1+'-'+DATE2+'_'+TBASE+'_'+BPERP+'.he5'

    print '-----------------------------------------'
    print 'writing >>> '+outName
    f = h5py.File(outName,'w')
    hdfeos = f.create_group('HDFEOS')
    if 'Y_FIRST' in meta_dict.keys():
        gg_coord = hdfeos.create_group('GRIDS')
    else:
        gg_coord = hdfeos.create_group('SWATHS')
    group = gg_coord.create_group('timeseries')

    ##### Write Attributes to the HDF File
    print 'write metadata to '+str(f)
    for key,value in meta_dict.iteritems():
        f.attrs[key] = value

    print 'write data to '+str(group)
    ##### Write Time Series Data
    print 'reading file: '+inps.timeseries
    print 'number of acquisitions: %d' % len(dateList)
    for date in dateList:
        print date
        data = h5_timeseries[k].get(date)[:,:]
        dset = group.create_dataset(date, data=data, compression='gzip')
        dset.attrs['Title'] = 'Time series displacement'
        dset.attrs['MissingValue'] = FLOAT_ZERO
        dset.attrs['Units'] = 'meters'
        dset.attrs['_FillValue'] = FLOAT_ZERO

    ##### Write Incidence_Angle
    if os.path.isfile(inps.incidence_angle):
        print 'reading file: '+inps.incidence_angle
        inc_angle, inc_angle_meta = readfile.read(inps.incidence_angle)
        dset = group.create_dataset('incidence_angle', data=inc_angle, compression='gzip')
        dset.attrs['Title'] = 'Incidence angle'
        dset.attrs['MissingValue'] = FLOAT_ZERO
        dset.attrs['Units'] = 'degrees'
        dset.attrs['_FillValue'] = FLOAT_ZERO

    ##### Write DEM
    if os.path.isfile(inps.dem):
        print 'reading file: '+inps.dem
        dem, dem_meta = readfile.read(inps.dem)
        dset = group.create_dataset('dem', data=dem, compression='gzip')
        dset.attrs['Title'] = 'Digital elevatino model'
        dset.attrs['MissingValue'] = INT_ZERO
        dset.attrs['Units'] = 'meters'
        dset.attrs['_FillValue'] = INT_ZERO

    ##### Write Coherence
    if os.path.isfile(inps.coherence):
        print 'reading file: '+inps.coherence
        coherence, coherence_meta = readfile.read(inps.coherence)
        dset = group.create_dataset('coherence', data=coherence, compression='gzip')
        dset.attrs['Title'] = 'Temporal Coherence'
        dset.attrs['MissingValue'] = FLOAT_ZERO
        dset.attrs['Units'] = 'None'
        dset.attrs['_FillValue'] = FLOAT_ZERO

    ##### Write Mask
    if os.path.isfile(inps.mask):
        print 'reading file: '+inps.mask
        mask, mask_meta = readfile.read(inps.mask)
        dset = group.create_dataset('mask', data=mask, compression='gzip')
        dset.attrs['Title'] = 'Mask'
        dset.attrs['MissingValue'] = INT_ZERO
        dset.attrs['Units'] = 'None'
        dset.attrs['_FillValue'] = INT_ZERO

    f.close()
    print 'Done.'
    return


################################################################
if __name__ == '__main__':
    main(sys.argv[:])
