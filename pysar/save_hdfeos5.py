#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2016, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################
# Based on ProductArchive/he5_converters/isce2he5.py written
# by Scott Baker
#


import os, sys, re
import argparse
import datetime as dt
import string
import h5py
import numpy as np
import pysar.utils.readfile as readfile
import pysar.info as info


################################################################
BOOL_ZERO = np.bool_(0)
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

    if 'mission' in meta_dict.keys():
        value = meta_dict['mission'].lower()
    elif 'PLATFORM' in meta_dict.keys():
        value = meta_dict['PLATFORM'].lower()
    else:
        print('No PLATFORM nor mission attribute found, can not identify mission name.')
        print('return None')
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
        print('Un-recognized PLATFORM attribute: '+value)
        print('return None')
    return mission


def metadata_pysar2unavco(pysar_meta_dict,dateList):
    ## Extract UNAVCO format metadata from PySAR attributes dictionary and dateList 

    for key in pysar_meta_dict.keys():
        if 'unavco.' in key:
            pysar_meta_dict[key.split('unavco.')[1]] = pysar_meta_dict[key]
        if 'hdfEos5.' in key:
            pysar_meta_dict[key.split('hdfEos5.')[1]] = pysar_meta_dict[key]

    unavco_meta_dict = dict()

    #################################
    ##### Required metadata
    #################################
    ##### Given manually
    ## mission
    # ERS,ENV,S1,RS1,RS2,CSK,TSX,JERS,ALOS,ALOS2
    try: unavco_meta_dict['mission'] = get_mission_name(pysar_meta_dict)
    except ValueError:
        print('Missing required attribute: mission')

    ## beam_mode/swath
    unavco_meta_dict['beam_mode']  = pysar_meta_dict['beam_mode']
    try:    unavco_meta_dict['beam_swath'] = int(pysar_meta_dict['beam_swath'])
    except: unavco_meta_dict['beam_swath'] = 0

    ## relative_orbit, or track number
    #atr_dict['relative_orbit'] = int(re.match(r'(\w+)T([0-9+])',atr['PROJECT_NAME']).groups()[1])
    unavco_meta_dict['relative_orbit'] = int(pysar_meta_dict['relative_orbit'])

    ## processing info
    try:    unavco_meta_dict['processing_type'] = pysar_meta_dict['processing_type']
    except: unavco_meta_dict['processing_type'] = 'LOS_TIMESERIES'
    #unavco_meta_dict['processing_software'] = pysar_meta_dict['processing_software']

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
    if 'frame' in pysar_meta_dict.keys():
        unavco_meta_dict['frame'] = int(pysar_meta_dict['frame'])
    elif 'first_frame' in pysar_meta_dict.keys():
        unavco_meta_dict['frame'] = int(pysar_meta_dict['first_frame'])
    else:
        unavco_meta_dict['frame'] = 0

    try:    unavco_meta_dict['atmos_correct_method'] = pysar_meta_dict['atmos_correct_method']
    except: pass
    try:    unavco_meta_dict['post_processing_method'] = pysar_meta_dict['post_processing_method']
    except: unavco_meta_dict['post_processing_method'] = 'PYSAR'
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
        lat1 = lat0 + float(pysar_meta_dict['Y_STEP'])*int(pysar_meta_dict['LENGTH'])
        lons = [str(lon0), str(lon1), str(lon1), str(lon0), str(lon0)]
        lats = [str(lat0), str(lat0), str(lat1), str(lat1), str(lat0)]
        unavco_meta_dict['data_footprint'] = "POLYGON((" + ",".join([lon+' '+lat for lon,lat in zip(lons,lats)]) + "))"
    else:
        print('Input file is not geocoded, no data_footprint without X/Y_FIRST/STEP info.')

    return unavco_meta_dict


def get_hdfeos5_filename(timeseriesFile):
    '''Get output file name of HDF-EOS5 time series file'''
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
    #outName = SAT+'_'+SW+'_'+RELORB+'_'+FRAME+'_'+DATE1+'-'+DATE2+'_'+TBASE+'_'+BPERP+'.he5'
    outName = SAT+'_'+SW+'_'+RELORB+'_'+FRAME+'_'+DATE1+'-'+DATE2+'.he5'

    return outName

def read_template2inps(template_file, inps=None):
    '''Read input template options into Namespace inps'''
    if not inps:
        inps = cmdLineParse()

    print('read options from template file: '+os.path.basename(template_file))
    template = readfile.read_template(template_file)

    # Coherence-based network modification
    prefix = 'pysar.save.hdfEos5.'

    key = prefix+'update'
    if key in template.keys() and template[key] == 'yes':
        inps.update = True

    key = prefix+'subset'
    if key in template.keys() and template[key] == 'yes':
        inps.subset = True

    return inps


################################################################
TEMPALTE='''
pysar.save.hdfEos5         = auto   #[yes / no], auto for no, save timeseries to HDF-EOS5 format
pysar.save.hdfEos5.update  = auto   #[yes / no], auto for no, put XXXXXXXX as endDate in output filename
pysar.save.hdfEos5.subset  = auto   #[yes / no], auto for no, put subset range info   in output filename
'''

EXAMPLE='''example:
  save_hdfeos5.py geo_timeseries_ECMWF_demErr_refDate_plane.h5
  save_hdfeos5.py timeseries_ECMWF_demErr_refDate_plane.h5
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Convert PySAR timeseries product into HDF-EOS5 format\n'+\
                                     'https://earthdata.nasa.gov/user-resources/standards-and-references/hdf-eos5',\
                                     formatter_class=argparse.RawDescriptionHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('timeseries_file', default='timeseries.h5', help='Timeseries file')
    parser.add_argument('--template','-t', dest='template_file', default='pysarApp_template.txt', help='Template file')

    parser.add_argument('-c','--coherence', dest='coherence_file', default='temporalCoherence.h5',
                        help='Coherence/correlation file, i.e. spatial_coherence.h5, temporal_coherence.h5')
    parser.add_argument('-m','--mask', dest='mask_file', default='maskTempCoh.h5',help='Mask file')

    geom = parser.add_argument_group('Geometry','Geometry files')
    geom.add_argument('-d','--dem', dest='dem_file', default='geometryGeo.h5', help='DEM file')
    geom.add_argument('--rg-coord', dest='rg_coord_file', default='geometryGeo.h5', help='DEM file')
    geom.add_argument('--az-coord', dest='az_coord_file', default='geometryGeo.h5', help='DEM file')

    geom.add_argument('-i','--incidence_angle', dest='inc_angle_file', default='geometryGeo.h5',\
                      help='Incidence angle file')
    geom.add_argument('--az-angle', dest='head_angle_file', default='geometryGeo.h5',\
                      help='Incidence angle file')

    geom.add_argument('--slant-range', dest='slant_range_dist_file', default='geometryGeo.h5',\
                      help='Slant range distance file')
    geom.add_argument('--water-mask', dest='water_mask_file', default='geometryGeo.h5', help='Water mask file')
    geom.add_argument('--shadow-mask', dest='shadow_mask_file', default='geometryGeo.h5', help='Shadow mask file')

    parser.add_argument('--update', action='store_true',\
                        help='Enable update mode, a.k.a. put XXXXXXXX as endDate in filename if endDate < 1 year')
    parser.add_argument('--subset', action='store_true',\
                        help='Enable subset mode, a.k.a. put suffix _N31700_N32100_E130500_E131100')

    inps = parser.parse_args()
    return inps


################################################################
def main(argv):
    inps = cmdLineParse()
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    ##### Prepare Metadata
    pysar_meta_dict = readfile.read_attribute(inps.timeseries_file)
    k = pysar_meta_dict['FILE_TYPE']
    length = int(pysar_meta_dict['LENGTH'])
    width = int(pysar_meta_dict['WIDTH'])
    h5_timeseries = h5py.File(inps.timeseries_file,'r')
    dateList = sorted(h5_timeseries[k].keys())
    dateNum = len(dateList)
    dateListStr = str(dateList).translate(str.maketrans("[],u'", "     ")).strip()
    pysar_meta_dict['DATE_TIMESERIES'] = dateListStr
    
    unavco_meta_dict = metadata_pysar2unavco(pysar_meta_dict, dateList)
    print('## UNAVCO Metadata:')
    print('-----------------------------------------')
    info.print_attributes(unavco_meta_dict)

    meta_dict = pysar_meta_dict.copy()
    meta_dict.update(unavco_meta_dict)
    print('-----------------------------------------')

    ##### Open HDF5 File
    #####Get output filename
    SAT = meta_dict['mission']
    SW  = meta_dict['beam_mode']
    if meta_dict['beam_swath']:
        SW += str(meta_dict['beam_swath'])
    RELORB = "%03d"%(int(meta_dict['relative_orbit']))

    ##Frist and/or Last Frame
    frame1 = int(meta_dict['frame'])
    key = 'first_frame'
    if key in meta_dict.keys():
        frame1 = int(meta_dict[key])
    FRAME  = "%04d"%(frame1)
    key = 'last_frame'
    if key in meta_dict.keys():
        frame2 = int(meta_dict[key])
        if frame2 != frame1:
            FRAME += "_%04d"%(frame2)

    TBASE = "%04d"%(0)
    BPERP = "%05d"%(0)
    DATE1 = dt.datetime.strptime(meta_dict['first_date'],'%Y-%m-%d').strftime('%Y%m%d')
    DATE2 = dt.datetime.strptime(meta_dict['last_date'], '%Y-%m-%d').strftime('%Y%m%d')
    #end_date = dt.datetime.strptime(meta_dict['last_date'], '%Y-%m-%d')
    #if inps.update and (dt.datetime.utcnow() - end_date) < dt.timedelta(days=365):
    if inps.update:
        print('Update mode is enabled, put endDate as XXXXXXXX.')
        DATE2 = 'XXXXXXXX'

    #outName = SAT+'_'+SW+'_'+RELORB+'_'+FRAME+'_'+DATE1+'-'+DATE2+'_'+TBASE+'_'+BPERP+'.he5'
    outName = SAT+'_'+SW+'_'+RELORB+'_'+FRAME+'_'+DATE1+'_'+DATE2+'.he5'

    if inps.subset:
        print('Subset mode is enabled, put subset range info in output filename.')
        lat1 = float(meta_dict['Y_FIRST'])
        lon0 = float(meta_dict['X_FIRST'])
        lat0 = lat1 + float(meta_dict['Y_STEP']) * length
        lon1 = lon0 + float(meta_dict['X_STEP']) * width

        lat0Str = 'N%05d' % (round(lat0*1e3))
        lat1Str = 'N%05d' % (round(lat1*1e3))
        lon0Str = 'E%06d' % (round(lon0*1e3))
        lon1Str = 'E%06d' % (round(lon1*1e3))
        if lat0 < 0.0: lat0Str = 'S%05d' % (round(abs(lat0)*1e3))
        if lat1 < 0.0: lat1Str = 'S%05d' % (round(abs(lat1)*1e3))
        if lon0 < 0.0: lon0Str = 'W%06d' % (round(abs(lon0)*1e3))
        if lon1 < 0.0: lon1Str = 'W%06d' % (round(abs(lon1)*1e3))

        SUB = '_%s_%s_%s_%s' % (lat0Str, lat1Str, lon0Str, lon1Str)
        outName = os.path.splitext(outName)[0] + SUB + os.path.splitext(outName)[1]


    ##### Open HDF5 File
    print('writing >>> '+outName)
    f = h5py.File(outName,'w')
    if 'Y_FIRST' in meta_dict.keys():
        group = f.create_group('HDFEOS/GRIDS/timeseries')
    else:
        group = f.create_group('HDFEOS/SWATHS/timeseries')

    ##### Write Attributes to the HDF File
    print('write metadata to '+str(f))
    for key,value in iter(meta_dict.items()):
        f.attrs[key] = value

    ##### Write Observation - Displacement
    groupObs = group.create_group('observation')
    print('write data to '+str(groupObs))

    disDset = np.zeros((dateNum, length, width), np.float32)
    for i in range(dateNum):
        sys.stdout.write('\rreading 3D displacement from file %s: %d/%d ...' % (inps.timeseries_file, i+1, dateNum))
        sys.stdout.flush()
        disDset[i] = h5_timeseries[k].get(dateList[i])[:]
    print(' ')

    dset = groupObs.create_dataset('displacement', data=disDset, dtype=np.float32)
    dset.attrs['DATE_TIMESERIES'] = dateListStr
    dset.attrs['Title'] = 'Displacement time-series'
    dset.attrs['MissingValue'] = FLOAT_ZERO
    dset.attrs['Units'] = 'meters'
    dset.attrs['_FillValue'] = FLOAT_ZERO


    ##### Write Quality
    groupQ = group.create_group('quality')
    print('write data to '+str(groupQ))

    ## 1 - temporalCoherence
    print('reading coherence       from file: '+inps.coherence_file)
    data = readfile.read(inps.coherence_file)[0]
    dset = groupQ.create_dataset('temporalCoherence', data=data, compression='gzip')
    dset.attrs['Title'] = 'Temporal Coherence'
    dset.attrs['MissingValue'] = FLOAT_ZERO
    dset.attrs['Units'] = '1'
    dset.attrs['_FillValue'] = FLOAT_ZERO

    ## 2 - mask
    print('reading mask            from file: '+inps.mask_file)
    data = readfile.read(inps.mask_file, datasetName='mask')[0]
    dset = groupQ.create_dataset('mask', data=data, compression='gzip')
    dset.attrs['Title'] = 'Mask'
    dset.attrs['MissingValue'] = BOOL_ZERO
    dset.attrs['Units'] = '1'
    dset.attrs['_FillValue'] = BOOL_ZERO


    ##### Write Geometry
    ## Required: height, incidenceAngle
    ## Optional: rangeCoord, azimuthCoord, headingAngle, slantRangeDistance, waterMask, shadowMask
    groupGeom = group.create_group('geometry')
    print('write data to '+str(groupGeom))

    ## 1 - height
    print('reading height          from file: '+inps.dem_file)
    data = readfile.read(inps.dem_file, datasetName='height')[0]
    dset = groupGeom.create_dataset('height', data=data, compression='gzip')
    dset.attrs['Title'] = 'Digital elevatino model'
    dset.attrs['MissingValue'] = INT_ZERO
    dset.attrs['Units'] = 'meters'
    dset.attrs['_FillValue'] = INT_ZERO

    ## 2 - incidenceAngle
    print('reading incidence angle from file: '+inps.inc_angle_file)
    data = readfile.read(inps.inc_angle_file, datasetName='incidenceAngle')[0]
    dset = groupGeom.create_dataset('incidenceAngle', data=data, compression='gzip')
    dset.attrs['Title'] = 'Incidence angle'
    dset.attrs['MissingValue'] = FLOAT_ZERO
    dset.attrs['Units'] = 'degrees'
    dset.attrs['_FillValue'] = FLOAT_ZERO

    ## 3 - rangeCoord
    try:
        data = readfile.read(inps.rg_coord_file, datasetName='rangeCoord', printMsg=False)[0]
        print('reading range coord     from file: '+inps.rg_coord_file)
        dset = groupGeom.create_dataset('rangeCoord', data=data, compression='gzip')
        dset.attrs['Title'] = 'Range Coordinates'
        dset.attrs['MissingValue'] = FLOAT_ZERO
        dset.attrs['Units'] = '1'
        dset.attrs['_FillValue'] = FLOAT_ZERO
    except:
        print('No rangeCoord found in file %s' % (inps.rg_coord_file))

    ## 4 - azimuthCoord
    try:
        data = readfile.read(inps.az_coord_file, datasetName='azimuthCoord', printMsg=False)[0]
        print('reading azimuth coord   from file: '+inps.az_coord_file)
        dset = groupGeom.create_dataset('azimuthCoord', data=data, compression='gzip')
        dset.attrs['Title'] = 'Azimuth Coordinates'
        dset.attrs['MissingValue'] = FLOAT_ZERO
        dset.attrs['Units'] = '1'
        dset.attrs['_FillValue'] = FLOAT_ZERO
    except:
        print('No azimuthCoord found in file %s' % (inps.az_coord_file))

    ## 5 - headingAngle
    try:
        data = readfile.read(inps.head_angle_file, datasetName='heandingAngle', printMsg=False)[0]
        print('reading azimuth coord   from file: '+inps.head_angle_file)
        dset = groupGeom.create_dataset('heandingAngle', data=data, compression='gzip')
        dset.attrs['Title'] = 'Heanding Angle'
        dset.attrs['MissingValue'] = FLOAT_ZERO
        dset.attrs['Units'] = 'degrees'
        dset.attrs['_FillValue'] = FLOAT_ZERO
    except:
        print('No headingAngle found in file %s' % (inps.head_angle_file))

    ## 6 - slantRangeDistance
    try:
        data = readfile.read(inps.slant_range_dist_file, datasetName='slantRangeDistance', printMsg=False)[0]
        print('reading slant range distance from file: '+inps.slant_range_dist_file)
        dset = groupGeom.create_dataset('slantRangeDistance', data=data, compression='gzip')
        dset.attrs['Title'] = 'Slant Range Distance'
        dset.attrs['MissingValue'] = FLOAT_ZERO
        dset.attrs['Units'] = 'meters'
        dset.attrs['_FillValue'] = FLOAT_ZERO
    except:
        print('No slantRangeDistance found in file %s' % (inps.slant_range_dist_file))

    ## 7 - waterMask
    try:
        data = readfile.read(inps.water_mask_file, datasetName='waterMask', printMsg=False)[0]
        print('reading water mask      from file: '+inps.water_mask_file)
        dset = groupGeom.create_dataset('waterMask', data=data, compression='gzip')
        dset.attrs['Title'] = 'Water Mask'
        dset.attrs['MissingValue'] = BOOL_ZERO
        dset.attrs['Units'] = '1'
        dset.attrs['_FillValue'] = BOOL_ZERO
    except:
        print('No waterMask found in file %s' % (inps.water_mask_file))

    ## 8 - shadowMask
    try:
        data = readfile.read(inps.shadow_mask_file, datasetName='shadowMask', printMsg=False)[0]
        print('reading shadow mask     from file: '+inps.shadow_mask_file)
        dset = groupGeom.create_dataset('shadowMask', data=data, compression='gzip')
        dset.attrs['Title'] = 'Shadow Mask'
        dset.attrs['MissingValue'] = BOOL_ZERO
        dset.attrs['Units'] = '1'
        dset.attrs['_FillValue'] = BOOL_ZERO
    except:
        print('No shadowMask found in file %s' % (inps.shadow_mask_file))

    f.close()
    print('Done.')
    return


################################################################
if __name__ == '__main__':
    main(sys.argv[:])
