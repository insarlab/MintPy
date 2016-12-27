#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2016, Yunjun Zhang                          #
# Author:  Yunjun Zhang                                    #
############################################################


import os
import sys
import argparse
import re
from datetime import datetime as dt

import h5py
import numpy as np

import pysar._readfile as readfile
import pysar.info as info


################################################################
izero = np.int16(0)
fzero = np.float32(0.0)
czero = np.complex64(0.0)

################################################################
def metadata_pysar2unavco(pysar_meta_dict,dateList):
    ## Extract UNAVCO format metadata from PySAR attributes dictionary and dateList 

    unavco_meta_dict = dict()
    
    #################################
    ##### Required metadata
    #################################

    ## mission
    unavco_meta_dict['mission'] = pysar_meta_dict['mission']
    #tmp = [atr['PLATFORM'],atr['PROJECT_NAME']]
    #if   'alos'  in tmp:  atr['mission'] = 'ALOS'
    #elif 'alos2' in tmp:  atr['mission'] = 'ALOS2'
    #elif 'csk'   in tmp:  atr['mission'] = 'CSK'
    #elif 'env'   in tmp:  atr['mission'] = 'ENV'
    #elif 'ers'   in tmp:  atr['mission'] = 'ERS'
    #elif 'rast'  in tmp:  atr['mission'] = 'RSAT'
    #elif 'rsat2' in tmp:  atr['mission'] = 'RSAT2'
    #elif 'sen'   in tmp:  atr['mission'] = 'SEN'
    #elif 's1'    in tmp:  atr['mission'] = 'SEN'
    #elif 'tdx'   in tmp:  atr['mission'] = 'TSX'
    #elif 'tsx'   in tmp:  atr['mission'] = 'TSX'
    #else: raise Exception('Cannot find mission name!')

    ## beam_mode/swath
    unavco_meta_dict['beam_mode']  = pysar_meta_dict['beam_mode']
    unavco_meta_dict['beam_swath'] = pysar_meta_dict['beam_swath']

    ## relative_orbit, or track number
    #atr_dict['relative_orbit'] = int(re.match(r'(\w+)T([0-9+])',atr['PROJECT_NAME']).groups()[1])
    unavco_meta_dict['relative_orbit'] = int(pysar_meta_dict['relative_orbit'])
    
    ## date info
    unavco_meta_dict['first_date'] = dt.strptime(dateList[0], '%Y%m%d').isoformat()[0:10]
    unavco_meta_dict['last_date']  = dt.strptime(dateList[-1],'%Y%m%d').isoformat()[0:10]

    ## footprint
    lats = [pysar_meta_dict['LAT_REF1'], pysar_meta_dict['LAT_REF3'], pysar_meta_dict['LAT_REF4'],\
            pysar_meta_dict['LAT_REF2'], pysar_meta_dict['LAT_REF1']]
    lons = [pysar_meta_dict['LON_REF1'], pysar_meta_dict['LON_REF3'], pysar_meta_dict['LON_REF4'],\
            pysar_meta_dict['LON_REF2'], pysar_meta_dict['LON_REF1']]
    unavco_meta_dict['scene_footprint'] = "POLYGON((" + ",".join([lon+' '+lat for lat,lon in zip(lats,lons)]) + "))"

    ## processing info
    unavco_meta_dict['processing_type']     = pysar_meta_dict['processing_type']
    unavco_meta_dict['processing_software'] = pysar_meta_dict['processing_software']
    unavco_meta_dict['history'] = dt.utcnow().isoformat()[0:10]

    #################################
    ##### Recommended metadata
    #################################
    unavco_meta_dict['frame'] = int(pysar_meta_dict['frame'])
    unavco_meta_dict['flight_direction'] = pysar_meta_dict['ORBIT_DIRECTION'][0].upper()
    if pysar_meta_dict['ANTENNA_SIDE'] == '-1':  unavco_meta_dict['look_direction'] = 'R'
    else:                                        unavco_meta_dict['look_direction'] = 'L'
    unavco_meta_dict['prf']         = float(pysar_meta_dict['PRF'])
    unavco_meta_dict['wavelength']  = float(pysar_meta_dict['WAVELENGTH'])

    return unavco_meta_dict

################################################################
def cmdLineParse():
    parser = argparse.ArgumentParser(description='Convert PySAR product into UNAVCO InSAR Archive format')

    parser.add_argument('timeseries', default='timeseries.h5', help='Timeseries file')

    parser.add_argument('-i','--incidence_angle', default='incidence_angle.h5', help='Incidence angle file')
    parser.add_argument('-d','--dem', default='dem.h5', help='DEM file')
    parser.add_argument('-c','--coherence', default='temporal_coherence.h5',
                        help='Coherence/correlation file, i.e. spatial_coherence.h5, temporal_coherence.h5')
    parser.add_argument('-m','--mask'default='mask.h5',help='Mask file')

    inps = parser.parse_args()
    
    return inps


################################################################
def main(argv):
    inps = cmdLineParse()

    print '\n**************** PySAR to UNAVCO **************'
    ##### Prepare Metadata
    pysar_meta_dict = readfile.read_attributes(inps.timeseries)
    k = pysar_meta_dict['FILE_TYPE']
    h5_timeseries = h5py.File(inps.timeseries,'r')
    dateList = h5_timeseries[k].keys();  dateList = sorted(dateList)
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
    DATE1 = dt.strptime(meta_dict['first_date'],'%Y-%m-%d').strftime('%Y%m%d')
    DATE2 = dt.strptime(meta_dict['last_date'], '%Y-%m-%d').strftime('%Y%m%d')
    TBASE = "%04d"%(0)
    BPERP = "%05d"%(0)
    outName = SAT+'_'+SW+'_'+RELORB+'_'+FRAME+'_'+DATE1+'-'+DATE2+'_'+TBASE+'_'+BPERP+'.he5'

    print '-----------------------------------------'
    print 'writing >>> '+outName
    f = h5py.File(outName,'w')
    group = f.create_group('timeseries')
    grid = group.create_group('GRIDS')

    ##### Write Attributes to the HDF File
    for key,value in meta_dict.iteritems():
        group.attrs[key] = value

    ##### Write Time Series Data
    print inps.timeseries
    for date in dateList:
        print date
        data = h5_timeseries[k].get(date)[:,:]
        dset = grid.create_dataset(date, data=data, compression='gzip')
        dset.attrs['Title'] = 'Time series displacement'
        dset.attrs['MissingValue'] = fzero
        dset.attrs['Units'] = 'meters'
        dset.attrs['_FillValue'] = fzero

    ##### Write Incidence_Angle
    if os.path.isfile(inps.incidence_angle):
        print inps.incidence_angle
        inc_angle, inc_angle_meta = readfile.read(inps.incidence_angle)
        dset = grid.create_dataset('incidence_angle', data=inc_angle, compression='gzip')
        dset.attrs['Title'] = 'Incidence angle'
        dset.attrs['MissingValue'] = fzero
        dset.attrs['Units'] = 'degrees'
        dset.attrs['_FillValue'] = fzero

    ##### Write DEM
    if os.path.isfile(inps.dem):
        print inps.dem
        dem, dem_meta = readfile.read(inps.dem)
        dset = grid.create_dataset('dem', data=dem, compression='gzip')
        dset.attrs['Title'] = 'Digital elevatino model'
        dset.attrs['MissingValue'] = izero
        dset.attrs['Units'] = 'meters'
        dset.attrs['_FillValue'] = izero

    ##### Write Coherence
    if os.path.isfile(inps.coherence):
        print inps.coherence
        coherence, coherence_meta = readfile.read(inps.coherence)
        dset = grid.create_dataset('coherence', data=coherence, compression='gzip')
        dset.attrs['Title'] = 'Temporal Coherence'
        dset.attrs['MissingValue'] = fzero
        dset.attrs['Units'] = 'None'
        dset.attrs['_FillValue'] = fzero

    ##### Write Mask
    if os.path.isfile(inps.mask):
        print inps.mask
        mask, mask_meta = readfile.read(inps.mask)
        dset = grid.create_dataset('mask', data=mask, compression='gzip')
        dset.attrs['Title'] = 'Mask'
        dset.attrs['MissingValue'] = izero
        dset.attrs['Units'] = 'None'
        dset.attrs['_FillValue'] = izero

    f.close()
    print 'Done.'
    return


################################################################
if __name__ == '__main__':
    main(sys.argv[:])
