#!/usr/bin/env python3
# Author: Piyush Agram, Nov 2019


import os
import errno
import argparse
import h5py
import numpy as np
from osgeo import ogr, gdal
from mintpy.objects import timeseries
from mintpy.utils import readfile


EXAMPLE = """example:
  save_qgis.py timeseries_ERA5_ramp_demErr.h5 -g inputs/geometryRadar.h5 -o ts.shp
  save_qgis.py geo/geo_timeseries_ERA5_ramp_demErr.h5 -g geo/geo_geometryRadar.h5 -o ts.shp
"""

def cmdLineParse():
    '''
    Command line parser.
    '''
    parser = argparse.ArgumentParser(description='Convert to QGIS compatible ps time-series',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)
    parser.add_argument('ts_file', type=str,
                        help='time-series HDF5 file')
    parser.add_argument('-g', '--geom', dest='geom_file', type=str, required=True,
                        help='geometry HDF5 file')
    parser.add_argument('-o', '--outshp', dest='outshp', required=True, type=str,
                        help='Output shape file')
    return parser.parse_args()


def addMetadata(feature, location, attrs):
    '''
    Create one point in compatible shape format.
    '''

    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(location[0], location[1]) #Lon, Lat
    feature.SetGeometry(point)

    for k, v in attrs.items():
        feature.SetField(k, v)
    return


def gatherMintpy(ts_file, geom_file):
    '''
    Gather mintpy outputs.
    '''
    # grab aux files name from ts_file
    if os.path.basename(ts_file).startswith('geo_'):
        vel_file = 'geo_velocity.h5'
        coh_file = 'geo_temporalCoherence.h5'
        msk_file = 'geo_maskTempCoh.h5'
    else:
        vel_file = 'velocity.h5'
        coh_file = 'temporalCoherence.h5'
        msk_file = 'maskTempCoh.h5'
    
    #Can also add DEM Error here: corrected_DEM = DEM + DEM_error
    ts_dir = os.path.dirname(ts_file)
    rdict = { 'TimeSeries' : ts_file,
              'Velocity'   : os.path.join(ts_dir, vel_file),
              'Coherence'  : os.path.join(ts_dir, coh_file),
              'Mask'       : os.path.join(ts_dir, msk_file),
              'Geometry'   : geom_file,
            }

    #Check if the files exists.
    for fname in rdict.values():
        if not os.path.isfile(fname):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), fname)

    return rdict


def main(mintpyDict, outshp):
    '''
    Main driver.
    '''

    shpDriver = ogr.GetDriverByName("ESRI Shapefile")

    ##Check if shape file already exists
    if os.path.exists(outshp):
        print('Output shape file {} exists. Will be overwritten ....'.format(outshp))
        shpDriver.DeleteDataSource(outshp)

    ##Start creating shapefile dataset and layer definition
    ds = shpDriver.CreateDataSource(outshp)
    srs = ogr.osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    layer = ds.CreateLayer( 'mintpy', srs, geom_type=ogr.wkbPoint)

    #Add code for each point
    fd = ogr.FieldDefn('CODE', ogr.OFTString)
    fd.SetWidth(8)
    layer.CreateField(fd)

    #Add DEM height for each point - this could be before / after DEM error correction
    fd =  ogr.FieldDefn('HEIGHT', ogr.OFTReal)
    fd.SetWidth(7)
    fd.SetPrecision(2)
    layer.CreateField(fd)

    #Supposed to represent DEM error estimation uncertainty
    fd = ogr.FieldDefn('H_STDEV', ogr.OFTReal)
    fd.SetWidth(5)
    fd.SetPrecision(2)
    layer.CreateField(fd)

    #Estimated LOS velocity
    fd = ogr.FieldDefn('VEL', ogr.OFTReal)
    fd.SetWidth(8)
    fd.SetPrecision(2)
    layer.CreateField(fd)

    #Estimated uncertainty in velocity
    fd = ogr.FieldDefn('V_STDEV', ogr.OFTReal)
    fd.SetWidth(6)
    fd.SetPrecision(2)
    layer.CreateField(fd)

    #Temporal coherence
    fd = ogr.FieldDefn('COHERENCE', ogr.OFTReal)
    fd.SetWidth(5)
    fd.SetPrecision(3)
    layer.CreateField(fd)

    #Effective area - SqueeSAR DS / PS
    layer.CreateField( ogr.FieldDefn('EFF_AREA', ogr.OFTInteger))

    ##Time to load the dates from time-series HDF5 field and create one attribute for each date
    ts_obj = timeseries(mintpyDict['TimeSeries'])
    ts_obj.open()
    for date in ts_obj.dateList:
        fd = ogr.FieldDefn('D{0}'.format(date), ogr.OFTReal)
        fd.SetWidth(8)
        fd.SetPrecision(2)
        layer.CreateField(fd)
    layerDefn = layer.GetLayerDefn()

    ####Total number of points
    mask = readfile.read(mintpyDict['Mask'])[0]
    nValid = np.sum(mask != 0)
    print('Number of points with time-series: ', nValid)

    gdal.TermProgress_nocb(0.0)
    ###Loop over all datasets in context managers to skip close statements
    with h5py.File(mintpyDict['TimeSeries'], 'r') as tsid:
        nLines = tsid['timeseries'].shape[1]
        nPixels = tsid['timeseries'].shape[2]
        with h5py.File(mintpyDict['Velocity'], 'r') as velid:
            with h5py.File(mintpyDict['Coherence'], 'r') as cohid:
                with h5py.File(mintpyDict['Geometry'], 'r') as geomid:

                    #Start counter
                    counter = 1

                    #For each line
                    for line in range(nLines):
                        coh = cohid['temporalCoherence'][line,:].astype(np.float64)
                        vel = velid['velocity'][line,:].astype(np.float64)
                        velstd = velid['velocityStd'][line,:].astype(np.float64)
                        ts = tsid['timeseries'][:,line,:].astype(np.float64)
                        lat = geomid['latitude'][line,:].astype(np.float64)
                        lon = geomid['longitude'][line,:].astype(np.float64)
                        hgt = geomid['height'][line,:].astype(np.float64)
                        for ii in range(nPixels):
                            #If velocity is zero, dont include. What about ref pixel?
                            #Reference point is included in maskTempCoh.h5
                            if mask[line, ii] == 0:
                                continue
                            
                            #Create metadata dict
                            rdict = { 'CODE'      : hex(counter)[2:].zfill(8),
                                      'HEIGHT'    : hgt[ii],
                                      'H_STDEV'   : 0.,
                                      'VEL'       : vel[ii]*1000,
                                      'V_STDEV'   : velstd[ii]*1000,
                                      'COHERENCE' : coh[ii],
                                      'EFF_AREA'  : 1}

                            for ind, date in enumerate(ts_obj.dateList):
                                rdict['D{0}'.format(date)] = ts[ind, ii] * 1000

                            #Create feature with definition
                            feature = ogr.Feature(layerDefn)
                            addMetadata(feature, [lon[ii], lat[ii]], rdict) 
                            layer.CreateFeature(feature)
                            feature = None
                            counter = counter + 1

                        gdal.TermProgress_nocb(counter / nValid)


if __name__ == '__main__':
    '''
    Exectuable
    '''

    #Parse command line
    inps = cmdLineParse()

    #Gather MintPy Outputs
    mintpyDict = gatherMintpy(inps.ts_file, inps.geom_file)

    #Process
    main(mintpyDict, inps.outshp)
