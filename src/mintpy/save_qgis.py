############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Piyush Agram, Zhang Yunjun, Nov 2019             #
############################################################


import errno
import os

import h5py
import numpy as np
from osgeo import ogr

from mintpy.objects import timeseries
from mintpy.utils import ptime, readfile, utils as ut


#########################################################################################
def add_metadata(feature, location, attrs):
    '''
    Create one point in compatible shape format.
    '''

    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(location[0], location[1]) #Lon, Lat
    feature.SetGeometry(point)

    for k, v in attrs.items():
        feature.SetField(k, v)
    return


def gather_files(ts_file, geom_file):
    '''
    Gather mintpy outputs.
    '''

    print('gather auxliary data files')
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
    fDict = {
        'TimeSeries' : ts_file,
        'Velocity'   : os.path.join(ts_dir, vel_file),
        'Coherence'  : os.path.join(ts_dir, coh_file),
        'Mask'       : os.path.join(ts_dir, msk_file),
        'Geometry'   : geom_file,
    }

    #Check if the files exists.
    for fname in fDict.values():
        if not os.path.isfile(fname):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), fname)

    for key, value in fDict.items():
        print(f'{key:<10}: {value}')
    return fDict


def read_bounding_box(pix_box, geo_box, geom_file):
    atr = readfile.read_attribute(geom_file)
    coord = ut.coordinate(atr, lookup_file=geom_file)

    if pix_box is None and geo_box is None:
        length, width = int(atr['LENGTH']), int(atr['WIDTH'])
        return (0, 0, width, length)

    if geo_box is not None:
        S, N, W, E = geo_box
        pix_box = coord.bbox_geo2radar((W, N, E, S))
        print(f'input bounding box in (S, N, W, E): {geo_box}')

    if pix_box is not None:
        pix_box = coord.check_box_within_data_coverage(pix_box)
        print(f'bounding box in (x0, y0, x1, y1): {pix_box}')

    return pix_box


def write_shape_file(fDict, shp_file, box=None, zero_first=False):
    '''Write time-series data to a shape file

    Parameters: fDict      - dict, with value for path of data files
                shp_file   - str, output filename
                box        - tuple of 4 int, in (x0, y0, x1, y1)
                zero_first - bool, set displacement at 1st acquisition to zero
    Returns:    shp_file   - str, output filename
    '''

    shpDriver = ogr.GetDriverByName("ESRI Shapefile")
    print(f'output shape file: {shp_file}')

    ##Check if shape file already exists
    if os.path.exists(shp_file):
        print(f'output shape file: {shp_file} exists, will be overwritten ....')
        shpDriver.DeleteDataSource(shp_file)

    ##Start creating shapefile dataset and layer definition
    ds = shpDriver.CreateDataSource(shp_file)
    srs = ogr.osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    layer = ds.CreateLayer('mintpy', srs, geom_type=ogr.wkbPoint)

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
    layer.CreateField(ogr.FieldDefn('EFF_AREA', ogr.OFTInteger))

    ##Time to load the dates from time-series HDF5 field and create one attribute for each date
    ts_obj = timeseries(fDict['TimeSeries'])
    ts_obj.open(print_msg=False)
    for date in ts_obj.dateList:
        fd = ogr.FieldDefn(f'D{date}', ogr.OFTReal)
        fd.SetWidth(8)
        fd.SetPrecision(2)
        layer.CreateField(fd)
    layerDefn = layer.GetLayerDefn()

    ####Total number of points
    mask = readfile.read(fDict['Mask'], box=box)[0]
    nValid = np.sum(mask != 0)
    print('number of points with time-series:', nValid)
    if zero_first:
        print('set displacement at the first acquisition to zero.')

    lats, lons = ut.get_lat_lon(ts_obj.metadata, geom_file=fDict['Geometry'], box=box)

    ###Loop over all datasets in context managers to skip close statements
    with h5py.File(fDict['TimeSeries'], 'r') as tsid:
        with h5py.File(fDict['Coherence'], 'r') as cohid:
            with h5py.File(fDict['Velocity'], 'r') as velid:
                with h5py.File(fDict['Geometry'], 'r') as geomid:

                    length = box[3] - box[1]
                    width = box[2] - box[0]

                    #Start counter
                    counter = 1
                    prog_bar = ptime.progressBar(maxValue=nValid)

                    #For each line
                    for i in range(length):
                        line = i + box[1]

                        # read data for the line
                        ts = tsid['timeseries'][:, line, box[0]:box[2]].astype(np.float64)
                        if zero_first:
                            ts -= np.tile(ts[0, :], (ts.shape[0], 1))

                        coh = cohid['temporalCoherence'][line, box[0]:box[2]].astype(np.float64)
                        vel = velid['velocity'][line, box[0]:box[2]].astype(np.float64)
                        vel_std = velid['velocityStd'][line, box[0]:box[2]].astype(np.float64)
                        hgt = geomid['height'][line, box[0]:box[2]].astype(np.float64)
                        lat = lats[i, :].astype(np.float64)
                        lon = lons[i, :].astype(np.float64)

                        for j in range(width):
                            if mask[i, j] == 0:
                                continue

                            #Create metadata dict
                            rdict = { 'CODE'      : hex(counter)[2:].zfill(8),
                                      'HEIGHT'    : hgt[j],
                                      'H_STDEV'   : 0.,
                                      'VEL'       : vel[j]*1000,
                                      'V_STDEV'   : vel_std[j]*1000,
                                      'COHERENCE' : coh[j],
                                      'EFF_AREA'  : 1}

                            for ind, date in enumerate(ts_obj.dateList):
                                rdict[f'D{date}'] = ts[ind, j] * 1000

                            #Create feature with definition
                            feature = ogr.Feature(layerDefn)
                            add_metadata(feature, [lon[j], lat[j]], rdict)
                            layer.CreateFeature(feature)
                            feature = None

                            # update counter / progress bar
                            counter += 1
                            prog_bar.update(counter, every=100, suffix=f'line {counter}/{nValid}')
                    prog_bar.close()

    print(f'finished writing to file: {shp_file}')
    return shp_file


#########################################################################################
def save_qgis(inps):

    # Read bounding box
    box = read_bounding_box(
        pix_box=inps.pix_bbox,
        geo_box=inps.geo_bbox,
        geom_file=inps.geom_file,
    )

    # Gather data files
    fDict = gather_files(inps.ts_file, inps.geom_file)

    # Write shape file
    write_shape_file(fDict, inps.shp_file, box=box, zero_first=inps.zero_first)

    return
