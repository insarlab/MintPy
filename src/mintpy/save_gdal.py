############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Forrest Williams, Zhang Yunjun, Jun 2020         #
############################################################


import os
import warnings

import numpy as np
from osgeo import gdal, osr

from mintpy.utils import plot as pp, readfile, utils0 as ut

# link: https://gdal.org/drivers/raster/index.html
GDAL_DRIVER2EXT = {
    'GTiff' : '.tif',
    'ENVI'  : '',
    'GMT'   : '.grd',
    'GRIB'  : '.grb',
    'JPEG'  : '.jpg',
    'PNG'   : '.png',
}


##############################################################################
def write_gdal(data, meta, out_file, out_fmt='GTiff'):
    """Write 2D matrix into GDAL raster.

    Parameters: data     - 2D np.ndarray
                meta     - dict, metadata used to calculate the transform / epsg info:
                           X/Y_FIRST/STEP
                           EPSG / UTM_ZONE
                out_file - str, output file path
                out_fmt  - str, output GDAL file format (driver name)
                           https://gdal.org/drivers/raster/index.html
    Returns:    out_file - str, output file path
    """

    # prepare geotransform to connect image- to the geo-coordinates.
    # link: https://gdal.org/tutorials/geotransforms_tut.html
    transform = (
        float(meta['X_FIRST']), float(meta['X_STEP']), 0,
        float(meta['Y_FIRST']), 0, float(meta['Y_STEP']),
    )

    # prepare EPSG code for the Spatial Reference System (SRS)
    if 'EPSG' in meta.keys():
        epsg = int(meta['EPSG'])
    elif 'UTM_ZONE' in meta.keys():
        epsg = int(ut.utm_zone2epsg_code(meta['UTM_ZONE']))
    else:
        epsg = 4326
        msg = 'No EPSG or UTM_ZONE metadata found! '
        msg += 'Assume EPSG = 4326 (WGS84) and continue.'
        warnings.warn(msg)

    # convert boolean to uint8, as GDAL does not have a direct analogue to boolean
    if data.dtype == 'bool':
        print('convert data from boolean to uint8, as GDAL does not support boolean')
        data = np.array(data, dtype=np.uint8)

    # write file
    driver = gdal.GetDriverByName(out_fmt)
    print(f'initiate GDAL driver: {driver.LongName}')

    rows, cols = data.shape
    dtype = readfile.DATA_TYPE_NUMPY2GDAL[str(data.dtype)]
    print('create raster band:')
    print(f'  raster row / column number: {rows}, {cols}')
    print(f'  raster data type: {dtype} ({data.dtype})')
    raster = driver.Create(
        out_file,
        xsize=cols,
        ysize=rows,
        bands=1,
        eType=dtype,
    )

    print(f'set transform info: {transform}')
    raster.SetGeoTransform(transform)

    print(f'set projection as: EPSG {epsg}')
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    raster.SetProjection(srs.ExportToWkt())

    print('write data to raster band')
    band = raster.GetRasterBand(1)
    band.WriteArray(data)

    band.FlushCache()
    print(f'finished writing to {out_file}')

    return out_file


def save_gdal(inps):

    ## read data
    ftype = readfile.read_attribute(inps.file)['FILE_TYPE']

    # grab ref_date from dset
    if ftype == 'timeseries' and inps.dset and '_' in inps.dset:
        inps.ref_date, inps.dset = inps.dset.split('_')
    else:
        inps.ref_date = None

    ds_name = inps.dset if inps.dset else 'data'
    print(f'read {ds_name} from file: {inps.file}')
    data, meta = readfile.read(inps.file, datasetName=inps.dset)

    if ftype == 'timeseries' and inps.ref_date:
        print(f'read {inps.ref_date} from file: {inps.file}')
        data -= readfile.read(inps.file, datasetName=inps.ref_date)[0]

    ## write file
    # output file name
    if not inps.outfile:
        fbase = pp.auto_figure_title(inps.file, inps.dset, vars(inps))
        fext = GDAL_DRIVER2EXT.get(inps.out_format, '')
        inps.outfile = fbase + fext
    inps.outfile = os.path.abspath(inps.outfile)

    write_gdal(data, meta, out_file=inps.outfile, out_fmt=inps.out_format)

    return
