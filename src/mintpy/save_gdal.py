############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Forrest Williams, Zhang Yunjun, Jun 2020         #
############################################################


import os
import warnings

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
def array2raster(array, out_file, transform, epsg, out_fmt='GTiff'):
    """Write 2D matrix into gdal raster.

    Parameters: array     - 2D np.ndarray
                out_file  - str, output file path
                transform - tuple(float), geotransform to connect image- to the geo-coordinates.
                            https://gdal.org/tutorials/geotransforms_tut.html
                epsg      - int, EPSG code for the Spatial Reference System (SRS)
                out_fmt   - str, gdal driver name
    Returns:    out_file  - str, output file path
    """

    driver = gdal.GetDriverByName(out_fmt)
    print(f'initiate GDAL driver: {driver.LongName}')

    rows, cols = array.shape
    print('create raster band')
    print(f'raster row / column number: {rows}, {cols}')
    print(f'raster transform info: {transform}')
    raster = driver.Create(out_file, cols, rows, 1, gdal.GDT_Float32)
    raster.SetGeoTransform(transform)

    print('write data to raster band')
    band = raster.GetRasterBand(1)
    band.WriteArray(array)

    print(f'set projection as: EPSG {epsg}')
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    raster.SetProjection(srs.ExportToWkt())

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
    data, atr = readfile.read(inps.file, datasetName=inps.dset)

    if ftype == 'timeseries' and inps.ref_date:
        print(f'read {inps.ref_date} from file: {inps.file}')
        data -= readfile.read(inps.file, datasetName=inps.ref_date)[0]

    ## prepare the output
    # output file name
    if not inps.outfile:
        fbase = pp.auto_figure_title(inps.file, inps.dset, vars(inps))
        inps.outfile = fbase + GDAL_DRIVER2EXT.get(inps.out_format, '')
    else:
        inps.outfile = os.path.abspath(inps.outfile)

    # geotransform
    transform = (
        float(atr['X_FIRST']), float(atr['X_STEP']), 0,
        float(atr['Y_FIRST']), 0, float(atr['Y_STEP']),
    )

    # epsg
    if 'EPSG' in atr.keys():
        epsg = int(atr['EPSG'])
    elif 'UTM_ZONE' in atr.keys():
        epsg = int(ut.utm_zone2epsg_code(atr['UTM_ZONE']))
    else:
        epsg = 4326
        msg = 'No EPSG or UTM_ZONE metadata found! '
        msg += 'Assume EPSG = 4326 (WGS84) and continue.'
        warnings.warn(msg)

    ## write gdal raster
    array2raster(
        data,
        out_file=inps.outfile,
        transform=transform,
        epsg=epsg,
        out_fmt=inps.out_format,
    )

    return
