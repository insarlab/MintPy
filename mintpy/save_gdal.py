############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Forrest Williams, Jun 2020                       #
############################################################


from osgeo import gdal, osr


##############################################################################
def array2raster(array, rasterName, rasterFormat, rasterOrigin, xStep, yStep, epsg=4326):

    # transform info
    cols = array.shape[1]
    rows = array.shape[0]
    originX = rasterOrigin[0]
    originY = rasterOrigin[1]
    transform = (originX, xStep, 0, originY, 0, yStep)

    # write
    driver = gdal.GetDriverByName(rasterFormat)
    print('initiate GDAL driver: {}'.format(driver.LongName))

    print('create raster band')
    print('raster row / column number: {}, {}'.format(rows, cols))
    print('raster transform info: {}'.format(transform))
    outRaster = driver.Create(rasterName, cols, rows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform(transform)

    print('write data to raster band')
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)

    print('set projection as: EPSG {}'.format(epsg))
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(epsg)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()
    print('finished writing to {}'.format(rasterName))

    return rasterName
