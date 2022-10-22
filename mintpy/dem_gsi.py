############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Mar 2019                           #
############################################################


import argparse
import os

import numpy as np

from mintpy.utils import writefile

# DEHM basic info
dehm = argparse.Namespace
dehm.step = 0.4 / 3600  #decimal degree
dehm.length = 6000      #40 mins in latitude  per grid
dehm.width  = 9000      #60 mins in longitude per grid
dehm.data_type = np.float32


##################################################################################################
def write_dem_file(SNWE, dem_file, grid_dir):
    S, N, W, E = SNWE

    # get the min/max xx/yy number of the corresponding grid files
    yy_min = int(np.floor(S * 1.5 - 1) + 1)
    yy_max = int(np.ceil(N * 1.5 - 1))
    xx_min = int(np.floor(W - 100))
    xx_max = int(np.ceil(E - 100) - 1)
    num_grid_xx = int(xx_max - xx_min + 1)  #number of grids in x/longitude direction
    num_grid_yy = int(yy_max - yy_min + 1)  #number of grids in y/latitude direction
    print('number of grid files included: {} in lat direction and {} in lon direction'.format(
        num_grid_yy, num_grid_xx))

    print('stitching DEM from grid data ...')
    dem0 = np.zeros((num_grid_yy * dehm.length,
                     num_grid_xx * dehm.width), dtype=dehm.data_type)
    for yy in range(yy_max, yy_min-1, -1):
        r0 = (yy_max - yy) * dehm.length
        r1 = r0 + dehm.length
        for xx in range(xx_min, xx_max+1, 1):
            c0 = (xx - xx_min) * dehm.width
            c1 = c0 + dehm.width

            grid_file = os.path.join(grid_dir, f'{yy}{xx}.dehm')
            if os.path.isfile(grid_file):
                print('read', grid_file)
                data = np.fromfile(grid_file,
                                   dtype=dehm.data_type,
                                   count=dehm.length*dehm.width).reshape(dehm.length,
                                                                         dehm.width)
            else:
                print('warning: {} does not exists, filled with zero value and continue.')
                data = np.zeros(dehm.length, dehm.width, dtype=dehm.data_type)
            dem0[r0:r1, c0:c1] = data

    print(f'cropping based on the input SNWE: {SNWE} ...')
    grids_N = (yy_max + 1) / 1.5
    grids_W = xx_min + 100
    x_step = dehm.step
    y_step = -dehm.step

    r0 = round((N - grids_N) / y_step)
    r1 = round((S - grids_N) / y_step)
    c0 = round((W - grids_W) / x_step)
    c1 = round((E - grids_W) / x_step)
    dem = np.array(dem0[r0:r1, c0:c1], dtype=np.int16)
    print(f'file size in (row, col): {dem.shape}')

    # write to binary file
    print(f'writing {dem_file}')
    dem.tofile(dem_file)

    # generate meta namespace
    meta = argparse.Namespace
    meta.length, meta.width = dem.shape
    meta.south = S
    meta.north = N
    meta.west = W
    meta.east = E
    meta.lat_step = -1. * dehm.step
    meta.lon_step = dehm.step
    meta.file_path = os.path.abspath(dem_file)
    return meta


def write_rsc_file(meta, fname):
    # initiate meta dict
    rsc = dict()
    rsc['FILE_LENGTH'] = meta.length
    rsc['WIDTH'] = meta.width
    rsc['XMIN'] = 0
    rsc['XMAX'] = meta.width - 1
    rsc['YMIN'] = 0
    rsc['YMAX'] = meta.length - 1
    rsc['X_FIRST'] = f'{meta.west:.12f}'
    rsc['Y_FIRST'] = f'{meta.north:.12f}'
    rsc['X_STEP'] = f'{meta.lon_step:.12f}'
    rsc['Y_STEP'] = f'{meta.lat_step:.12f}'
    rsc['X_UNIT'] = 'degrees'
    rsc['Y_UNIT'] = 'degrees'
    rsc['RLOOKS'] = 1
    rsc['ALOOKS'] = 1
    rsc['Z_OFFSET'] = 0
    rsc['Z_SCALE'] = 1
    rsc['PROJECTION'] = 'LATLON'
    rsc['DATE12'] = '111111-222222'
    rsc['DATA_TYPE'] = 'int16'

    # write rsc file
    rsc_file = fname + '.rsc'
    writefile.write_roipac_rsc(rsc, rsc_file, print_msg=True)

    return rsc_file


def write_vrt_file(meta, fname):
    # initiate vrt string
    vrt_str = """<VRTDataset rasterXSize="{w}" rasterYSize="{l}">
    <SRS>EPSG:4326</SRS>
    <GeoTransform>{x0}, {xs}, 0.0, {y0}, 0.0, {ys}</GeoTransform>
    <VRTRasterBand band="1" dataType="Int16" subClass="VRTRawRasterBand">
        <SourceFilename relativeToVRT="1">{f}</SourceFilename>
        <ByteOrder>LSB</ByteOrder>
        <ImageOffset>0</ImageOffset>
        <PixelOffset>2</PixelOffset>
        <LineOffset>{lo}</LineOffset>
    </VRTRasterBand>
</VRTDataset>
""".format(w=meta.width,
           l=meta.length,
           x0=meta.west,
           xs=meta.lon_step,
           y0=meta.north,
           ys=meta.lat_step,
           f=os.path.basename(meta.file_path),
           lo=2*meta.width)

    # write to vrt file
    vrt_file = fname + '.vrt'
    with open(vrt_file, 'w') as f:
        f.write(vrt_str)
    print(f'write {vrt_file}')

    return vrt_file


def write_isce_metadata(meta, fname):
    """
    Write metadata files in ISCE format (.vrt and .xml files)
    """
    import isce
    import isceobj

    # create isce object for xml file
    img = isceobj.createDemImage()
    img.setFilename(os.path.abspath(fname))
    img.setWidth(meta.width)
    img.setLength(meta.length)
    img.setAccessMode('READ')
    img.bands = 1
    img.dataType = 'SHORT'
    img.scheme = 'BIP'
    img.reference = 'WGS84'

    img.firstLatitude  = meta.north + meta.lat_step / 2.
    img.firstLongitude = meta.west + meta.lon_step / 2.
    img.deltaLatitude  = meta.lat_step
    img.deltaLongitude = meta.lon_step

    # write to xml file
    xml_file = fname + '.xml'
    img.dump(xml_file)

    return xml_file


def add_reference_datum(xml_file):
    """
    Example of modifying an existing XML file
    """

    import xml.etree.ElementTree as ET
    from xml.dom import minidom
    print(f'add <reference> info to xml file: {os.path.basename(xml_file)}')

    # get property element for reference
    ref = ET.Element("property", attrib={'name': 'reference'})

    val = ET.SubElement(ref, "value")
    val.text = "WGS84"

    doc = ET.SubElement(ref, "doc")
    doc.text = "Geodetic datum"

    # pretty xml
    ref_str = minidom.parseString(ET.tostring(ref)).toprettyxml(indent="    ")
    ref = ET.fromstring(ref_str)

    # write back to xml file
    tree = ET.parse(xml_file)
    root = tree.getroot()
    root.append(ref)
    tree.write(xml_file)
    return xml_file


def prep_gsi_dem(inps):
    """Prepare the GSI DEM for InSAR processing."""

    meta = write_dem_file(
        inps.SNWE,
        dem_file=inps.outfile,
        grid_dir=inps.grid_dir,
    )

    # rsc file for roipac
    write_rsc_file(meta, inps.outfile)

    # vrt/xml file for isce
    write_isce_metadata(meta, inps.outfile)

    return inps.outfile
