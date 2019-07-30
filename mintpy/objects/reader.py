# Author: Heresh Fattahi

import os
import gdal
from gdalconst import GA_ReadOnly
import numpy as np
from lxml import objectify


standardMetadatKeys = {'width': 'WIDTH', 'Width': 'WIDTH', 'length': 'LENGTH', 'FILE_LENGTH': 'LENGTH',
                       'wavelength': 'WAVELENGTH', 'Wavelength': 'WAVELENGTH', 'prf': 'PRF'
                       }

GDAL2NUMPY_DATATYPE = {
    1: np.uint8,
    2: np.uint16,
    3: np.int16,
    4: np.uint32,
    5: np.int32,
    6: np.float32,
    7: np.float64,
    10: np.complex64,
    11: np.complex128,
}


def read(file, processor='ISCE', bands=None, dataType=None):
    ''' raeder based on GDAL.

    Args:

        * file      -> File name to be read

    Kwargs:

        * processor -> the processor used for the InSAR processing. default: ISCE
        * bands     -> a list of bands to be extracted. If not specified all bands will be extracted. 
        * dataType  -> if not specified, it will be extracted from the data itself
    Returns:
        * data : A numpy array with dimensions : number_of_bands * length * width
    '''

    # if processor == 'ISCE':
    #    cmd = 'isce2gis.py envi -i ' + file
    #    os.system(cmd)
    print('reading ', file)
    dataset = gdal.Open(file, GA_ReadOnly)

    ######################################
    # if the bands have not been specified, all bands will be extracted
    if bands is None:
        bands = range(1, dataset.RasterCount+1)
    ######################################
    # if dataType is not known let's get it from the data:
    if dataType is None:
        band = dataset.GetRasterBand(1)
        dataType = GDAL2NUMPY_DATATYPE[band.DataType]

    ######################################
    # Form a numpy array of zeros with the the shape of (number of bands * length * width) and a given data type
    data = np.zeros((len(bands), dataset.RasterYSize, dataset.RasterXSize), dtype=dataType)
    ######################################
    # Fill the array with the Raster bands
    idx = 0
    for i in bands:
        band = dataset.GetRasterBand(i)
        data[idx, :, :] = band.ReadAsArray()
        idx += 1

    dataset = None
    return data


def read_metadata(file, processor):

    if processor == 'ISCE':
        metadataDict = read_isce_xml(file + '.xml')
    elif processor == 'ROI_PAC':
        metadataDict = read_rsc(file + '.rsc')

    metadataDict = standardize_metadat(metadataDict, standardMetadatKeys)
    return metadataDict


def read_isce_xml(file):
    xmlDict = {}
    fObj = objectify.parse(file)
    root = fObj.getroot()
    for c in root.property:
        xmlDict[c.attrib['name']] = str(c.value)

    return xmlDict


def read_rsc(inname):
    '''Reading a ROI-PAC style RSC file.

    Args:

        * inname (str): Path to the RSC file.

    Returns:

        * rdict (dict): Dictionaty of values in RSC file.
    '''

    logging.info("PROGRESS: READING %s RSC FILE" % (inname))
    rsc_dict = dict(np.loadtxt(file, dtype=str))

    return rsc_dict


def standardize_metadat(xmlDict, standardMetadatKeys):
    keys = xmlDict.keys()
    standardKeys = standardMetadatKeys.keys()
    xmlDict_standard = {}
    for k in keys:
        if k in standardKeys:
            xmlDict_standard[standardMetadatKeys[k]] = xmlDict[k]
        else:
            xmlDict_standard[k] = xmlDict[k]

    return xmlDict_standard
