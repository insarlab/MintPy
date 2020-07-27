#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Forrest Williams, Apr 2020         #
############################################################


import os
import sys
import glob
import argparse
import h5py
import numpy as np
import defusedxml.ElementTree as ET

try:
    from osgeo import gdal
except ImportError:
    raise ImportError("Can not import gdal!")

from mintpy.utils import (
    ptime,
    readfile,
    writefile,
    utils0 as ut,
    isce_utils,
)


####################################################################################
EXAMPLE = """example:
  prep_fringe.py -u './PS_DS/unwrap/*.unw' -c ./PS_DS/tcorr_ds_ps.bin -g ./geometry -m ../reference/IW1.xml -b ../baselines -o ./mintpy

  ## example commands after prep_fringe.py
  reference_point.py timeseries.h5 -y 500 -x 1150
  generate_mask.py temporalCoherence.h5 -m 0.7 -o maskTempCoh.h5
  tropo_pyaps3.py -f timeseries.h5 -g inputs/geometryRadar.h5
  remove_ramp.py timeseries_ERA5.h5 -m maskTempCoh.h5 -s linear
  dem_error.py timeseries_ERA5_ramp.h5 -g inputs/geometryRadar.h5
  timeseries2velocity.py timeseries_ERA5_ramp_demErr.h5
  geocode.py velocity.h5 -l inputs/geometryRadar.h5
"""

def create_parser():
    """Command Line Parser"""
    parser = argparse.ArgumentParser(description="Prepare FRInGE products for MintPy",
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('-u', '--unw-file', dest='unwFile', type=str, default='./PS_DS/unwrap/*.unw',
                        help='path pattern of unwrapped interferograms (default: %(default)s).')
    parser.add_argument('-c', '--coh-file', dest='cohFile', type=str, default='./PS_DS/tcorr_ds_ps.bin',
                        help='temporal coherence file (default: %(default)s).')
    parser.add_argument('-g', '--geom-dir', dest='geomDir', type=str, default='./geometry',
                        help='FRInGE geometry directory (default: %(default)s).\n'
                             'This is used to grab 1) bounding box\n'
                             '                 AND 2) geometry source directory where the binary files are.')

    parser.add_argument('-m', '--meta-file', dest='metaFile', type=str, default='../reference/IW*.xml',
                        help='metadata file (default: %(default)s).\n'
                             'e.g.: ./reference/IW1.xml        for ISCE/topsStack OR\n'
                             '      ./referenceShelve/data.dat for ISCE/stripmapStack')
    parser.add_argument('-b', '--baseline-dir', dest='baselineDir', type=str, default='../baselines',
                        help='baseline directory (default: %(default)s).')

    parser.add_argument('-o', '--out-dir', dest='outDir', type=str, default='./mintpy',
                        help='output directory (default: %(default)s).')

    return parser


def cmd_line_parse(iargs = None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # in case meta_file is input as wildcard
    inps.metaFile = sorted(glob.glob(inps.metaFile))[0]

    #if all(not i for i in [inps.ifgDir, inps.geomDir, inps.cohFile, inps.baselineDir, inps.bbox]):
    #    parser.print_usage()
    #    raise SystemExit("ERROR: all of the following arguments are required: -i, -f, -g, -c, -p, -b")
    return inps


####################################################################################
def read_vrt_info(vrt_file):
    '''
    Read info from VRT file

    Parameters: vrt_file - str, geometry vrt file
    Returns:    box      - tuple of 4 int, bounding box in (x0, y0, x1, y1)
                src_dir  - str, path of geometry directory with binary data files
    '''
    root = ET.parse(vrt_file).getroot()

    # get VRT tag structure
    prefix_cand = ['VRTRasterBand/SimpleSource', 'VRTRasterBand']
    prefix_list = [prefix for prefix  in prefix_cand
                   if root.find(prefix + '/SourceFilename') is not None]
    if len(prefix_list) > 0:
        prefix = prefix_list[0]
    else:
        msg = 'No pre-defined tag structure found in file: {}!'.format(vrt_file)
        msg += '\nPre-defined tag structure candidates:'
        for prefix in prefix_cand:
            msg += '\n    {}/SourceFilename'.format(prefix)
        raise ValueError(msg)

    # box
    type_tag = root.find(prefix + '/SrcRect')
    xmin = int(type_tag.get('xOff'))
    ymin = int(type_tag.get('yOff'))
    xsize = int(type_tag.get('xSize'))
    ysize = int(type_tag.get('ySize'))
    xmax = xmin + xsize
    ymax = ymin + ysize
    box = (xmin, ymin, xmax, ymax)
    print('read bounding box from VRT file: {} as (x0, y0, x1, y1): {}'.format(vrt_file, box))

    # source dir
    type_tag = root.find(prefix + '/SourceFilename')
    src_dir = os.path.dirname(type_tag.text)

    return box, src_dir


def prepare_metadata(meta_file, geom_src_dir, box=None):
    print('-'*50)

    # extract metadata from ISCE to MintPy (ROIPAC) format
    meta = isce_utils.extract_isce_metadata(meta_file, update_mode=False)[0]

    if 'Y_FIRST' in meta.keys():
        geom_ext = '.geo.full'
    else:
        geom_ext = '.rdr.full'

    # add LAT/LON_REF1/2/3/4, HEADING, A/RLOOKS
    meta = isce_utils.extract_geometry_metadata(geom_src_dir,
                                                metadata=meta,
                                                box=box,
                                                fext_list=[geom_ext])

    # add LENGTH / WIDTH
    atr = readfile.read_attribute(os.path.join(geom_src_dir, 'lat{}'.format(geom_ext)))
    meta['LENGTH'] = atr['LENGTH']
    meta['WIDTH'] = atr['WIDTH']

    ## update metadata due to subset
    print('update metadata due to subset with bounding box')
    meta = ut.subset_attribute(meta, box)

    return meta


def prepare_timeseries(outfile, unw_file, metadata, processor, baseline_dir=None):
    print('-'*50)
    print('preparing timeseries file: {}'.format(outfile))

    # copy metadata to meta
    meta = {key : value for key, value in metadata.items()}
    phase2range = -1. * float(meta['WAVELENGTH']) / (4. * np.pi)

    # grab date list from the filename
    unw_files = sorted(glob.glob(unw_file))
    date12_list = [os.path.splitext(os.path.basename(i))[0] for i in unw_files]
    num_file = len(unw_files)
    print('number of unwrapped interferograms: {}'.format(num_file))

    ref_date = date12_list[0].split('_')[0]
    date_list = [ref_date] + [date12.split('_')[1] for date12 in date12_list]
    num_date = len(date_list)
    print('number of acquisitions: {}\n{}'.format(num_date, date_list))

    # define dataset structure
    length, width = int(meta['LENGTH']), int(meta['WIDTH'])
    dsNameDict = {
        "date"       : (np.dtype("S8"), (num_date,)),
        "timeseries" : (np.float32,     (num_date, length, width))
    }

    # baseline info
    baseline_dict = {}
    if baseline_dir is not None:
        # read baseline data
        baseline_dict = isce_utils.read_baseline_timeseries(baseline_dir,
                                                            processor=processor,
                                                            ref_date=ref_date)
        # dict to array
        pbase = np.zeros(num_date, dtype=np.float32)
        for i in range(num_date):
            pbase_top, pbase_bottom = baseline_dict[date_list[i]]
            pbase[i] = (pbase_top + pbase_bottom) / 2.0

        # update dataset structure
        dsNameDict["bperp"] = (np.float32, (num_date,))

    # initiate HDF5 file
    meta["FILE_TYPE"] = "timeseries"
    meta["UNIT"] = "m"
    meta['REF_DATE'] = ref_date
    writefile.layout_hdf5(outfile, dsNameDict, meta)

    # writing data to HDF5 file
    print('writing data to HDF5 file {} with a mode ...'.format(outfile))
    with h5py.File(outfile, "a") as f:
        f["date"][:,] = np.array([np.string_(i) for i in date_list])
        f["bperp"][:,] = pbase

        prog_bar = ptime.progressBar(maxValue=num_file)
        for i in range(num_file):
            # read data using gdal
            ds = gdal.Open(unw_files[i], gdal.GA_ReadOnly)
            data = np.array(ds.GetRasterBand(2).ReadAsArray())

            f["timeseries"][i+1] = data * phase2range
            prog_bar.update(i+1, suffix=date12_list[i])
        prog_bar.close()

        print('set value at the first acquisition to ZERO.')
        f["timeseries"][0] = 0.

    print('finished writing to HDF5 file: {}'.format(outfile))
    return outfile


def prepare_temporal_coherence(outfile, coh_file, metadata):
    print('-'*50)
    print('preparing temporal coherence file: {}'.format(outfile))

    # copy metadata to meta
    meta = {key : value for key, value in metadata.items()}
    meta["FILE_TYPE"] = "temporalCoherence"
    meta["UNIT"] = "1"

    # read data using gdal
    ds = gdal.Open(coh_file, gdal.GA_ReadOnly)
    data = np.array(ds.GetRasterBand(1).ReadAsArray(), dtype = "float32")

    print('set all data less than 0 to 0.')
    data[data < 0] = 0

    # write to HDF5 file
    writefile.write(data, outfile, metadata=meta)
    return outfile


def prepare_geometry(outfile, geom_dir, box, metadata):
    print('-'*50)
    print('preparing geometry file: {}'.format(outfile))

    # copy metadata to meta
    meta = {key : value for key, value in metadata.items()}
    meta["FILE_TYPE"] = "temporalCoherence"

    fDict = {
        'height'         : os.path.join(geom_dir, 'hgt.rdr.full'),
        'latitude'       : os.path.join(geom_dir, 'lat.rdr.full'),
        'longitude'      : os.path.join(geom_dir, 'lon.rdr.full'),
        'incidenceAngle' : os.path.join(geom_dir, 'los.rdr.full'),
        'azimuthAngle'   : os.path.join(geom_dir, 'los.rdr.full'),
        'shadowMask'     : os.path.join(geom_dir, 'shadowMask.rdr.full'),
    }

    # initiate dsDict
    dsDict = {}
    for dsName, fname in fDict.items():
        dsDict[dsName] = readfile.read(fname, datasetName=dsName, box=box)[0]

    dsDict['slantRangeDistance'] = ut.range_distance(meta, dimension=2)

    # write data to HDF5 file
    writefile.write(dsDict, outfile, metadata=meta)

    return outfile


####################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    # translate input options
    processor = isce_utils.get_processor(inps.metaFile)
    box, geom_src_dir = read_vrt_info(os.path.join(inps.geomDir, 'lat.vrt'))

    # metadata
    meta = prepare_metadata(inps.metaFile, geom_src_dir, box=box)

    ## output directory
    for dname in [inps.outDir, os.path.join(inps.outDir, 'inputs')]:
        os.makedirs(dname, exist_ok=True)

    ## output filename
    ts_file   = os.path.join(inps.outDir, 'timeseries.h5')
    tcoh_file = os.path.join(inps.outDir, 'temporalCoherence.h5')
    if 'Y_FIRST' in meta.keys():
        geom_file = os.path.join(inps.outDir, 'inputs/geometryGeo.h5')
    else:
        geom_file = os.path.join(inps.outDir, 'inputs/geometryRadar.h5')

    ## file 1 - time-series
    prepare_timeseries(outfile=ts_file,
                       unw_file=inps.unwFile,
                       metadata=meta,
                       processor=processor,
                       baseline_dir=inps.baselineDir)


    ## file 2 - temporal coherence
    prepare_temporal_coherence(outfile=tcoh_file,
                               coh_file=inps.cohFile,
                               metadata=meta)

    ## file 3 - geometry
    prepare_geometry(outfile=geom_file,
                     geom_dir=geom_src_dir,
                     box=box,
                     metadata=meta)

    return ts_file, tcoh_file, geom_file


####################################################################################
if __name__=="__main__":
    main(sys.argv[1:])
