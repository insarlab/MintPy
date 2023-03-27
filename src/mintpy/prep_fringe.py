############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Forrest Williams, Apr 2020         #
############################################################


import glob
import os
import xml.etree.ElementTree as ET

import h5py
import numpy as np

try:
    from osgeo import gdal
except ImportError:
    raise ImportError("Can not import gdal!")

from mintpy import subset
from mintpy.utils import (
    attribute as attr,
    isce_utils,
    ptime,
    readfile,
    utils as ut,
    writefile,
)


####################################################################################
def read_vrt_info(vrt_file):
    '''Read info from VRT file.
    Parameters: vrt_file - str, geometry vrt file
    Returns:    src_box  - tuple of 4 int, bounding box in (x0, y0, x1, y1)
                           indicating the area processed by FRInGE.
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
        msg = f'No pre-defined tag structure found in file: {vrt_file}!'
        msg += '\nPre-defined tag structure candidates:'
        for prefix in prefix_cand:
            msg += f'\n    {prefix}/SourceFilename'
        raise ValueError(msg)

    # src_box
    type_tag = root.find(prefix + '/SrcRect')
    xmin = int(type_tag.get('xOff'))
    ymin = int(type_tag.get('yOff'))
    xsize = int(type_tag.get('xSize'))
    ysize = int(type_tag.get('ySize'))
    xmax = xmin + xsize
    ymax = ymin + ysize
    src_box = (xmin, ymin, xmax, ymax)
    print(f'read bounding box from VRT file: {vrt_file} as (x0, y0, x1, y1): {src_box}')

    # source dir
    type_tag = root.find(prefix + '/SourceFilename')
    src_dir = os.path.dirname(type_tag.text)

    # in case of a (usually multilooked) vrt file missing SourceFilename field
    if not src_dir:
        src_dir = os.path.dirname(vrt_file)

    return src_box, src_dir


def prepare_metadata(meta_file, geom_src_dir, box=None, nlks_x=1, nlks_y=1):
    print('-'*50)

    # extract metadata from ISCE to MintPy (ROIPAC) format
    meta = isce_utils.extract_isce_metadata(meta_file, update_mode=False)[0]

    if 'Y_FIRST' in meta.keys():
        geom_ext = '.geo.full'
    else:
        geom_ext = '.rdr.full'

    # add LAT/LON_REF1/2/3/4, HEADING, A/RLOOKS
    meta = isce_utils.extract_geometry_metadata(
        geom_src_dir,
        meta=meta,
        box=box,
        fext_list=[geom_ext],
    )

    # apply optional user multilooking
    if nlks_x > 1:
        meta['RANGE_PIXEL_SIZE'] = str(float(meta['RANGE_PIXEL_SIZE']) * nlks_x)
        meta['RLOOKS'] = str(float(meta['RLOOKS']) * nlks_x)

    if nlks_y > 1:
        meta['AZIMUTH_PIXEL_SIZE'] = str(float(meta['AZIMUTH_PIXEL_SIZE']) * nlks_y)
        meta['ALOOKS'] = str(float(meta['ALOOKS']) * nlks_y)

    return meta


def prepare_timeseries(outfile, unw_file, metadata, processor, baseline_dir=None, box=None):
    print('-'*50)
    print(f'preparing timeseries file: {outfile}')

    # copy metadata to meta
    meta = {key : value for key, value in metadata.items()}
    phase2range = float(meta['WAVELENGTH']) / (4. * np.pi)

    # grab date list from the filename
    unw_files = sorted(glob.glob(unw_file))
    date12_list = [os.path.splitext(os.path.basename(i))[0] for i in unw_files]
    num_file = len(unw_files)
    print(f'number of unwrapped interferograms: {num_file}')

    ref_date = date12_list[0].split('_')[0]
    date_list = [ref_date] + [date12.split('_')[1] for date12 in date12_list]
    num_date = len(date_list)
    print(f'number of acquisitions: {num_date}\n{date_list}')

    # baseline info
    if baseline_dir is not None:
        # read baseline data
        baseline_dict = isce_utils.read_baseline_timeseries(
            baseline_dir,
            processor=processor,
            ref_date=ref_date,
        )

        # dict to array
        pbase = np.zeros(num_date, dtype=np.float32)
        for i in range(num_date):
            pbase_top, pbase_bottom = baseline_dict[date_list[i]]
            pbase[i] = (pbase_top + pbase_bottom) / 2.0

    # size info
    box = box if box else (0, 0, int(meta['WIDTH']), int(meta['LENGTH']))
    kwargs = dict(
        xoff=box[0],
        yoff=box[1],
        win_xsize=box[2]-box[0],
        win_ysize=box[3]-box[1],
    )

    # define dataset structure
    dates = np.array(date_list, dtype=np.string_)
    ds_name_dict = {
        "date"       : [dates.dtype, (num_date,), dates],
        "bperp"      : [np.float32,  (num_date,), pbase],
        "timeseries" : [np.float32,  (num_date, box[3]-box[1], box[2]-box[0]), None],
    }

    # initiate HDF5 file
    meta["FILE_TYPE"] = "timeseries"
    meta["UNIT"] = "m"
    meta['REF_DATE'] = ref_date
    writefile.layout_hdf5(outfile, ds_name_dict, metadata=meta)

    # writing data to HDF5 file
    print(f'writing data to HDF5 file {outfile} with a mode ...')
    with h5py.File(outfile, "a") as f:
        prog_bar = ptime.progressBar(maxValue=num_file)
        for i, unw_file in enumerate(unw_files):
            # read data using gdal
            ds = gdal.Open(unw_file, gdal.GA_ReadOnly)
            data = np.array(ds.GetRasterBand(2).ReadAsArray(**kwargs), dtype=np.float32)

            f["timeseries"][i+1] = data * phase2range
            prog_bar.update(i+1, suffix=date12_list[i])
        prog_bar.close()

        print('set value at the first acquisition to ZERO.')
        f["timeseries"][0] = 0.

    print(f'finished writing to HDF5 file: {outfile}')
    return outfile


def prepare_temporal_coherence(outfile, infile, metadata, box=None):
    print('-'*50)
    print(f'preparing temporal coherence file: {outfile}')

    # copy metadata to meta
    meta = {key : value for key, value in metadata.items()}
    meta["FILE_TYPE"] = "temporalCoherence"
    meta["UNIT"] = "1"

    # size info
    box = box if box else (0, 0, int(meta['WIDTH']), int(meta['LENGTH']))
    kwargs = dict(
        xoff=box[0],
        yoff=box[1],
        win_xsize=box[2]-box[0],
        win_ysize=box[3]-box[1],
    )

    # read data using gdal
    ds = gdal.Open(infile, gdal.GA_ReadOnly)
    data = np.array(ds.GetRasterBand(1).ReadAsArray(**kwargs), dtype=np.float32)

    print('set all data less than 0 to 0.')
    data[data < 0] = 0

    # write to HDF5 file
    writefile.write(data, outfile, metadata=meta)
    return outfile


def prepare_ps_mask(outfile, infile, metadata, box=None):
    print('-'*50)
    print(f'preparing PS mask file: {outfile}')

    # copy metadata to meta
    meta = {key : value for key, value in metadata.items()}
    meta["FILE_TYPE"] = "mask"
    meta["UNIT"] = "1"

    # size info
    box = box if box else (0, 0, int(meta['WIDTH']), int(meta['LENGTH']))
    kwargs = dict(
        xoff=box[0],
        yoff=box[1],
        win_xsize=box[2]-box[0],
        win_ysize=box[3]-box[1],
    )

    # read data using gdal
    ds = gdal.Open(infile, gdal.GA_ReadOnly)
    data = np.array(ds.GetRasterBand(1).ReadAsArray(**kwargs), dtype=np.float32)

    # write to HDF5 file
    writefile.write(data, outfile, metadata=meta)
    return outfile


def prepare_geometry(outfile, geom_dir, metadata, box, water_mask_file=None):
    print('-'*50)
    print(f'preparing geometry file: {outfile}')

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
    if water_mask_file:
        fDict['waterMask'] = water_mask_file

    # initiate dsDict
    dsDict = {}
    for dsName, fname in fDict.items():
        dsDict[dsName] = readfile.read(fname, datasetName=dsName, box=box)[0]

    dsDict['slantRangeDistance'] = ut.range_distance(meta, dimension=2)

    # write data to HDF5 file
    writefile.write(dsDict, outfile, metadata=meta)

    return outfile


def prepare_stack(outfile, unw_file, metadata, processor, baseline_dir=None, box=None):
    print('-'*50)
    print(f'preparing ifgramStack file: {outfile}')
    # copy metadata to meta
    meta = {key : value for key, value in metadata.items()}

    # get list of *.unw file
    unw_files = sorted(glob.glob(unw_file))
    num_pair = len(unw_files)
    print('number of interferograms:', num_pair)

    # get list of *.unw.conncomp file
    cc_files = [f'{x}.conncomp' for x in unw_files]
    cc_files = [x for x in cc_files if os.path.isfile(x)]
    print(f'number of connected components files: {len(cc_files)}')

    if len(cc_files) != len(unw_files):
        print('the number of *.unw and *.unw.conncomp files are NOT consistent')
        print('skip creating ifgramStack.h5 file.')
        return

    # get date info: date12_list
    date12_list = ptime.yyyymmdd_date12([os.path.basename(x).split('.')[0] for x in unw_files])

    # prepare baseline info
    if baseline_dir is not None:
        # read baseline timeseries
        baseline_dict = isce_utils.read_baseline_timeseries(baseline_dir, processor=processor)

        # calc baseline for each pair
        print('calc perp baseline pairs from time-series')
        pbase = np.zeros(num_pair, dtype=np.float32)
        for i, date12 in enumerate(date12_list):
            [date1, date2] = date12.split('_')
            pbase[i] = np.subtract(baseline_dict[date2], baseline_dict[date1]).mean()

    # size info
    box = box if box else (0, 0, int(meta['WIDTH']), int(meta['LENGTH']))
    kwargs = dict(
        xoff=box[0],
        yoff=box[1],
        win_xsize=box[2]-box[0],
        win_ysize=box[3]-box[1],
    )

    # define (and fill out some) dataset structure
    date12_arr = np.array([x.split('_') for x in date12_list], dtype=np.string_)
    drop_ifgram = np.ones(num_pair, dtype=np.bool_)
    ds_name_dict = {
        "date"             : [date12_arr.dtype, (num_pair, 2), date12_arr],
        "bperp"            : [np.float32,       (num_pair,),   pbase],
        "dropIfgram"       : [np.bool_,         (num_pair,),   drop_ifgram],
        "unwrapPhase"      : [np.float32,       (num_pair, box[3]-box[1], box[2]-box[0]), None],
        "connectComponent" : [np.float32,       (num_pair, box[3]-box[1], box[2]-box[0]), None],
    }

    # initiate HDF5 file
    meta["FILE_TYPE"] = "ifgramStack"
    writefile.layout_hdf5(outfile, ds_name_dict, metadata=meta)

    # writing data to HDF5 file
    print(f'writing data to HDF5 file {outfile} with a mode ...')
    with h5py.File(outfile, "a") as f:
        prog_bar = ptime.progressBar(maxValue=num_pair)
        for i, (unw_file, cc_file) in enumerate(zip(unw_files, cc_files)):

            # read/write *.unw file
            ds   = gdal.Open(unw_file, gdal.GA_ReadOnly)
            data = np.array(ds.GetRasterBand(2).ReadAsArray(**kwargs), dtype=np.float32)
            f["unwrapPhase"][i] = data

            # read/write *.unw.conncomp file
            ds   = gdal.Open(cc_file, gdal.GA_ReadOnly)
            data = np.array(ds.GetRasterBand(1).ReadAsArray(**kwargs), dtype=np.float32)
            f["connectComponent"][i] = data

            prog_bar.update(i+1, suffix=date12_list[i])
        prog_bar.close()

    print(f'finished writing to HDF5 file: {outfile}')
    return outfile


####################################################################################
def load_fringe(inps):
    """Load FRInGE products into MintPy."""

    # translate input options
    processor = isce_utils.get_processor(inps.metaFile)
    src_box, geom_src_dir = read_vrt_info(os.path.join(inps.geomDir, 'lat.vrt'))

    # metadata
    meta = prepare_metadata(
        inps.metaFile,
        geom_src_dir,
        src_box,
        nlks_x=inps.lks_x,
        nlks_y=inps.lks_y,
    )

    # subset - read pix_box for fringe file
    pix_box = subset.subset_input_dict2box(vars(inps), meta)[0]
    pix_box = ut.coordinate(meta).check_box_within_data_coverage(pix_box)
    print(f'input subset in y/x: {pix_box}')

    # subset - update src_box for isce file and meta
    src_box = (pix_box[0] + src_box[0],
               pix_box[1] + src_box[1],
               pix_box[2] + src_box[0],
               pix_box[3] + src_box[1])
    meta = attr.update_attribute4subset(meta, pix_box)
    print(f'input subset in y/x with respect to the VRT file: {src_box}')


    ## output directory
    for dname in [inps.outDir, os.path.join(inps.outDir, 'inputs')]:
        os.makedirs(dname, exist_ok=True)

    ## output filename
    ts_file      = os.path.join(inps.outDir, 'timeseries.h5')
    tcoh_file    = os.path.join(inps.outDir, 'temporalCoherence.h5')
    ps_mask_file = os.path.join(inps.outDir, 'maskPS.h5')
    stack_file   = os.path.join(inps.outDir, 'inputs/ifgramStack.h5')
    if 'Y_FIRST' in meta.keys():
        geom_file = os.path.join(inps.outDir, 'inputs/geometryGeo.h5')
    else:
        geom_file = os.path.join(inps.outDir, 'inputs/geometryRadar.h5')

    ## 1 - geometry (from SLC stacks before fringe, e.g. ISCE2)
    prepare_geometry(
        outfile=geom_file,
        geom_dir=geom_src_dir,
        metadata=meta,
        box=src_box,
        water_mask_file=inps.water_mask_file,
    )

    if inps.geom_only:
        return ts_file, tcoh_file, ps_mask_file, geom_file

    ## 2 - time-series (from fringe)
    prepare_timeseries(
        outfile=ts_file,
        unw_file=inps.unwFile,
        metadata=meta,
        processor=processor,
        baseline_dir=inps.baselineDir,
        box=pix_box,
    )

    ## 3 - temporal coherence and mask for PS (from fringe)
    prepare_temporal_coherence(
        outfile=tcoh_file,
        infile=inps.cohFile,
        metadata=meta,
        box=pix_box,
    )

    prepare_ps_mask(
        outfile=ps_mask_file,
        infile=inps.psMaskFile,
        metadata=meta,
        box=pix_box,
    )

    ## 4 - ifgramStack for unwrapped phase and connected components
    prepare_stack(
        outfile=stack_file,
        unw_file=inps.unwFile,
        metadata=meta,
        processor=processor,
        baseline_dir=inps.baselineDir,
        box=pix_box,
    )

    print('Done.')
    return
