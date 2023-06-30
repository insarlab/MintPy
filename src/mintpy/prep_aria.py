############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, Zhang Yunjun, Emre Havazli, 2019 #
############################################################


import os
import time

import h5py
import numpy as np

try:
    from osgeo import gdal
except ImportError:
    raise ImportError('Can not import gdal [version>=3.0]!')

from mintpy.multilook import multilook_data
from mintpy.objects import geometry, ifgramStack, sensor
from mintpy.subset import read_subset_template2box
from mintpy.utils import attribute as attr, ptime, utils as ut, writefile


####################################################################################
def run_or_skip(inps, ds_name_dict, out_file):
    flag = 'run'

    # check 1 - update mode status
    if not inps.updateMode:
        return flag

    # check 2 - output file existence
    if ut.run_or_skip(out_file, readable=True) == 'run':
        return flag

    # check 3 - output dataset info
    key = [i for i in ['unwrapPhase', 'height'] if i in ds_name_dict.keys()][0]
    ds_shape = ds_name_dict[key][1]
    in_shape = ds_shape[-2:]

    if 'unwrapPhase' in ds_name_dict.keys():
        # compare date12 and size
        ds = gdal.Open(inps.unwFile, gdal.GA_ReadOnly)
        in_date12_list = [ds.GetRasterBand(i+1).GetMetadata("unwrappedPhase")['Dates']
                          for i in range(ds_shape[0])]
        in_date12_list = ['_'.join(d.split('_')[::-1]) for d in in_date12_list]

        try:
            out_obj = ifgramStack(out_file)
            out_obj.open(print_msg=False)
            out_shape = (out_obj.length, out_obj.width)
            out_date12_list = out_obj.get_date12_list(dropIfgram=False)

            if out_shape == in_shape and set(in_date12_list).issubset(set(out_date12_list)):
                print('All date12   exists in file {} with same size as required,'
                      ' no need to re-load.'.format(os.path.basename(out_file)))
                flag = 'skip'
        except:
            pass

    elif 'height' in ds_name_dict.keys():
        # compare dataset names and size
        in_dsNames = list(ds_name_dict.keys())
        in_size = in_shape[0] * in_shape[1] * 4 * len(in_dsNames)

        out_obj = geometry(out_file)
        out_obj.open(print_msg=False)
        out_dsNames = out_obj.datasetNames
        out_shape = (out_obj.length, out_obj.width)
        out_size = os.path.getsize(out_file)

        if (set(in_dsNames).issubset(set(out_dsNames))
                and out_shape == in_shape
                and out_size > in_size * 0.3):
            print('All datasets exists in file {} with same size as required,'
                  ' no need to re-load.'.format(os.path.basename(out_file)))
            flag = 'skip'

    return flag


def read_subset_box(template_file, meta):
    """Read subset info from template file

    Parameters: template_file - str, path of template file
                meta          - dict, metadata
    Returns:    pix_box       - tuple of 4 int in (x0, y0, x1, y1)
                meta          - dict, metadata
    """

    if template_file and os.path.isfile(template_file):

        # read subset info from template file
        pix_box, geo_box = read_subset_template2box(template_file)

        # geo_box --> pix_box
        if geo_box is not None:
            coord = ut.coordinate(meta)
            pix_box = coord.bbox_geo2radar(geo_box)
            pix_box = coord.check_box_within_data_coverage(pix_box)
            print(f'input bounding box in lalo: {geo_box}')

    else:
        pix_box = None

    if pix_box is not None:
        # update metadata against the new bounding box
        print(f'input bounding box in yx: {pix_box}')
        meta = attr.update_attribute4subset(meta, pix_box)
    else:
        # translate box of None to tuple of 4 int
        length, width = int(meta['LENGTH']), int(meta['WIDTH'])
        pix_box = (0, 0, width, length)

    # ensure all index are in int16
    pix_box = tuple(int(i) for i in pix_box)

    return pix_box, meta


####################################################################################
def extract_metadata(stack):
    """Extract ARIA metadata for MintPy."""

    meta = {}
    ds = gdal.Open(stack, gdal.GA_ReadOnly)
    if not ds:
        raise RuntimeError(f'Failed to open file {stack} with GDAL.')

    # read metadata from unwrapStack.vrt file
    print(f'extract metadata from {stack}')
    metaUnw = ds.GetRasterBand(1).GetMetadata("unwrappedPhase")

    # copy over all metadata from unwrapStack
    for key, value in metaUnw.items():
        if key not in ["Dates", "perpendicularBaseline"]:
            meta[key] = value

    meta["ANTENNA_SIDE"] = -1
    meta["PROCESSOR"] = "isce"
    meta["FILE_LENGTH"] = ds.RasterYSize
    meta["LENGTH"] = ds.RasterYSize
    meta["ORBIT_DIRECTION"] = meta["orbitDirection"].upper()
    meta["PLATFORM"] = "Sen"
    meta["WAVELENGTH"] = float(meta["Wavelength (m)"])
    meta["WIDTH"] = ds.RasterXSize
    meta["NUMBER_OF_PAIRS"] = ds.RasterCount
    meta["STARTING_RANGE"] = float(meta["startRange"])

    # Note from YZ, 2019-07-25
    # convert isce azimuth angle to roipac orbit heading angle
    # This value is not consistent with band2 of los.rdr from ISCE/topsStack
    # need to check with ARIA-tools team.
    # use hardwired value for now
    #az_angle = float(meta["azimuthAngle"])
    #head_angle = -1 * (270 + az_angle)
    #head_angle -= np.round(head_angle / 360.) * 360.
    #meta['HEADING'] = head_angle
    if meta["ORBIT_DIRECTION"].startswith("D"):
        meta["HEADING"] = -168
    else:
        meta["HEADING"] = -12

    # ARIA standard products currently don't have number of range and
    # azimuth looks. They are however fixed to the following values
    meta["ALOOKS"] = 7
    meta["RLOOKS"] = 19
    meta["RANGE_PIXEL_SIZE"] = float(meta["slantRangeSpacing"]) * meta["RLOOKS"]

    # number of independent looks
    sen_dict = sensor.SENSOR_DICT['sen']
    rgfact = sen_dict['IW2']['range_resolution'] / sen_dict['range_pixel_size']
    azfact = sen_dict['IW2']['azimuth_resolution'] / sen_dict['azimuth_pixel_size']
    meta['NCORRLOOKS'] = meta['RLOOKS'] * meta['ALOOKS'] / (rgfact * azfact)

    # geo transformation
    transform = ds.GetGeoTransform()

    x_step = abs(transform[1])
    y_step = abs(transform[5]) * -1.

    W = transform[0]
    N = transform[3]
    E = W + x_step * ds.RasterXSize
    S = N + y_step * ds.RasterYSize

    meta["X_FIRST"] = f'{W:.9f}'
    meta["Y_FIRST"] = f'{N:.9f}'
    meta["X_STEP"] = f'{x_step:.9f}'
    meta["Y_STEP"] = f'{y_step:.9f}'
    meta["X_UNIT"] = "degrees"
    meta["Y_UNIT"] = "degrees"

    utc = meta["UTCTime (HH:MM:SS.ss)"]
    utc = time.strptime(utc, "%H:%M:%S.%f")
    meta["CENTER_LINE_UTC"] = utc.tm_hour*3600.0 + utc.tm_min*60.0 + utc.tm_sec

    # following values probably won't be used anywhere for the geocoded data
    # earth radius
    meta["EARTH_RADIUS"] = 6337286.638938101
    # nominal altitude of Sentinel1 orbit
    meta["HEIGHT"] = 693000.0

    if meta["ORBIT_DIRECTION"].startswith("ASC"):
        meta["LAT_REF1"] = str(S)
        meta["LAT_REF2"] = str(S)
        meta["LAT_REF3"] = str(N)
        meta["LAT_REF4"] = str(N)
        meta["LON_REF1"] = str(W)
        meta["LON_REF2"] = str(E)
        meta["LON_REF3"] = str(W)
        meta["LON_REF4"] = str(E)
    else:
        meta["LAT_REF1"] = str(N)
        meta["LAT_REF2"] = str(N)
        meta["LAT_REF3"] = str(S)
        meta["LAT_REF4"] = str(S)
        meta["LON_REF1"] = str(E)
        meta["LON_REF2"] = str(W)
        meta["LON_REF3"] = str(E)
        meta["LON_REF4"] = str(W)

    ds = None
    return meta


def write_geometry(outfile, demFile, incAngleFile, azAngleFile=None, waterMaskFile=None,
                   box=None, xstep=1, ystep=1):
    """Write geometry HDF5 file from list of VRT files."""

    print('-'*50)
    # box to gdal arguments
    # link: https://gdal.org/python/osgeo.gdal.Band-class.html#ReadAsArray
    if box is not None:
        kwargs = dict(
            xoff=box[0],
            yoff=box[1],
            win_xsize=box[2]-box[0],
            win_ysize=box[3]-box[1])
    else:
        kwargs = dict()

    print(f'writing data to HDF5 file {outfile} with a mode ...')
    with h5py.File(outfile, 'a') as f:

        # height
        ds = gdal.Open(demFile, gdal.GA_ReadOnly)
        bnd = ds.GetRasterBand(1)
        data = np.array(bnd.ReadAsArray(**kwargs), dtype=np.float32)
        data = multilook_data(data, ystep, xstep, method='nearest')
        data[data == bnd.GetNoDataValue()] = np.nan
        f['height'][:,:] = data

        # slantRangeDistance
        f['slantRangeDistance'][:,:] = float(f.attrs['STARTING_RANGE'])

        # incidenceAngle
        ds = gdal.Open(incAngleFile, gdal.GA_ReadOnly)
        bnd = ds.GetRasterBand(1)
        data = bnd.ReadAsArray(**kwargs)
        data = multilook_data(data, ystep, xstep, method='nearest')
        data[data == bnd.GetNoDataValue()] = np.nan
        f['incidenceAngle'][:,:] = data

        # azimuthAngle
        if azAngleFile is not None:
            ds = gdal.Open(azAngleFile, gdal.GA_ReadOnly)
            bnd = ds.GetRasterBand(1)
            data = bnd.ReadAsArray(**kwargs)
            data = multilook_data(data, ystep, xstep, method='nearest')
            data[data == bnd.GetNoDataValue()] = np.nan
            # azimuth angle of the line-of-sight vector:
            # ARIA: vector from target to sensor measured from the east  in counterclockwise direction
            # ISCE: vector from sensor to target measured from the north in counterclockwise direction
            # convert ARIA format to ISCE format, which is used in mintpy
            data -= 90
            f['azimuthAngle'][:,:] = data

        # waterMask
        if waterMaskFile is not None:
            # read
            ds = gdal.Open(waterMaskFile, gdal.GA_ReadOnly)
            bnd = ds.GetRasterBand(1)
            water_mask = bnd.ReadAsArray(**kwargs)
            water_mask = multilook_data(water_mask, ystep, xstep, method='nearest')
            water_mask[water_mask == bnd.GetNoDataValue()] = False

            # assign False to invalid pixels based on incAngle data
            ds = gdal.Open(incAngleFile, gdal.GA_ReadOnly)
            bnd = ds.GetRasterBand(1)
            data = bnd.ReadAsArray(**kwargs)
            data = multilook_data(data, ystep, xstep, method='nearest')
            water_mask[data == bnd.GetNoDataValue()] = False

            # write
            f['waterMask'][:,:] = water_mask

    print(f'finished writing to HD5 file: {outfile}')
    return outfile


def write_ifgram_stack(outfile, unwStack, cohStack, connCompStack, ampStack=None,
                       box=None, xstep=1, ystep=1, mli_method='nearest'):
    """Write ifgramStack HDF5 file from stack VRT files
    """

    print('-'*50)
    stackFiles = [unwStack, cohStack, connCompStack, ampStack]
    max_digit = max(len(os.path.basename(str(i))) for i in stackFiles)
    for stackFile in stackFiles:
        if stackFile is not None:
            print('open {f:<{w}} with gdal ...'.format(f=os.path.basename(stackFile), w=max_digit))

    dsUnw = gdal.Open(unwStack, gdal.GA_ReadOnly)
    dsCoh = gdal.Open(cohStack, gdal.GA_ReadOnly)
    dsComp = gdal.Open(connCompStack, gdal.GA_ReadOnly)
    if ampStack is not None:
        dsAmp = gdal.Open(ampStack, gdal.GA_ReadOnly)
    else:
        dsAmp = None

    # extract NoDataValue (from the last */date2_date1.vrt file for example)
    ds = gdal.Open(dsUnw.GetFileList()[-1], gdal.GA_ReadOnly)
    noDataValueUnw = ds.GetRasterBand(1).GetNoDataValue()
    print(f'grab NoDataValue for unwrapPhase     : {noDataValueUnw:<5} and convert to 0.')

    ds = gdal.Open(dsCoh.GetFileList()[-1], gdal.GA_ReadOnly)
    noDataValueCoh = ds.GetRasterBand(1).GetNoDataValue()
    print(f'grab NoDataValue for coherence       : {noDataValueCoh:<5} and convert to 0.')

    ds = gdal.Open(dsComp.GetFileList()[-1], gdal.GA_ReadOnly)
    noDataValueComp = ds.GetRasterBand(1).GetNoDataValue()
    print(f'grab NoDataValue for connectComponent: {noDataValueComp:<5} and convert to 0.')
    ds = None

    if dsAmp is not None:
        ds = gdal.Open(dsAmp.GetFileList()[-1], gdal.GA_ReadOnly)
        noDataValueAmp = ds.GetRasterBand(1).GetNoDataValue()
        print(f'grab NoDataValue for magnitude       : {noDataValueAmp:<5} and convert to 0.')
        ds = None

    # sort the order of interferograms based on date1_date2 with date1 < date2
    nPairs = dsUnw.RasterCount
    d12BandDict = {}
    for ii in range(nPairs):
        bnd = dsUnw.GetRasterBand(ii+1)
        d12 = bnd.GetMetadata("unwrappedPhase")["Dates"]
        d12 = sorted(d12.split("_"))
        d12 = f'{d12[0]}_{d12[1]}'
        d12BandDict[d12] = ii+1
    d12List = sorted(d12BandDict.keys())
    print(f'number of interferograms: {len(d12List)}')

    # box to gdal arguments
    # link: https://gdal.org/python/osgeo.gdal.Band-class.html#ReadAsArray
    if box is not None:
        kwargs = dict(
            xoff=box[0],
            yoff=box[1],
            win_xsize=box[2]-box[0],
            win_ysize=box[3]-box[1])
    else:
        kwargs = dict()

    if xstep * ystep > 1:
        msg = f'apply {xstep} x {ystep} multilooking/downsampling via {mli_method} to: unwrapPhase, coherence'
        msg += ', magnitude' if dsAmp is not None else ''
        msg += f'\napply {xstep} x {ystep} multilooking/downsampling via nearest to: connectComponent'
        print(msg)
    print(f'writing data to HDF5 file {outfile} with a mode ...')
    with h5py.File(outfile, "a") as f:

        prog_bar = ptime.progressBar(maxValue=nPairs)
        for ii in range(nPairs):
            d12 = d12List[ii]
            bndIdx = d12BandDict[d12]
            prog_bar.update(ii+1, suffix=f'{d12} {ii+1}/{nPairs}')

            f["date"][ii,0] = d12.split("_")[0].encode("utf-8")
            f["date"][ii,1] = d12.split("_")[1].encode("utf-8")
            f["dropIfgram"][ii] = True

            bnd = dsUnw.GetRasterBand(bndIdx)
            data = bnd.ReadAsArray(**kwargs)
            data = multilook_data(data, ystep, xstep, method=mli_method)
            data[data == noDataValueUnw] = 0      #assign pixel with no-data to 0
            data *= -1.0                          #date2_date1 -> date1_date2
            f["unwrapPhase"][ii,:,:] = data

            bperp = float(bnd.GetMetadata("unwrappedPhase")["perpendicularBaseline"])
            bperp *= -1.0                         #date2_date1 -> date1_date2
            f["bperp"][ii] = bperp

            bnd = dsCoh.GetRasterBand(bndIdx)
            data = bnd.ReadAsArray(**kwargs)
            data = multilook_data(data, ystep, xstep, method=mli_method)
            data[data == noDataValueCoh] = 0      #assign pixel with no-data to 0
            f["coherence"][ii,:,:] = data

            bnd = dsComp.GetRasterBand(bndIdx)
            data = bnd.ReadAsArray(**kwargs)
            data = multilook_data(data, ystep, xstep, method='nearest')
            data[data == noDataValueComp] = 0     #assign pixel with no-data to 0
            f["connectComponent"][ii,:,:] = data

            if dsAmp is not None:
                bnd = dsAmp.GetRasterBand(bndIdx)
                data = bnd.ReadAsArray(**kwargs)
                data = multilook_data(data, ystep, xstep, method=mli_method)
                data[data == noDataValueAmp] = 0  #assign pixel with no-data to 0
                f["magnitude"][ii,:,:] = data

        prog_bar.close()

        # add MODIFICATION_TIME metadata to each 3D dataset
        for dsName in ['unwrapPhase','coherence','connectComponent']:
            f[dsName].attrs['MODIFICATION_TIME'] = str(time.time())

    print(f'finished writing to HD5 file: {outfile}')
    dsUnw = None
    dsCoh = None
    dsComp = None
    dsAmp = None
    return outfile


####################################################################################
def load_aria(inps):
    """Prepare and load ARIA data and metadata into HDF5/MintPy format."""

    start_time = time.time()
    print(f'update mode: {inps.updateMode}')

    # extract metadata
    meta = extract_metadata(inps.unwFile)
    box, meta = read_subset_box(inps.template_file, meta)
    if inps.xstep * inps.ystep > 1:
        meta = attr.update_attribute4multilook(
            meta,
            lks_y=inps.ystep,
            lks_x=inps.xstep,
        )

    length = int(meta["LENGTH"])
    width = int(meta["WIDTH"])
    num_pair = int(meta["NUMBER_OF_PAIRS"])

    # prepare output directory
    out_dir = os.path.dirname(inps.outfile[0])
    os.makedirs(out_dir, exist_ok=True)

    ########## output file 1 - ifgramStack
    # define dataset structure for ifgramStack
    ds_name_dict = {
        "date"             : (np.dtype('S8'), (num_pair, 2)),
        "dropIfgram"       : (np.bool_,       (num_pair,)),
        "bperp"            : (np.float32,     (num_pair,)),
        "unwrapPhase"      : (np.float32,     (num_pair, length, width)),
        "coherence"        : (np.float32,     (num_pair, length, width)),
        "connectComponent" : (np.int16,       (num_pair, length, width)),
    }
    if inps.magFile is not None:
        ds_name_dict['magnitude'] = (np.float32, (num_pair, length, width))

    if run_or_skip(inps, ds_name_dict, out_file=inps.outfile[0]) == 'run':
        # initiate h5 file with defined structure
        meta['FILE_TYPE'] = 'ifgramStack'
        writefile.layout_hdf5(
            inps.outfile[0],
            ds_name_dict,
            metadata=meta,
            compression=inps.compression,
        )

        # write data to h5 file in disk
        write_ifgram_stack(
            inps.outfile[0],
            unwStack=inps.unwFile,
            cohStack=inps.corFile,
            connCompStack=inps.connCompFile,
            ampStack=inps.magFile,
            box=box,
            xstep=inps.xstep,
            ystep=inps.ystep,
            mli_method=inps.method,
        )

    ########## output file 2 - geometryGeo
    # define dataset structure for geometry
    ds_name_dict = {
        "height"             : (np.float32, (length, width)),
        "incidenceAngle"     : (np.float32, (length, width)),
        "slantRangeDistance" : (np.float32, (length, width)),
    }
    if inps.azAngleFile is not None:
        ds_name_dict["azimuthAngle"] = (np.float32, (length, width))
    if inps.waterMaskFile is not None:
        ds_name_dict["waterMask"] = (np.bool_, (length, width))

    if run_or_skip(inps, ds_name_dict, out_file=inps.outfile[1]) == 'run':
        # initiate h5 file with defined structure
        meta['FILE_TYPE'] = 'geometry'
        writefile.layout_hdf5(
            inps.outfile[1],
            ds_name_dict,
            metadata=meta,
            compression=inps.compression,
        )

        # write data to disk
        write_geometry(
            inps.outfile[1],
            demFile=inps.demFile,
            incAngleFile=inps.incAngleFile,
            azAngleFile=inps.azAngleFile,
            waterMaskFile=inps.waterMaskFile,
            box=box,
            xstep=inps.xstep,
            ystep=inps.ystep,
        )

    print('-'*50)

    # used time
    m, s = divmod(time.time() - start_time, 60)
    print(f'time used: {m:02.0f} mins {s:02.1f} secs.')

    return
