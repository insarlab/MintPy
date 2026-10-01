############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, Zhang Yunjun, Emre Havazli, 2019 #
############################################################

import datetime as dt
import os
import time

import h5py
import numpy as np
from osgeo import gdal

gdal.UseExceptions()

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
    if not metaUnw:
        domains = ds.GetRasterBand(1).GetMetadataDomainList()
        metaUnw = ds.GetRasterBand(1).GetMetadata(domains[0]) if domains else {}

    # copy over all metadata from unwrapStack
    for key, value in metaUnw.items():
        if key not in ["Dates", "perpendicularBaseline"]:
            meta[key] = value

    meta["ANTENNA_SIDE"] = -1
    meta["PROCESSOR"] = "isce"
    meta["FILE_LENGTH"] = ds.RasterYSize
    meta["LENGTH"] = ds.RasterYSize
    meta["ORBIT_DIRECTION"] = meta.get("orbitDirection", "UNKNOWN").upper()
    meta["PLATFORM"] = meta.get("PLATFORM", "Sen")
    meta["WAVELENGTH"] = float(meta.get("Wavelength (m)", 0.05546576))
    meta["WIDTH"] = ds.RasterXSize
    meta["NUMBER_OF_PAIRS"] = ds.RasterCount
    meta["STARTING_RANGE"] = float(meta.get("startRange", 0.0))

    if meta["ORBIT_DIRECTION"].startswith("D"):
        meta["HEADING"] = -168
    else:
        meta["HEADING"] = -12

    # Set number of looks and correlation looks based on the platform
    if meta["PLATFORM"].upper().startswith("NISAR"):
        meta["ALOOKS"] = 1
        meta["RLOOKS"] = 1
        meta["RANGE_PIXEL_SIZE"] = float(meta.get("slantRangeSpacing", 1.0))
        meta['NCORRLOOKS'] = 1.0 
    else:
        meta["ALOOKS"] = 7
        meta["RLOOKS"] = 19
        meta["RANGE_PIXEL_SIZE"] = float(meta.get("slantRangeSpacing", 2.33)) * meta["RLOOKS"]

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

    utc = meta.get("UTCTime (HH:MM:SS.ss)", "00:00:00.00")
    try:
        utc_obj = time.strptime(utc, "%H:%M:%S.%f")
        meta["CENTER_LINE_UTC"] = utc_obj.tm_hour*3600.0 + utc_obj.tm_min*60.0 + utc_obj.tm_sec
    except ValueError:
        meta["CENTER_LINE_UTC"] = 0.0

    meta["EARTH_RADIUS"] = 6337286.638938101
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
            data -= 90
            f['azimuthAngle'][:,:] = data

        # waterMask
        if waterMaskFile is not None:
            ds = gdal.Open(waterMaskFile, gdal.GA_ReadOnly)
            bnd = ds.GetRasterBand(1)
            water_mask = bnd.ReadAsArray(**kwargs)
            water_mask = multilook_data(water_mask, ystep, xstep, method='nearest')
            water_mask[water_mask == bnd.GetNoDataValue()] = False

            ds = gdal.Open(incAngleFile, gdal.GA_ReadOnly)
            bnd = ds.GetRasterBand(1)
            data = bnd.ReadAsArray(**kwargs)
            data = multilook_data(data, ystep, xstep, method='nearest')
            water_mask[data == bnd.GetNoDataValue()] = False

            f['waterMask'][:,:] = water_mask

            if azAngleFile is not None:
                f['azimuthAngle'][:,:] *= water_mask

    print(f'finished writing to HD5 file: {outfile}\n')
    return outfile


def write_ifgram_stack(outfile, stackFiles, box=None, xstep=1, ystep=1, mli_method='nearest'):
    """Write stacks to HDF5 files from stack VRT files
    """

    print('-'*50)

    # remove None entries
    stackFiles = {key:val for key, val in stackFiles.items() if val is not None}

    # check all files exist
    for dsName, fname in stackFiles.items():
        if not os.path.exists(fname):
            raise Exception("%s does not exist" % fname)

    # determine field length for printing
    max_digit = max(len(os.path.basename(str(i))) for i in stackFiles.values())
    for stackFile in stackFiles.values():
        if stackFile is not None:
            print('open {f:<{w}} with gdal ...'.format(f=os.path.basename(stackFile), w=max_digit))

    # extract NoDataValue for each stack
    noDataValues = {}
    for dsName in stackFiles.keys():
        dsStack = gdal.Open(stackFiles[dsName], gdal.GA_ReadOnly)
        ds = gdal.Open(dsStack.GetFileList()[-1], gdal.GA_ReadOnly)
        noDataValues[dsName] = ds.GetRasterBand(1).GetNoDataValue()

        fileName = os.path.basename(stackFiles[dsName])
        nodata_str = str(noDataValues[dsName])
        print(f'grab NoDataValue for {fileName:<{max_digit}}: '
              f'{nodata_str:<5} and convert to 0.')
        ds = None

    # sort the order of interferograms based on date1_date2 with date1 < date2
    nPairs = dsStack.RasterCount
    d12BandDict = {}
    for ii in range(nPairs):
        bnd = dsStack.GetRasterBand(ii+1)
        domains = bnd.GetMetadataDomainList()
        meta = bnd.GetMetadata(domains[0]) if domains else {}
        d12 = meta.get("Dates", f"pair_{ii}")
        d12 = sorted(d12.split("_"))
        d12 = f'{d12[0]}_{d12[1]}' if len(d12) > 1 else d12[0]
        d12BandDict[d12] = ii+1
    d12List = sorted(d12BandDict.keys())
    print(f'number of interferograms: {len(d12List)}')

    if box is not None:
        kwargs = dict(
            xoff=box[0],
            yoff=box[1],
            win_xsize=box[2]-box[0],
            win_ysize=box[3]-box[1])
    else:
        kwargs = dict()

    # write to HDF5 file
    print(f'writing data to HDF5 file {outfile} with a mode ...')
    with h5py.File(outfile, "a") as f:
        prog_bar = ptime.progressBar(maxValue=nPairs)
        for ii in range(nPairs):
            d12 = d12List[ii]
            bndIdx = d12BandDict[d12]
            prog_bar.update(ii+1, suffix=f'{d12} {ii+1}/{nPairs}')

            dates_parts = d12.split("_")
            f["date"][ii,0] = dates_parts[0].encode("utf-8")
            f["date"][ii,1] = dates_parts[1].encode("utf-8") if len(dates_parts) > 1 else dates_parts[0].encode("utf-8")
            f["dropIfgram"][ii] = True

            # loop through stacks
            for dsName in stackFiles.keys():
                dsStack = gdal.Open(stackFiles[dsName], gdal.GA_ReadOnly)
                bnd = dsStack.GetRasterBand(bndIdx)
                data = bnd.ReadAsArray(**kwargs)
                if xstep * ystep > 1:
                    mli_method_spec = mli_method if dsName not in \
                        ['connCompStack'] else 'nearest'
                    data = multilook_data(data, ystep, xstep, method=mli_method_spec)

                nodata = noDataValues[dsName]
                if nodata is not None:
                    if np.isnan(nodata):
                        data[np.isnan(data)] = 0
                    else:
                        data[data == nodata] = 0
                data[np.isnan(data)] = 0

                domains = bnd.GetMetadataDomainList()
                bnd_meta = bnd.GetMetadata(domains[0]) if domains else {}

                if dsName == 'unwrapPhase':
                    data *= -1  # date2_date1 -> date1_date2
                    f['unwrapPhase'][ii,:,:] = data

                    bperp = float(bnd_meta.get("perpendicularBaseline", 0.0))
                    bperp *= -1.0  # date2_date1 -> date1_date2
                    f["bperp"][ii] = bperp

                elif dsName == 'coherence':
                    f["coherence"][ii,:,:] = data

                elif dsName == 'connectComponent':
                    f["connectComponent"][ii,:,:] = data

                elif dsName == 'magnitude':
                    f["magnitude"][ii,:,:] = data

                elif dsName == 'ionosphere':
                    platform = f.attrs.get("PLATFORM", "Sen")
                    if hasattr(platform, "decode"):
                        platform = platform.decode("utf-8")
                    is_nisar = str(platform).upper().startswith("NISAR")
                    
                    data *= -1.0
                    f["unwrapPhase"][ii,:,:] = data

                    bperp = float(bnd_meta.get("perpendicularBaseline", 0.0))
                    if not is_nisar:
                        bperp *= -1.0
                    f["bperp"][ii] = bperp

        prog_bar.close()

        # add MODIFICATION_TIME metadata to each 3D dataset
        for dsName in stackFiles.keys():
            dsName = 'unwrapPhase' if dsName == 'ionosphere' else dsName
            f[dsName].attrs['MODIFICATION_TIME'] = str(time.time())

    print(f'finished writing to HD5 file: {outfile}\n')
    return outfile


# OPTIONAL - ARIA model-based corrections troposphereTotal, solidearthtides
def write_timeseries(outfile, corrStack, box=None,
                      xstep=1, ystep=1, mli_method='nearest'):
    """Write per-epoch SET and tropospheric phase corrections in meters.

    ARIA correction layers use the same phase convention as the unwrapped
    interferograms.  MintPy displacement is therefore ``-phase * lambda / 4pi``.
    The resulting time series can be subtracted directly by ``diff.py``.
    """

    print('-'*50)

    max_digit = len(os.path.basename(str(corrStack)))

    if corrStack is not None:
        print('open {f:<{w}} with gdal ...'.format(f=os.path.basename(corrStack), w=max_digit))

        if not os.path.exists(corrStack):
            raise Exception("%s does not exist" % corrStack)

    dsCor = gdal.Open(corrStack, gdal.GA_ReadOnly)
    ds = gdal.Open(dsCor.GetFileList()[-1], gdal.GA_ReadOnly)
    noDataValue = ds.GetRasterBand(1).GetNoDataValue()
    ds = None

    domains = dsCor.GetRasterBand(1).GetMetadataDomainList()
    layer = domains[0] if domains else "UNKNOWN"

    bnd_meta = dsCor.GetRasterBand(1).GetMetadata(layer)
    wavelength = np.float64(bnd_meta.get("Wavelength (m)", 0.05546576))
    phase2range = -wavelength / (4.0 * np.pi)

    nDate = dsCor.RasterCount
    dateDict = {}
    sensingDict = {}
    for ii in range(nDate):
        bnd = dsCor.GetRasterBand(ii+1)
        meta = bnd.GetMetadata(layer)
        date = meta.get("Dates", f"date_{ii}")
        utc_str = meta.get("UTCTime (HH:MM:SS.ss)", "00:00:00.00")
        try:
            utc = dt.datetime.strptime(date + ',' + utc_str, "%Y%m%d,%H:%M:%S.%f")
        except ValueError:
            utc = f"{date} 00:00:00"
        dateDict[date] = ii+1
        sensingDict[ii+1] = str(utc)
    dateList = sorted(dateDict.keys())
    print(f'number of {layer} datasets: {len(dateList)}')

    if box is not None:
        kwargs = dict(
            xoff=box[0],
            yoff=box[1],
            win_xsize=box[2]-box[0],
            win_ysize=box[3]-box[1])
    else:
        kwargs = dict()

    if xstep * ystep > 1:
        msg = f'apply {xstep} x {ystep} multilooking/downsampling via {mli_method} to {layer}'
        print(msg)

    print(f'writing data to HDF5 file {outfile} with a mode ...')
    with h5py.File(outfile, "a") as f:
        prog_bar = ptime.progressBar(maxValue=nDate)
        for ii in range(nDate):
            date = dateList[ii]
            bndIdx = dateDict[date]
            utc = sensingDict[bndIdx]
            prog_bar.update(ii+1, suffix=f'{date} {ii+1}/{nDate}')

            f["date"][ii] = date.encode("utf-8")
            f["sensingMid"][ii] = utc.encode("utf-8")

            bnd = dsCor.GetRasterBand(bndIdx)
            data = bnd.ReadAsArray(**kwargs)
            data = multilook_data(data, ystep, xstep, method=mli_method)
            if noDataValue is not None:
                data[data == noDataValue] = 0
            data[np.isnan(data)] = 0
            f["timeseries"][ii,:,:] = data * phase2range

        prog_bar.close()

        for dsName in ['timeseries']:
            f[dsName].attrs['MODIFICATION_TIME'] = str(time.time())

    print(f'finished writing to HD5 file: {outfile}\n')
    dsCor = None

    return outfile


def get_nisar_dates(corr_stack):
    """Extract unique acquisition dates and UTC sensing times from NISAR VRT."""
    ds_cor = gdal.Open(corr_stack, gdal.GA_ReadOnly)
    domains = ds_cor.GetRasterBand(1).GetMetadataDomainList()
    layer = domains[0] if domains else "UNKNOWN"
    n_pairs = ds_cor.RasterCount

    unique_dates = set()
    date_utc_dict = {}

    for ii in range(n_pairs):
        bnd = ds_cor.GetRasterBand(ii + 1)
        meta_dict = bnd.GetMetadata(layer)
        dates = meta_dict.get("Dates", "00000000_00000000").split("_")
        utc_str = meta_dict.get("UTCTime (HH:MM:SS.ss)", "00:00:00.000000")

        for d in dates:
            unique_dates.add(d)
            if d not in date_utc_dict:
                try:
                    utc = dt.datetime.strptime(f"{d},{utc_str}", "%Y%m%d,%H:%M:%S.%f")
                    date_utc_dict[d] = str(utc)
                except ValueError:
                    date_utc_dict[d] = f"{d} 00:00:00"

    date_list = sorted(list(unique_dates))
    ds_cor = None
    return date_list, date_utc_dict


def invert_diff_corrections(outfile, corrStack, box=None, xstep=1, ystep=1, mli_method='nearest'):
    """Invert NISAR pairwise phase corrections into a meter time series.

    ARIA-tools normalizes NISAR external corrections to its ``date2_date1``
    phase convention.  Use the same phase-to-range sign as the Sentinel
    per-epoch path so the output is ready for subtraction by MintPy ``diff.py``.
    """
    print('-' * 50)
    max_digit = len(os.path.basename(str(corrStack)))

    if corrStack is not None:
        print(f'open {os.path.basename(corrStack):<{max_digit}} with gdal ...')
        if not os.path.exists(corrStack):
            raise FileNotFoundError(f"{corrStack} does not exist")

    ds_cor = gdal.Open(corrStack, gdal.GA_ReadOnly)
    ds = gdal.Open(ds_cor.GetFileList()[-1], gdal.GA_ReadOnly)
    no_data_val = ds.GetRasterBand(1).GetNoDataValue()
    ds = None

    layer_name, _ = get_correction_layer(corrStack)
    bnd_meta = ds_cor.GetRasterBand(1).GetMetadata(layer_name)
    wavelength = np.float64(bnd_meta.get("Wavelength (m)", 0.05546576))
    phase2range = -wavelength / (4.0 * np.pi)

    n_pairs = ds_cor.RasterCount
    date_list, date_utc_dict = get_nisar_dates(corrStack)
    n_dates = len(date_list)
    date2idx = {d: i for i, d in enumerate(date_list)}

    print(f'number of differential pairs: {n_pairs}')
    print(f'number of unique {layer_name} dates: {n_dates}')

    A = np.zeros((n_pairs, n_dates), dtype=np.float32)

    kwargs = dict()
    if box is not None:
        kwargs = dict(
            xoff=box[0],
            yoff=box[1],
            win_xsize=box[2] - box[0],
            win_ysize=box[3] - box[1],
        )

    bnd1 = ds_cor.GetRasterBand(1)
    sample_data = bnd1.ReadAsArray(**kwargs)
    if xstep * ystep > 1:
        sample_data = multilook_data(sample_data, ystep, xstep, method=mli_method)
    length, width = sample_data.shape
    num_pixels = length * width

    diff_data = np.zeros((n_pairs, length, width), dtype=np.float32)

    for ii in range(n_pairs):
        bnd = ds_cor.GetRasterBand(ii + 1)
        data = bnd.ReadAsArray(**kwargs)
        if xstep * ystep > 1:
            data = multilook_data(data, ystep, xstep, method=mli_method)

        d12 = bnd.GetMetadata(layer_name).get("Dates", "00000000_00000000")
        d_parts = d12.split("_")
        date2 = d_parts[0]
        date1 = d_parts[1] if len(d_parts) > 1 else d_parts[0]

        idx_date2 = date2idx[date2]
        idx_date1 = date2idx[date1]

        A[ii, idx_date2] = 1.0
        A[ii, idx_date1] = -1.0

        if no_data_val is not None:
            data[data == no_data_val] = np.nan
        diff_data[ii, :, :] = data * phase2range

    A_sub = A[:, 1:]
    A_pinv = np.linalg.pinv(A_sub)

    diff_data_flat = diff_data.reshape(n_pairs, num_pixels)
    nan_mask = np.isnan(diff_data_flat)
    diff_data_flat[nan_mask] = 0.0

    ts_sub = np.matmul(A_pinv, diff_data_flat)

    ts_flat = np.zeros((n_dates, num_pixels), dtype=np.float32)
    ts_flat[1:, :] = ts_sub

    has_nan_pixel = np.any(nan_mask, axis=0)
    ts_flat[:, has_nan_pixel] = 0.0

    timeseries = ts_flat.reshape(n_dates, length, width)

    if xstep * ystep > 1:
        print(f'apply {xstep} x {ystep} multilooking/downsampling via {mli_method} to {layer_name}')

    print(f'writing data to HDF5 file {outfile} with a mode ...')
    with h5py.File(outfile, "a") as f:
        prog_bar = ptime.progressBar(maxValue=n_dates)
        for ii in range(n_dates):
            date = date_list[ii]
            utc = date_utc_dict.get(date, f"{date} 00:00:00")
            prog_bar.update(ii + 1, suffix=f'{date} {ii + 1}/{n_dates}')

            f["date"][ii] = date.encode("utf-8")
            f["sensingMid"][ii] = utc.encode("utf-8")
            f["timeseries"][ii, :, :] = timeseries[ii, :, :]

        prog_bar.close()

        f["timeseries"].attrs['MODIFICATION_TIME'] = str(time.time())

    print(f'finished writing to HD5 file: {outfile}\n')
    ds_cor = None

    return outfile


def get_number_of_epochs(vrtfile):
    ds = gdal.Open(vrtfile, gdal.GA_ReadOnly)
    count = ds.RasterCount
    ds = None
    return count


def get_correction_layer(correction_filename):
    ds = gdal.Open(correction_filename, gdal.GA_ReadOnly)
    domains = ds.GetRasterBand(1).GetMetadataDomainList()
    layer_name = domains[0] if domains else "UNKNOWN"

    tropo_models = [
        'GMAO', 'HRES', 'HRRR', 'ERA5',
        'troposphereTotal', 'troposphereWet', 'troposphereHydrostatic',
    ]
    if layer_name in tropo_models or 'troposphere' in layer_name.lower():
        layer_type = 'tropo'
    else:
        layer_type = layer_name

    ds = None
    return layer_name, layer_type

####################################################################################
def load_aria(inps):
    """Prepare and load ARIA data and metadata into HDF5/MintPy format."""

    start_time = time.time()
    print(f'update mode: {inps.updateMode}')

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

    out_dir = os.path.dirname(inps.outfile[0])
    os.makedirs(out_dir, exist_ok=True)

    ########## output file 1 - ifgramStack
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
        meta['FILE_TYPE'] = 'ifgramStack'
        writefile.layout_hdf5(
            inps.outfile[0],
            ds_name_dict,
            metadata=meta,
            compression=None if inps.compression == 'default' else inps.compression,
        )

        write_ifgram_stack(
            inps.outfile[0],
            stackFiles={'unwrapPhase': inps.unwFile,
                         'coherence': inps.corFile,
                         'connectComponent': inps.connCompFile,
                         'magnitude': inps.magFile},
            box=box,
            xstep=inps.xstep,
            ystep=inps.ystep,
            mli_method=inps.method,
        )

    ########## output file 2 - geometryGeo
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
        meta['FILE_TYPE'] = 'geometry'
        writefile.layout_hdf5(
            inps.outfile[1],
            ds_name_dict,
            metadata=meta,
            compression='lzf' if inps.compression == 'default' else inps.compression,
        )

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

    ########## output file 3 - correction layers

    # 3.1 - ionosphere
    if inps.ionoFile:
        ds_name_dict = {
            'date'             : (np.dtype('S8'), (num_pair, 2)),
            'dropIfgram'       : (np.bool_,       (num_pair,)),
            'bperp'            : (np.float32,     (num_pair,)),
            'unwrapPhase'      : (np.float32,     (num_pair, length, width)),
            "coherence"        : (np.float32,     (num_pair, length, width)),
        }
        meta['FILE_TYPE'] = 'ifgramStack'

        layer_name, _ = get_correction_layer(inps.ionoFile)

        outname = f'{out_dir}/ionStack.h5'
        if run_or_skip(inps, ds_name_dict, out_file=outname) == 'run':

            writefile.layout_hdf5(
                outname,
                ds_name_dict,
                metadata=meta,
                compression=None if inps.compression == 'default' else inps.compression,
                )

            write_ifgram_stack(
                outname,
                stackFiles={'ionosphere': inps.ionoFile,
                            'coherence': inps.corFile,},
                box=box,
                xstep=inps.xstep,
                ystep=inps.ystep,
                mli_method=inps.method,
            )

    # 3.2 - model based corrections: SolidEarthTides and Troposphere
    if inps.tropoFile is None:
        inps.tropoFile = [None]
    correction_layers = inps.tropoFile + [inps.setFile]
    for layer in correction_layers:
        if layer:
            layer_name, _ = get_correction_layer(layer)
            is_nisar = meta.get("PLATFORM", "").upper().startswith("NISAR")

            if is_nisar:
                date_list, _ = get_nisar_dates(layer)
                num_dates = len(date_list)
            else:
                num_dates = get_number_of_epochs(layer)

            meta['FILE_TYPE'] = 'timeseries'
            meta['UNIT'] = 'm'
            meta['DATA_TYPE'] = 'float32'
            for key in ['REF_Y', 'REF_X', 'REF_DATE']:
                if key in meta.keys():
                    meta.pop(key)

            ds_name_dict = {
                'date'       : (np.dtype('S8'),  (num_dates, )),
                'sensingMid' : (np.dtype('S15'), (num_dates, )),
                'timeseries' : (np.float32,      (num_dates, length, width)),
            }

            out_file = f'{out_dir}/{layer_name}_ARIA.h5'
            if run_or_skip(inps, ds_name_dict, out_file=out_file) == 'run':
                writefile.layout_hdf5(
                    out_file,
                    ds_name_dict,
                    metadata=meta,
                    compression=None if inps.compression == 'default' else inps.compression,
                )

                if is_nisar:
                    invert_diff_corrections(
                        out_file,
                        corrStack=layer,
                        box=box,
                        xstep=inps.xstep,
                        ystep=inps.ystep,
                        mli_method=inps.method,
                    )
                else:
                    write_timeseries(
                        out_file,
                        corrStack=layer,
                        box=box,
                        xstep=inps.xstep,
                        ystep=inps.ystep,
                    )
    print('-'*50)

    m, s = divmod(time.time() - start_time, 60)
    print(f'time used: {m:02.0f} mins {s:02.1f} secs.')
