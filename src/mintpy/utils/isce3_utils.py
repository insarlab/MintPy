import os
import re
import json
import tempfile
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional, Union

import h5py
import subprocess
import numpy as np
from osgeo import gdal, osr

from mintpy.utils import ptime, readfile
from mintpy.objects import sensor



def extract_isce3_metadata(meta_file: str, update_mode: bool = True) -> dict:
    """Extract common metadata from an ISCE3/Dolphin burst XML file.

    Parameters
    ----------
    meta_file : str
        Path to the reference burst XML file.
    update_mode : bool
        Not used here (kept for consistency).

    Returns
    -------
    dict
        Common metadata dictionary with keys required by MintPy.
    """
    # Parse XML file
    tree = ET.parse(meta_file)
    root = tree.getroot()
    burst_elem = root.find('burst_attributes')
    if burst_elem is None:
        raise ValueError(f'Missing <burst_attributes> in {meta_file}')

    # Helper to extract text from element
    def get_value(tag):
        elem = burst_elem.find(tag)
        if elem is not None:
            return elem.text.strip()
        return None

    meta = {}

    # Basic radar parameters
    meta['prf'] = get_value('prf')
    meta['startUTC'] = get_value('burstStartUTC')
    meta['stopUTC'] = get_value('burstStopUTC')
    meta['radarWavelength'] = get_value('radarWavelength')
    meta['startingRange'] = get_value('startingRange')
    meta['passDirection'] = get_value('passDirection')
    meta['polarization'] = get_value('polarization')
    meta['trackNumber'] = get_value('trackNumber')
    meta['orbitNumber'] = get_value('orbitNumber')

    # Platform name (Sentinel-1)
    meta['PLATFORM'] = 'sen'

    # Sensing mid time and center line UTC
    sensing_mid = get_value('sensingMid')
    if sensing_mid:
        # Parse time to compute seconds of day (optional)
        # For simplicity, we may skip CENTER_LINE_UTC calculation or use a placeholder
        # Here we set a default
        meta['CENTER_LINE_UTC'] = '0'
    else:
        meta['CENTER_LINE_UTC'] = '0'

    # Pixel sizes
    az_time_interval = get_value('azimuthTimeInterval')
    range_pixel_size = get_value('rangePixelSize')
    satellite_speed = get_value('satelliteSpeed')

    if az_time_interval and satellite_speed:
        try:
            az_pixel_size = float(satellite_speed) * float(az_time_interval)
            meta['azimuthPixelSize'] = str(az_pixel_size)
        except ValueError:
            meta['azimuthPixelSize'] = '0.0'
    else:
        meta['azimuthPixelSize'] = '0.0'

    if range_pixel_size:
        meta['rangePixelSize'] = range_pixel_size
    else:
        meta['rangePixelSize'] = '0.0'

    # Spatial resolution from sensor database (Sentinel-1 IW)
    # Determine swath (e.g., 'IW2')
    swath_num = get_value('swathNumber')
    if swath_num:
        iw_str = f'IW{swath_num}'
    else:
        # fallback: try to infer from filename
        base = os.path.basename(meta_file)
        if base.startswith('IW'):
            iw_str = base.split('.')[0]
        else:
            iw_str = 'IW2'
    try:
        meta['azimuthResolution'] = str(sensor.SENSOR_DICT['sen'][iw_str]['azimuth_resolution'])
        meta['rangeResolution']   = str(sensor.SENSOR_DICT['sen'][iw_str]['range_resolution'])
    except KeyError:
        meta['azimuthResolution'] = '0.0'
        meta['rangeResolution']   = '0.0'

    # Heading, earth radius, altitude
    meta['HEADING'] = get_value('HEADING')
    meta['earthRadius'] = get_value('earthRadius')
    meta['altitude'] = get_value('altitude')

    # Beam mode and swath
    meta['beam_mode'] = 'IW'
    meta['swathNumber'] = swath_num if swath_num else '2'

    # Frame numbers (may be 0 for ISCE3 products)
    meta['firstFrameNumber'] = get_value('firstFrameNumber') or '0'
    meta['lastFrameNumber'] = get_value('lastFrameNumber') or '0'

    # Antenna side (default to -1, right-looking)
    meta['ANTENNA_SIDE'] = '-1'

    # Processor
    meta['PROCESSOR'] = 'isce3'

    # Initial looks
    meta['ALOOKS'] = '1'
    meta['RLOOKS'] = '1'

    # Convert all values to strings (MintPy standard)
    for key, value in meta.items():
        if value is not None:
            meta[key] = str(value)
        else:
            meta[key] = ''

    # Standardize keys
    meta = readfile.standardize_metadata(meta)

    # Note: LENGTH and WIDTH are not in this XML; they will be added later
    # from geometry or interferogram files.

    return meta


def read_baseline_timeseries_isce3(baseline_dir: str, processor: str = 'tops') -> Dict:
    """Read baseline time series from ISCE3/Dolphin baseline directory.

    Expected structure: baseline_dir/*.txt where each filename is YYYYMMDD_YYYYMMDD.txt
    File content example:
        Bperp average (m): -19.081015753493432
        Bpar average (m): -6.823039248934357

    Parameters
    ----------
    baseline_dir : str
        Path to the baseline directory.
    processor : str
        Processor name (unused, kept for compatibility).

    Returns
    -------
    dict
        Dictionary of baseline values keyed by date (YYYYMMDD).
        Each value is [bperp_top, bperp_bottom] (identical for both).
    """
    import glob
    import os
    from mintpy.utils import ptime

    baseline_dict = {}
    
    # Find all .txt files matching the pattern YYYYMMDD_YYYYMMDD.txt
    pattern = os.path.join(baseline_dir, '[0-9]*_[0-9]*.txt')
    txt_files = sorted(glob.glob(pattern))
    
    if not txt_files:
        print(f'WARNING: no baseline text files found in {os.path.abspath(baseline_dir)}')
        return baseline_dict

    # Identify the common reference date from filenames (first part before underscore)
    ref_date_candidates = [os.path.basename(f).split('_')[0] for f in txt_files]
    from collections import Counter
    ref_date = Counter(ref_date_candidates).most_common(1)[0][0]
    
    # Filter files that start with the reference date
    bFiles = [f for f in txt_files if os.path.basename(f).startswith(ref_date)]
    
    # Read each file
    for bFile in bFiles:
        filename = os.path.basename(bFile)
        date_pair = filename.replace('.txt', '')
        dates = date_pair.split('_')
        if len(dates) != 2:
            continue
        date1, date2 = dates
        
        # Parse file content
        bperp = 0.0
        with open(bFile, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('Bperp average'):
                    # Extract numeric value after colon
                    parts = line.split(':')
                    if len(parts) >= 2:
                        try:
                            bperp = float(parts[1].strip())
                        except ValueError:
                            pass
                    break
        
        # For tops, top and bottom are assumed equal (average)
        baseline_dict[date2] = [bperp, bperp]
    
    # Set reference date baseline to [0, 0]
    baseline_dict[ref_date] = [0.0, 0.0]
    
    return baseline_dict

def update_attribute4multilook(atr_in, ref_shape):
    import numpy as np
    atr = atr_in.copy()
    orig_length = int(atr.get('LENGTH', ref_shape[0]))
    orig_width  = int(atr.get('WIDTH', ref_shape[1]))
    yscale = orig_length / ref_shape[0]
    xscale = orig_width / ref_shape[1]

    atr['LENGTH'] = str(ref_shape[0])
    atr['WIDTH'] = str(ref_shape[1])

    if 'ALOOKS' in atr:
        atr['ALOOKS'] = str(int(np.rint(int(float(atr['ALOOKS'])) * yscale)))
    if 'RLOOKS' in atr:
        atr['RLOOKS'] = str(int(np.rint(int(float(atr['RLOOKS'])) * xscale)))

    if 'AZIMUTH_PIXEL_SIZE' in atr:
        atr['AZIMUTH_PIXEL_SIZE'] = str(float(atr['AZIMUTH_PIXEL_SIZE']) * yscale)
    if 'RANGE_PIXEL_SIZE' in atr:
        atr['RANGE_PIXEL_SIZE'] = str(float(atr['RANGE_PIXEL_SIZE']) * xscale)

    if 'NCORRLOOKS' in atr:
        atr['NCORRLOOKS'] = str(float(atr['NCORRLOOKS']) * yscale * xscale)

    return atr

def multilook_geometry_files(geometry_dict, lks_y, lks_x, output_dir=None, overwrite=True):
    """Multilook merged geometry GeoTIFF files to a coarser resolution.

    Parameters
    ----------
    geometry_dict : dict
        Dictionary mapping geometry type to full-resolution file path.
    lks_y : int
        Number of looks in Y / row direction.
    lks_x : int
        Number of looks in X / column direction.
    output_dir : str or Path, optional
        Output directory. If None, files are overwritten in place.
    overwrite : bool
        If True and output_dir is the same as source, overwrite original files.

    Returns
    -------
    dict
        Updated dictionary with paths to multilooked files.
    """
    from osgeo import gdal
    from pathlib import Path

    if output_dir is None:
        output_dir = None   # overwrite in place
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    multilooked = {}
    for gtype, src_file in geometry_dict.items():
        if not src_file or not os.path.isfile(src_file):
            continue

        src_file = Path(src_file)
        # Determine output path
        if output_dir is None:
            out_file = src_file
        else:
            out_file = output_dir / src_file.name

        print(f'Multilooking {gtype} by {lks_y} x {lks_x} -> {out_file}')

        # If overwriting, we need to write to a temporary file first and then replace
        if out_file == src_file:
            temp_file = out_file.with_suffix('.tmp.tif')
        else:
            temp_file = out_file

        src_ds = gdal.Open(str(src_file))
        if src_ds is None:
            print(f'  ERROR: Cannot open {src_file}')
            continue

        gt = src_ds.GetGeoTransform()
        new_gt = (
            gt[0],
            gt[1] * lks_x,
            gt[2],
            gt[3],
            gt[4],
            gt[5] * lks_y
        )

        out_cols = int(src_ds.RasterXSize / lks_x)
        out_rows = int(src_ds.RasterYSize / lks_y)

        driver = gdal.GetDriverByName('GTiff')
        dst_ds = driver.Create(
            str(temp_file),
            out_cols, out_rows,
            src_ds.RasterCount,
            src_ds.GetRasterBand(1).DataType,
            options=['COMPRESS=LZW']
        )
        dst_ds.SetGeoTransform(new_gt)
        dst_ds.SetProjection(src_ds.GetProjection())

        # Use average resampling for floating-point data, nearest for others
        dtype = src_ds.GetRasterBand(1).DataType
        if dtype in (gdal.GDT_Float32, gdal.GDT_Float64):
            resample_alg = gdal.GRA_Average
        else:
            resample_alg = gdal.GRA_NearestNeighbour

        gdal.ReprojectImage(src_ds, dst_ds, None, None, resample_alg)

        # Copy nodata
        for i in range(1, src_ds.RasterCount + 1):
            src_band = src_ds.GetRasterBand(i)
            dst_band = dst_ds.GetRasterBand(i)
            nodata = src_band.GetNoDataValue()
            if nodata is not None:
                dst_band.SetNoDataValue(nodata)
            dst_band.FlushCache()

        src_ds = None
        dst_ds = None

        # Replace original if overwriting
        if out_file == src_file:
            temp_file.replace(src_file)
            out_file = src_file

        multilooked[gtype] = out_file

    return multilooked

def update_metadata_for_multilook(metadata, lks_y, lks_x):
    if metadata is None:
        return None
    meta = metadata.copy()

    # Update look factors (ensure integer strings)
    alooks = int(float(meta.get('ALOOKS', 1))) * lks_y
    rlooks = int(float(meta.get('RLOOKS', 1))) * lks_x
    meta['ALOOKS'] = str(int(alooks))
    meta['RLOOKS'] = str(int(rlooks))

    # Update pixel sizes (keep as float strings)
    if 'AZIMUTH_PIXEL_SIZE' in meta:
        meta['AZIMUTH_PIXEL_SIZE'] = str(float(meta['AZIMUTH_PIXEL_SIZE']) * lks_y)
    if 'RANGE_PIXEL_SIZE' in meta:
        meta['RANGE_PIXEL_SIZE'] = str(float(meta['RANGE_PIXEL_SIZE']) * lks_x)

    if 'NCORRLOOKS' in meta:
        meta['NCORRLOOKS'] = str(float(meta['NCORRLOOKS']) * lks_y * lks_x)

    return meta

def extract_h5_geometry(
    h5_file: Union[str, Path],
    output_dir: Path,
    geom_types: List[str],
    dataset_mapping: Optional[Dict[str, str]] = None,
    x_coords_name: str = "x_coordinates",
    y_coords_name: str = "y_coordinates",
    projection_name: str = "projection"
) -> Dict[str, Dict]:
    """Extract geometry datasets from a static_layers HDF5 file to GeoTIFF.

    Parameters
    ----------
    h5_file : Path or str
        Path to the HDF5 file.
    output_dir : Path
        Directory to save extracted GeoTIFFs.
    geom_types : list of str
        Desired output filenames (e.g., ['height.tif', 'los_east.tif']).
    dataset_mapping : dict, optional
        Mapping from output filename to HDF5 dataset path.
    x_coords_name, y_coords_name, projection_name : str
        Names of coordinate and projection datasets.

    Returns
    -------
    dict
        Keys are geometry types, values are dicts with 'file_list' and 'nodata'.
    """
    if dataset_mapping is None:
        dataset_mapping = {
            "height.tif": "z",
            "layover_shadow_mask.tif": "layover_shadow_mask",
            "local_incidence_angle.tif": "local_incidence_angle",
            "los_east.tif": "los_east",
            "los_north.tif": "los_north",
        }

    extracted = defaultdict(lambda: {'file_list': None, 'nodata': None})
    h5_file = Path(h5_file)
    if not h5_file.exists():
        return extracted

    try:
        with h5py.File(h5_file, 'r') as h5f:
            data_group = h5f['/data']

            # Get EPSG from projection dataset
            epsg = 4326
            if projection_name in data_group:
                proj_ds = data_group[projection_name]
                if 'epsg_code' in proj_ds.attrs:
                    epsg = int(proj_ds.attrs['epsg_code'])
                elif 'spatial_ref' in proj_ds.attrs:
                    wkt = proj_ds.attrs['spatial_ref'].decode('utf-8')
                    srs = osr.SpatialReference()
                    srs.ImportFromWkt(wkt)
                    if srs.GetAuthorityCode(None):
                        epsg = int(srs.GetAuthorityCode(None))

            # Read coordinates
            x_coords = data_group[x_coords_name][:] if x_coords_name in data_group else None
            y_coords = data_group[y_coords_name][:] if y_coords_name in data_group else None
            geotransform = None
            if x_coords is not None and y_coords is not None and len(x_coords) > 1 and len(y_coords) > 1:
                dx = abs(x_coords[1] - x_coords[0])
                dy = abs(y_coords[1] - y_coords[0])
                left = x_coords[0] - dx / 2
                top = y_coords[0] + dy / 2
                geotransform = (left, dx, 0.0, top, 0.0, -dy)

            # Extract each geometry type
            for geom_type in geom_types:
                ds_path = dataset_mapping.get(geom_type)
                if ds_path is None:
                    continue
                if ds_path not in data_group:
                    continue
                dataset = data_group[ds_path]
                data = dataset[:]

                # Determine nodata value
                nodata = None
                if '_FillValue' in dataset.attrs:
                    nodata = dataset.attrs['_FillValue'].item()
                elif np.issubdtype(data.dtype, np.floating):
                    nodata = np.nan
                elif np.issubdtype(data.dtype, np.integer):
                    nodata = np.iinfo(data.dtype).max

                # Write GeoTIFF
                h5_stem = h5_file.stem
                file_output_dir = output_dir / h5_stem
                file_output_dir.mkdir(parents=True, exist_ok=True)
                out_file = file_output_dir / geom_type

                driver = gdal.GetDriverByName('GTiff')
                ds_out = driver.Create(
                    str(out_file),
                    data.shape[1], data.shape[0],
                    1,
                    gdal.GDT_Float32 if data.dtype == np.float32 else gdal.GDT_Float64,
                    options=['COMPRESS=LZW']
                )
                if geotransform:
                    ds_out.SetGeoTransform(geotransform)
                srs = osr.SpatialReference()
                srs.ImportFromEPSG(epsg)
                ds_out.SetProjection(srs.ExportToWkt())
                band = ds_out.GetRasterBand(1)
                band.WriteArray(data)
                if nodata is not None:
                    band.SetNoDataValue(nodata)
                band.FlushCache()
                ds_out = None

                extracted[geom_type]['file_list'] = out_file
                extracted[geom_type]['nodata'] = nodata

    except Exception as e:
        print(f"Error processing {h5_file}: {e}")

    return extracted

def merge_geometry_files(
    burst_ids: List[str],
    geom_dir: Path,
    output_dir: Path,
    geom_types: List[str],
    ref_int_file: Optional[str] = None,
    keep_temp: bool = False
) -> Dict[str, Path]:
    """Merge geometry files across multiple burst IDs and compute incidence angle.

    Parameters
    ----------
    burst_ids : list of str
        List of burst IDs (subdirectory names).
    geom_dir : Path
        Base directory containing burst subdirectories.
    output_dir : Path
        Output directory for merged files.
    geom_types : list of str
        List of geometry file basenames to merge.
    ref_int_file : str, optional
        Path to a reference interferogram GeoTIFF to set output grid.
    keep_temp : bool
        If True, keep temporary extracted files.

    Returns
    -------
    dict
        Mapping from geometry type to merged file path.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Temporary directory for extracted per-burst GeoTIFFs
    if keep_temp:
        temp_dir = output_dir / "temp_extracted"
        temp_dir.mkdir(exist_ok=True)
    else:
        temp_dir = Path(tempfile.mkdtemp())

    # Collect files per geometry type
    geometry_files = defaultdict(list)
    nodata_dict = {}

    for burst_id in burst_ids:
        burst_path = geom_dir / burst_id
        if not burst_path.exists():
            continue
        h5_files = list(burst_path.glob("static_layers*.h5"))
        if not h5_files:
            continue

        for h5_file in h5_files:
            extracted = extract_h5_geometry(h5_file, temp_dir, geom_types)
            for gtype, info in extracted.items():
                if info['file_list']:
                    geometry_files[gtype].append(info['file_list'])
                    if gtype not in nodata_dict and info['nodata'] is not None:
                        nodata_dict[gtype] = info['nodata']

    # Determine output bounds and resolution from reference interferogram if provided
    warp_options = []
    if ref_int_file and os.path.exists(ref_int_file):
        ds = gdal.Open(ref_int_file)
        gt = ds.GetGeoTransform()
        xmin = gt[0]
        ymax = gt[3]
        xmax = xmin + gt[1] * ds.RasterXSize
        ymin = ymax + gt[5] * ds.RasterYSize
        warp_options = [
            '-te', str(xmin), str(ymin), str(xmax), str(ymax),
            '-tr', str(abs(gt[1])), str(abs(gt[5]))
        ]
        ds = None

    # Merge each geometry type using gdalwarp
    merged = {}
    for gtype, file_list in geometry_files.items():
        if not file_list:
            continue
        out_file = output_dir / gtype

        # Remove existing file to avoid gdalwarp error (or use -overwrite)
        if out_file.exists():
            out_file.unlink()

        # Build gdalwarp command
        cmd = [
            'gdalwarp',
            '-overwrite',           # Allow overwriting intermediate temp files if any
            '-of', 'GTiff',
            '-co', 'COMPRESS=LZW'
        ]
        if warp_options:
            cmd.extend(warp_options)
        cmd.extend([str(f) for f in file_list])
        cmd.append(str(out_file))

        print(f"Merging {gtype} from {len(file_list)} bursts using gdalwarp...")
        try:
            # Run gdalwarp; capture output for debugging but print on error
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"Error merging {gtype}:")
            print(f"STDOUT: {e.stdout}")
            print(f"STDERR: {e.stderr}")
            raise
        else:
            if out_file.exists():
                merged[gtype] = out_file

    # Compute incidence angle if los_east and los_north are present
    los_east = merged.get('los_east.tif')
    los_north = merged.get('los_north.tif')
    if los_east and los_north:
        inc_file = output_dir / 'incidenceAngle.tif'
        compute_incidence_angle(
            los_east, los_north, inc_file,
            nodata=nodata_dict.get('los_east.tif')
        )
        merged['incidenceAngle.tif'] = inc_file

    # Clean up temporary directory if not keeping
    if not keep_temp:
        shutil.rmtree(temp_dir, ignore_errors=True)

    return merged


def compute_incidence_angle(los_east_file: Path, los_north_file: Path, output_file: Path, nodata: float = None):
    """Compute incidence angle from LOS east and north components.

    Parameters
    ----------
    los_east_file, los_north_file : Path
        Paths to GeoTIFF files of LOS vector components.
    output_file : Path
        Output path for incidenceAngle.tif.
    nodata : float, optional
        No-data value to use.
    """
    ds_east = gdal.Open(str(los_east_file))
    ds_north = gdal.Open(str(los_north_file))
    east = ds_east.GetRasterBand(1).ReadAsArray()
    north = ds_north.GetRasterBand(1).ReadAsArray()

    # incidence = arccos(up), where up = sqrt(1 - east^2 - north^2)
    up_sq = 1.0 - east**2 - north**2
    up_sq = np.clip(up_sq, 0, 1)
    up = np.sqrt(up_sq)
    inc_angle = np.rad2deg(np.arccos(up))

    # Mask nodata
    if nodata is not None:
        mask = np.isnan(east) | np.isnan(north)
        inc_angle[mask] = nodata

    driver = gdal.GetDriverByName('GTiff')
    ds_out = driver.Create(str(output_file), ds_east.RasterXSize, ds_east.RasterYSize,
                           1, gdal.GDT_Float32, options=['COMPRESS=LZW'])
    ds_out.SetGeoTransform(ds_east.GetGeoTransform())
    ds_out.SetProjection(ds_east.GetProjection())
    band = ds_out.GetRasterBand(1)
    band.WriteArray(inc_angle)
    if nodata is not None:
        band.SetNoDataValue(nodata)
    band.FlushCache()
    ds_out = None
    ds_east = None
    ds_north = None


def extract_merge_geometry(
    geom_dir: str,
    output_dir: str,
    geom_types: List[str],
    ref_int_file: Optional[str] = None,
    metadata: Optional[Dict] = None
) -> Dict[str, Path]:
    """High-level function to extract, merge, and prepare geometry.

    Parameters
    ----------
    geom_dir : str
        Base geometry directory containing burst subdirs.
    output_dir : str
        Output directory.
    geom_types : list
        Desired geometry file names.
    ref_int_file : str, optional
        Reference interferogram GeoTIFF for extent.
    metadata : dict, optional
        Metadata dictionary to be updated with LENGTH and WIDTH.

    Returns
    -------
    dict
        Merged geometry file paths.
    """
    geom_path = Path(geom_dir)
    # Find burst IDs (subdirectories containing static_layers*.h5)
    burst_ids = []
    for subdir in geom_path.iterdir():
        if subdir.is_dir() and list(subdir.glob("static_layers*.h5")):
            burst_ids.append(subdir.name)
    if not burst_ids:
        if list(geom_path.glob("static_layers*.h5")):
            burst_ids = ['.']
        else:
            raise FileNotFoundError(f"No static_layers HDF5 found in {geom_dir}")

    output_path = Path(output_dir)
    merged = merge_geometry_files(
        burst_ids=burst_ids,
        geom_dir=geom_path,
        output_dir=output_path,
        geom_types=geom_types,
        ref_int_file=ref_int_file,
        keep_temp=False
    )

    # Update metadata with dimensions from the merged height file
    if metadata is not None:
        height_file = merged.get('height.tif')
        if height_file and height_file.exists():
            from mintpy.utils import readfile
            geom_atr = readfile.read_attribute(str(height_file))
            metadata['LENGTH'] = geom_atr['LENGTH']
            metadata['WIDTH'] = geom_atr['WIDTH']
            print(f"Updated metadata with LENGTH={metadata['LENGTH']}, WIDTH={metadata['WIDTH']}")

    return merged