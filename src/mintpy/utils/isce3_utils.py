import math
import os
import re
import shutil
import tempfile
import xml.dom.minidom
import xml.etree.ElementTree as ET
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import h5py
import numpy as np
from osgeo import gdal, osr

from mintpy.objects import sensor
from mintpy.utils import readfile

try:
    from scipy.interpolate import CubicHermiteSpline
except ImportError:
    CubicHermiteSpline = None
    print("Warning: scipy not available. Hermite interpolation will fall back to linear.")


# Ensure PROJ/GDAL data paths are set when running without conda activation,
# otherwise osr.ImportFromEPSG() fails with "proj_create_from_database" errors.
def _setup_gdal_proj_data():
    import sys
    for env_key, subdir, probe in [('PROJ_DATA', 'share/proj', 'proj.db'),
                                   ('GDAL_DATA', 'share/gdal', None)]:
        if env_key in os.environ or (env_key == 'PROJ_DATA' and 'PROJ_LIB' in os.environ):
            continue
        data_dir = os.path.join(sys.prefix, *subdir.split('/'))
        if os.path.isdir(data_dir) and (probe is None or os.path.isfile(os.path.join(data_dir, probe))):
            os.environ[env_key] = data_dir

_setup_gdal_proj_data()


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
    # Parse XML file with entity resolution disabled (XXE protection)
    parser = ET.XMLParser(resolve_entities=False)
    tree = ET.parse(meta_file, parser=parser)
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
        try:
            dt = datetime.strptime(sensing_mid, '%Y-%m-%d %H:%M:%S.%f')
            seconds_of_day = dt.hour * 3600.0 + dt.minute * 60.0 + dt.second + dt.microsecond / 1e6
            meta['CENTER_LINE_UTC'] = str(seconds_of_day)
        except (ValueError, TypeError):
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

    # Height / Earth radius (map from lowercase keys used in the XML)
    if meta.get('altitude') and not meta.get('HEIGHT'):
        meta['HEIGHT'] = meta['altitude']
    if meta.get('earthRadius') and not meta.get('EARTH_RADIUS'):
        meta['EARTH_RADIUS'] = meta['earthRadius']

    # Compute center incidence angle if possible
    if meta.get('HEIGHT') and meta.get('EARTH_RADIUS') and meta.get('startingRange'):
        try:
            H = float(meta['HEIGHT'])
            R = float(meta['EARTH_RADIUS'])
            look_angle = np.arcsin(R / (R + H))
            inc_angle = np.arcsin((R + H) / R * np.sin(look_angle))
            meta['CENTER_INCIDENCE_ANGLE'] = str(np.rad2deg(inc_angle))
        except (ValueError, TypeError):
            pass

    # Looks will be computed later from burst vs interferogram dimensions
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
        _, date2 = dates

        # Parse file content (robust: handles Bperp average (m), Bperp (m), bperp, etc.)
        bperp = 0.0
        with open(bFile) as f:
            for line in f:
                line = line.strip()
                match = re.match(r'[Bb]perp\b.*?:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', line)
                if match:
                    try:
                        bperp = float(match.group(1))
                    except ValueError:
                        pass
                    break

        # For tops, top and bottom are assumed equal (average)
        baseline_dict[date2] = [bperp, bperp]

    # Set reference date baseline to [0, 0]
    baseline_dict[ref_date] = [0.0, 0.0]

    return baseline_dict


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
                elif proj_ds.shape == () and np.issubdtype(proj_ds.dtype, np.integer):
                    epsg = int(proj_ds[()])

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
                try:
                    srs = osr.SpatialReference()
                    if srs.ImportFromEPSG(epsg) != 0:
                        raise RuntimeError(f'ImportFromEPSG({epsg}) failed')
                    ds_out.SetProjection(srs.ExportToWkt())
                except Exception as crs_err:
                    print(f'WARNING: could not set CRS EPSG:{epsg} for {geom_type} '
                          f'({crs_err}). Writing GeoTIFF without projection. '
                          f'Check PROJ_DATA environment variable / conda activation.')
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

_GDAL_DTYPE_MAP = {
    'float32': 'Float32',
    'float64': 'Float64',
}


def build_vrt_from_h5(
    h5_file: Union[str, Path],
    output_dir: Path,
    geom_types: List[str],
    dataset_mapping: Optional[Dict[str, str]] = None,
) -> Dict[str, Dict]:
    """Build tiny VRT files referencing HDF5 subdatasets directly (no raster copy).

    Fast alternative to extract_h5_geometry(): only the 1D coordinate arrays and
    dataset attributes are read; the actual raster data stays in the HDF5 file and
    is streamed block-wise by gdalwarp later.

    Parameters
    ----------
    h5_file : Path or str
        Path to the static_layers HDF5 file.
    output_dir : Path
        Directory to save the VRT files (one subdir per HDF5 stem).
    geom_types : list of str
        Desired output filenames (e.g., ['height.tif', 'los_east.tif']).
    dataset_mapping : dict, optional
        Mapping from output filename to HDF5 dataset name under /data.

    Returns
    -------
    dict
        Keys are geometry types, values are dicts with 'file_list' (VRT path)
        and 'nodata'.

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
    h5_file = Path(h5_file).resolve()
    if not h5_file.exists():
        return extracted

    with h5py.File(h5_file, 'r') as h5f:
        data_group = h5f['/data']

        # EPSG code
        epsg = 4326
        if 'projection' in data_group:
            proj_ds = data_group['projection']
            if 'epsg_code' in proj_ds.attrs:
                epsg = int(proj_ds.attrs['epsg_code'])
            elif proj_ds.shape == () and np.issubdtype(proj_ds.dtype, np.integer):
                epsg = int(proj_ds[()])

        # Geotransform from 1D coordinate arrays (cell centers -> top-left corner)
        x_coords = data_group['x_coordinates'][:]
        y_coords = data_group['y_coordinates'][:]
        dx = abs(x_coords[1] - x_coords[0])
        dy = abs(y_coords[1] - y_coords[0])
        left = x_coords[0] - dx / 2
        top = y_coords[0] + dy / 2
        geotransform = (left, dx, 0.0, top, 0.0, -dy)

        srs = osr.SpatialReference()
        if srs.ImportFromEPSG(epsg) != 0:
            raise RuntimeError(f'ImportFromEPSG({epsg}) failed (check PROJ_DATA)')
        srs_wkt = srs.ExportToWkt()

        vrt_dir = Path(output_dir) / h5_file.stem
        vrt_dir.mkdir(parents=True, exist_ok=True)

        for geom_type in geom_types:
            ds_name = dataset_mapping.get(geom_type)
            if ds_name is None or ds_name not in data_group:
                continue
            dataset = data_group[ds_name]
            length, width = dataset.shape

            # Match legacy extract_h5_geometry dtype behavior: float32 or float64
            gdal_dtype = _GDAL_DTYPE_MAP.get(dataset.dtype.name, 'Float32')

            # Determine nodata (same rules as extract_h5_geometry)
            if '_FillValue' in dataset.attrs:
                nodata = dataset.attrs['_FillValue'].item()
            elif np.issubdtype(dataset.dtype, np.floating):
                nodata = np.nan
            elif np.issubdtype(dataset.dtype, np.integer):
                nodata = np.iinfo(dataset.dtype).max
            else:
                nodata = None

            src_fname = f'HDF5:"{h5_file}"://data/{ds_name}'
            nodata_xml = ''
            if nodata is not None:
                nodata_xml = f'    <NoDataValue>{nodata}</NoDataValue>\n'
            vrt_xml = (
                f'<VRTDataset rasterXSize="{width}" rasterYSize="{length}">\n'
                f'  <SRS>{srs_wkt}</SRS>\n'
                f'  <GeoTransform>{", ".join(repr(float(v)) for v in geotransform)}</GeoTransform>\n'
                f'  <VRTRasterBand dataType="{gdal_dtype}" band="1">\n'
                f'{nodata_xml}'
                f'    <SimpleSource>\n'
                f'      <SourceFilename relativeToVRT="0">{src_fname}</SourceFilename>\n'
                f'      <SourceBand>1</SourceBand>\n'
                f'      <SrcRect xOff="0" yOff="0" xSize="{width}" ySize="{length}"/>\n'
                f'      <DstRect xOff="0" yOff="0" xSize="{width}" ySize="{length}"/>\n'
                f'    </SimpleSource>\n'
                f'  </VRTRasterBand>\n'
                f'</VRTDataset>\n'
            )
            vrt_file = vrt_dir / (geom_type + '.vrt')
            with open(vrt_file, 'w') as f:
                f.write(vrt_xml)

            # Sanity check: VRT (and the HDF5 subdataset behind it) must be readable
            ds_check = gdal.Open(str(vrt_file))
            if ds_check is None:
                raise RuntimeError(f'GDAL cannot open generated VRT: {vrt_file}')
            if ds_check.RasterXSize != width or ds_check.RasterYSize != length:
                raise RuntimeError(f'VRT size mismatch for {vrt_file}')
            ds_check = None

            extracted[geom_type]['file_list'] = vrt_file
            extracted[geom_type]['nodata'] = nodata

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

    # Temporary directory for per-burst VRT (fast path) or GeoTIFF (legacy path)
    # files. Created next to the output dir to avoid exhausting the system /tmp.
    if keep_temp:
        temp_dir = output_dir / "temp_extracted"
        temp_dir.mkdir(exist_ok=True)
    else:
        temp_dir = Path(tempfile.mkdtemp(prefix='geom_merge_', dir=str(output_dir)))

    def _collect(use_vrt):
        """Collect per-burst sources (VRT or extracted GeoTIFF) per geometry type."""
        geometry_files = defaultdict(list)
        nodata_dict = {}
        for burst_id in burst_ids:
            burst_path = geom_dir / burst_id
            if not burst_path.exists():
                continue
            for h5_file in sorted(burst_path.glob("static_layers*.h5")):
                if use_vrt:
                    extracted = build_vrt_from_h5(h5_file, temp_dir, geom_types)
                else:
                    extracted = extract_h5_geometry(h5_file, temp_dir, geom_types)
                for gtype, info in extracted.items():
                    if info['file_list']:
                        geometry_files[gtype].append(info['file_list'])
                        if gtype not in nodata_dict and info['nodata'] is not None:
                            nodata_dict[gtype] = info['nodata']
        return geometry_files, nodata_dict

    # Fast path: VRTs referencing HDF5 subdatasets directly, no full-res raster
    # copy. Fall back to the legacy GeoTIFF extraction on any failure.
    try:
        geometry_files, nodata_dict = _collect(use_vrt=True)
    except Exception as e:
        print(f'WARNING: fast VRT-based merge failed: {e}')
        print('         Falling back to legacy per-burst GeoTIFF extraction ...')
        geometry_files, nodata_dict = _collect(use_vrt=False)

    if not any(geometry_files.values()):
        print('#' * 60)
        print('WARNING: no geometry layers extracted from any static_layers*.h5!')
        print('         Merged geometry files will NOT be (re)generated;')
        print('         stale files in the output directory may be reused downstream.')
        print('         Check the "Error processing ..." messages above.')
        print('#' * 60)

    # Determine output bounds and resolution from reference interferogram if provided
    bounds, xres, yres, expected_size = None, None, None, None
    if ref_int_file and os.path.exists(ref_int_file):
        ds = gdal.Open(str(ref_int_file))
        gt = ds.GetGeoTransform()
        xmin = gt[0]
        ymax = gt[3]
        xmax = xmin + gt[1] * ds.RasterXSize
        ymin = ymax + gt[5] * ds.RasterYSize
        bounds = (xmin, ymin, xmax, ymax)
        xres, yres = abs(gt[1]), abs(gt[5])
        expected_size = (ds.RasterXSize, ds.RasterYSize)
        ds = None

    # Merge each geometry type using gdal.Warp (python API, streams block-wise)
    merged = {}
    for gtype, file_list in geometry_files.items():
        if not file_list:
            continue
        out_file = output_dir / gtype
        if out_file.exists():
            out_file.unlink()

        print(f"Merging {gtype} from {len(file_list)} bursts using gdal.Warp...")
        warp_kwargs = dict(
            format='GTiff',
            creationOptions=['COMPRESS=LZW'],
            resampleAlg='near',
            multithread=True,
            warpOptions=['NUM_THREADS=ALL_CPUS'],
        )
        if bounds is not None:
            warp_kwargs.update(outputBounds=bounds, xRes=xres, yRes=yres)
        nodata = nodata_dict.get(gtype)
        if nodata is not None:
            warp_kwargs.update(dstNodata=nodata)

        ds_out = gdal.Warp(str(out_file), [str(f) for f in file_list], **warp_kwargs)
        if ds_out is None:
            raise RuntimeError(f'gdal.Warp failed for {gtype} '
                               f'(inputs: {[str(f) for f in file_list]})')

        # Validate output grid and report coverage
        out_size = (ds_out.RasterXSize, ds_out.RasterYSize)
        if expected_size is not None and out_size != expected_size:
            raise RuntimeError(f'merged {gtype} size {out_size} does not match '
                               f'reference interferogram size {expected_size}')
        arr = ds_out.GetRasterBand(1).ReadAsArray()
        if nodata is not None and not (isinstance(nodata, float) and np.isnan(nodata)):
            valid_ratio = float(np.mean(arr != nodata))
        else:
            valid_ratio = float(np.mean(np.isfinite(arr)))
        print(f'    size: {out_size[0]} x {out_size[1]}, valid pixels: {valid_ratio:.1%}')
        ds_out = None
        merged[gtype] = out_file

    # Compute incidence angle and azimuth angle if los_east and los_north are present
    los_east = merged.get('los_east.tif')
    los_north = merged.get('los_north.tif')
    if los_east and los_north:
        inc_file = output_dir / 'incidenceAngle.tif'
        compute_incidence_angle(
            los_east, los_north, inc_file,
            nodata=nodata_dict.get('los_east.tif')
        )
        merged['incidenceAngle.tif'] = inc_file

        az_file = output_dir / 'azimuthAngle.tif'
        compute_azimuth_angle(
            los_east, los_north, az_file,
            nodata=nodata_dict.get('los_east.tif')
        )
        merged['azimuthAngle.tif'] = az_file

    # Clean up temporary directory if not keeping
    if not keep_temp:
        shutil.rmtree(temp_dir, ignore_errors=True)

    return merged

def compute_azimuth_angle(los_east_file: Path, los_north_file: Path, output_file: Path, nodata: float = None):
    """Compute azimuth angle from LOS east and north components.

    The azimuth angle is defined as the angle from the North, measured
    anti‑clockwise as positive (standard mathematical convention).

    Parameters
    ----------
    los_east_file, los_north_file : Path
        Paths to GeoTIFF files of LOS vector components.
    output_file : Path
        Output path for azimuthAngle.tif.
    nodata : float, optional
        No-data value to use.

    """
    ds_east = gdal.Open(str(los_east_file))
    ds_north = gdal.Open(str(los_north_file))
    east = ds_east.GetRasterBand(1).ReadAsArray()
    north = ds_north.GetRasterBand(1).ReadAsArray()

    # Azimuth from east/north components:
    # az = -arctan2(east, north) * 180/pi  (mod 360)
    az_angle = -1 * np.rad2deg(np.arctan2(east, north)) % 360.0

    # Mask nodata
    if nodata is not None:
        mask = np.isnan(east) | np.isnan(north)
        az_angle[mask] = nodata

    driver = gdal.GetDriverByName('GTiff')
    ds_out = driver.Create(str(output_file), ds_east.RasterXSize, ds_east.RasterYSize,
                           1, gdal.GDT_Float32, options=['COMPRESS=LZW'])
    ds_out.SetGeoTransform(ds_east.GetGeoTransform())
    ds_out.SetProjection(ds_east.GetProjection())
    band = ds_out.GetRasterBand(1)
    band.WriteArray(az_angle)
    if nodata is not None:
        band.SetNoDataValue(nodata)
    band.FlushCache()
    ds_out = None
    ds_east = None
    ds_north = None

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
    metadata: Optional[Dict] = None,
    extra_dirs: Optional[List[str]] = None
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
    extra_dirs : list of str, optional
        Additional directories containing static_layers*.h5 files
        (e.g. individual burst dirs from glob expansion).

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

    # When extra dirs are provided (multi-burst from glob expansion),
    # build a temporary directory tree with symlinks so that
    # merge_geometry_files can process all bursts together.
    base_dir = geom_path
    temp_base_dir = None
    if extra_dirs:
        temp_base_dir = Path(tempfile.mkdtemp())
        new_burst_ids = []
        for bid in burst_ids:
            src = (geom_path / bid).resolve()
            dst_name = f"_m_{bid}".replace('.', 'r') if bid == '.' else f"_m_{bid}"
            os.symlink(src, temp_base_dir / dst_name, target_is_directory=True)
            new_burst_ids.append(dst_name)
        for i, edir in enumerate(extra_dirs):
            edir_path = Path(edir)
            if edir_path.is_dir() and list(edir_path.glob("static_layers*.h5")):
                dst_name = f"_x_{i}"
                os.symlink(edir_path.resolve(), temp_base_dir / dst_name, target_is_directory=True)
                new_burst_ids.append(dst_name)
        burst_ids = new_burst_ids
        base_dir = temp_base_dir

    if not burst_ids:
        raise FileNotFoundError(f"No static_layers HDF5 found in {geom_dir}")

    output_path = Path(output_dir)
    merged = merge_geometry_files(
        burst_ids=burst_ids,
        geom_dir=base_dir,
        output_dir=output_path,
        geom_types=geom_types,
        ref_int_file=ref_int_file,
        keep_temp=False
    )

    if temp_base_dir:
        shutil.rmtree(temp_base_dir, ignore_errors=True)

    # Update metadata with dimensions from the merged height file
    if metadata is not None:
        height_file = merged.get('height.tif')
        if height_file and height_file.exists():
            geom_atr = readfile.read_attribute(str(height_file))
            metadata['LENGTH'] = geom_atr['LENGTH']
            metadata['WIDTH'] = geom_atr['WIDTH']
            print(f"Updated metadata with LENGTH={metadata['LENGTH']}, WIDTH={metadata['WIDTH']}")

    return merged

###############################################################################
# XML generation for ISCE3 burst metadata (optional)
###############################################################################
def _to_seconds(t_str: str, ref_epoch_str: str) -> float:
    """Convert datetime string to seconds relative to reference epoch."""
    t = datetime.strptime(t_str, '%Y-%m-%d %H:%M:%S.%f')
    ref = datetime.strptime(ref_epoch_str, '%Y-%m-%d %H:%M:%S.%f')
    return (t - ref).total_seconds()


def _compute_heading(state_vector):

    """
    Calculate ENU heading angle from satellite state vectors.

    Parameters
    ----------
    state_vector : tuple
        (position, velocity) in ECEF coordinates.
        position: [x, y, z] in meters.
        velocity: [vx, vy, vz] in m/s.

    Returns
    -------
    heading : float
        ENU heading angle in degrees (0-360), clockwise from North.

    """
    # Unpack state vector
    state_pos, state_vel = state_vector

    # Convert ECEF position to LLH
    # Using WGS84 ellipsoid parameters
    a = 6378137.0  # semi-major axis
    f = 1.0 / 298.257223563  # flattening
    b = a * (1 - f)  # semi-minor axis
    e2 = 1 - (b**2 / a**2)  # eccentricity squared

    x, y, z = state_pos

    # Longitude
    lon = np.arctan2(y, x)

    # Latitude using iterative method
    p = np.sqrt(x**2 + y**2)
    lat = np.arctan2(z, p * (1 - e2))

    # Iterate to improve latitude accuracy
    for _ in range(10):
        N = a / np.sqrt(1 - e2 * np.sin(lat)**2)
        h = p / np.cos(lat) - N
        lat_new = np.arctan2(z, p * (1 - e2 * N / (N + h)))
        if np.abs(lat_new - lat) < 1e-12:
            break
        lat = lat_new

    # Altitude
    N = a / np.sqrt(1 - e2 * np.sin(lat)**2)

    # Calculate ENU basis vectors at satellite position
    sin_lat = np.sin(lat)
    cos_lat = np.cos(lat)
    sin_lon = np.sin(lon)
    cos_lon = np.cos(lon)

    # ECEF to ENU rotation matrix
    R = np.array([
        [-sin_lon,           cos_lon,           0.0],
        [-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat],
        [cos_lat * cos_lon,  cos_lat * sin_lon, sin_lat]
    ])

    # Rotate velocity from ECEF to ENU
    vel_enu = np.dot(R, state_vel)

    # Extract East and North components
    v_east = vel_enu[0]
    v_north = vel_enu[1]

    # Calculate heading angle
    heading_rad = np.arctan2(v_east, v_north)
    heading_deg = np.degrees(heading_rad)

    # Normalize to 0-360
    if heading_deg < 0:
        heading_deg += 360.0

    return heading_deg


def _orbit_interp_hermite(metadata, time):

    """
    Interpolate orbit state vectors at given time using Hermite interpolation.

    Parameters
    ----------
    metadata : dict
        Dictionary containing orbit metadata with keys:
        - 'time': array of times in seconds relative to reference_epoch
        - 'position_x', 'position_y', 'position_z': arrays of position components
        - 'velocity_x', 'velocity_y', 'velocity_z': arrays of velocity components
        - 'reference_epoch': reference epoch datetime string
    time : float or array_like
        Time(s) in seconds relative to reference_epoch for interpolation

    Returns
    -------
    dict
        Dictionary containing interpolated position and velocity at given time(s)

    """
    # Extract time array
    t_array = np.array(metadata['time'])

    # Create Hermite interpolators for each component
    # Position interpolation with velocity as derivatives
    pos_x_interp = CubicHermiteSpline(t_array, metadata['position_x'], metadata['velocity_x'])
    pos_y_interp = CubicHermiteSpline(t_array, metadata['position_y'], metadata['velocity_y'])
    pos_z_interp = CubicHermiteSpline(t_array, metadata['position_z'], metadata['velocity_z'])

    # For velocity, we can either use the derivative of position interpolators
    # or create separate velocity interpolators. Using derivative ensures consistency.

    # Calculate interpolated values
    if np.isscalar(time):
        position = np.array([
            float(pos_x_interp(time)),
            float(pos_y_interp(time)),
            float(pos_z_interp(time))
        ])
        # Get velocity from derivative of position interpolators
        velocity = np.array([
            float(pos_x_interp.derivative()(time)),
            float(pos_y_interp.derivative()(time)),
            float(pos_z_interp.derivative()(time))
        ])
    else:
        time_array = np.asarray(time)
        position = np.column_stack([
            pos_x_interp(time_array),
            pos_y_interp(time_array),
            pos_z_interp(time_array)
        ])
        # Get velocity from derivative of position interpolators
        velocity = np.column_stack([
            pos_x_interp.derivative()(time_array),
            pos_y_interp.derivative()(time_array),
            pos_z_interp.derivative()(time_array)
        ])

    return position, velocity


def read_burst_metadata_h5(
    h5_file: Path,
    layer_names: List[str] = None,
    group_path: str = "/metadata/processing_information/input_burst_metadata/"
) -> Dict[str, Any]:

    """
    Read metadata for burst attributes from an HDF5 file.

    Parameters
    ----------
    h5_file : Path
        Path to the HDF5 file
    layer_names : List[str], optional
        List of burst metadata attribute names to read.
        If None, all datasets directly under the group will be read.
        Defaults to None.
    group_path : str, optional
        Path to the burst metadata group in the HDF5 file.
        Defaults to '/metadata/processing_information/input_burst_metadata/'

    Returns
    -------
    Dict[str, Any]
        Metadata dictionary for the burst attributes

    """
    metadata = {}

    try:
        with h5py.File(h5_file, 'r') as f:
            # Check if group exists
            if group_path not in f:
                print(f"Group {group_path} not found in {h5_file}")
                return metadata

            # Get the group object
            group = f[group_path]

            # If layer_names is not provided, get all direct datasets in the group
            if layer_names is None:
                # Get all items in the group, filter only datasets (not subgroups)
                layer_names = []
                for name, item in group.items():
                    if isinstance(item, h5py.Dataset):
                        layer_names.append(name)

            # Read each requested layer
            for layer_name in layer_names:
                dataset_path = f"{group_path.rstrip('/')}/{layer_name}"

                if dataset_path in f:
                    # Get the dataset
                    dataset = f[dataset_path]

                    # Read the value
                    value = dataset[()]

                    # Handle different data types
                    if isinstance(value, np.ndarray):
                        # For string arrays, decode bytes to string
                        if value.dtype.kind == 'S' or value.dtype.kind == 'O':
                            if value.size == 1:
                                metadata[layer_name] = value.item().decode('utf-8') if isinstance(value.item(), bytes) else value.item()
                            else:
                                metadata[layer_name] = [v.decode('utf-8') if isinstance(v, bytes) else v for v in value.tolist()]
                        else:
                            # For shape array, extract length and width
                            if layer_name == 'shape' and value.shape == (2,):
                                metadata['length'] = int(value[0])
                                metadata['width'] = int(value[1])
                                metadata['shape'] = value.tolist()
                            else:
                                metadata[layer_name] = value.tolist()
                    elif isinstance(value, bytes):
                        metadata[layer_name] = value.decode('utf-8')
                    else:
                        # Convert scalar numpy types to Python types
                        metadata[layer_name] = value.item() if hasattr(value, 'item') else value

                    # Optionally add dataset attributes if they exist
                    if dataset.attrs:
                        metadata[f"{layer_name}_attrs"] = {
                            key: (val.tolist() if isinstance(val, np.ndarray) else val)
                            for key, val in dataset.attrs.items()
                        }
                else:
                    print(f"Dataset {dataset_path} not found in {h5_file}")

    except Exception as e:
        print(f"Failed to read burst metadata from {h5_file}: {e}")
        return metadata

    # Calculate sensing_mid if we have sensing_start and sensing_stop
    if 'sensing_start' in metadata and 'sensing_stop' in metadata:
        try:
            # Parse datetime strings
            start_dt = datetime.strptime(metadata['sensing_start'], '%Y-%m-%d %H:%M:%S.%f')
            stop_dt = datetime.strptime(metadata['sensing_stop'], '%Y-%m-%d %H:%M:%S.%f')

            # Calculate time difference and mid time
            time_diff = stop_dt - start_dt
            mid_dt = start_dt + (time_diff / 2)

            # Format back to string
            metadata['sensing_mid'] = mid_dt.strftime('%Y-%m-%d %H:%M:%S.%f')

        except Exception as e:
            print(f"Failed to calculate sensing_mid: {e}")
            # If calculation fails, try to read it directly from file
            try:
                dataset_path = f"{group_path.rstrip('/')}/sensing_mid"
                with h5py.File(h5_file, 'r') as f:
                    if dataset_path in f:
                        value = f[dataset_path][()]
                        if isinstance(value, bytes):
                            metadata['sensing_mid'] = value.decode('utf-8')
                        elif hasattr(value, 'item'):
                            metadata['sensing_mid'] = value.item()
                        else:
                            metadata['sensing_mid'] = value
            except Exception as read_e:
                print(f"Failed to read sensing_mid from file: {read_e}")

    # Calculate mid_range if we have starting_range, width, and range_pixel_spacing
    if all(key in metadata for key in ['starting_range', 'width', 'range_pixel_spacing']):
        try:
            mid_range = (metadata['starting_range'] +
                        (metadata['width'] / 2) *
                        metadata['range_pixel_spacing'])
            metadata['mid_range'] = mid_range
        except Exception as e:
            print(f"Failed to calculate mid_range: {e}")

    return metadata


def prepare_mintpy_metadata(metafile: Path) -> Dict[str, Any]:

    """
    Prepare metadata dictionary from MintPy HDF5 file for processing.

    This function extracts specific metadata groups from a MintPy HDF5 file,
    combines them into a single dictionary, and extracts additional derived
    metadata fields such as swath number.

    Parameters
    ----------
    metafile : Path
        Path to the MintPy HDF5 metadata file.

    Returns
    -------
    Dict[str, Any]
        Combined metadata dictionary containing:
        - All metadata from specified HDF5 groups
        - Derived fields like 'swathNumber'

    Notes
    -----
    The function reads metadata from three predefined HDF5 group paths:
    1. Processing information and input burst metadata
    2. Orbit information
    3. Identification information

    The swath number is extracted from the 'burst_id' field using regex
    pattern matching for IW1, IW2, or IW3 swaths.

    """
    # Define HDF5 group paths to extract metadata from
    group_path_list = [
        '/metadata/processing_information/input_burst_metadata/',
        '/metadata/orbit',
        '/identification'
    ]

    # Initialize empty metadata dictionary
    metadata = {}

    # Iterate through each group path and read metadata
    for group_path in group_path_list:
        # Update metadata dictionary with contents from current group
        metadata.update(read_burst_metadata_h5(metafile, group_path=group_path))

    # Extract swath number from burst_id using regex pattern matching
    # Pattern matches iw0, iw1, or iw2 (case-insensitive) and extracts the digit
    if re.search(r"iw([012])", metadata.get('burst_id', ''), re.IGNORECASE):
        metadata['swathNumber'] = int(
            re.search(r"iw([012])", metadata['burst_id'], re.IGNORECASE).group(1)
        )
    else:
        metadata['swathNumber'] = None

    return metadata


def extract_required_attributes(metadata):
    """Extract only the burst attributes needed by extract_tops_metadata."""
    isce3_available = True
    try:
        import isce3
    except ImportError:
        isce3_available = False
        print("WARNING: isce3 not available. Some metadata fields will use defaults.")

    meta = {}

    # Direct burst attributes used in original function (with fallbacks)
    meta['prf'] = metadata.get('prf_raw_data',
                    metadata.get('prf', 1717.0))
    meta['burstStartUTC'] = metadata.get('sensing_start',
                             metadata.get('zero_doppler_start_time', ''))
    meta['burstStopUTC'] = metadata.get('sensing_stop',
                            metadata.get('zero_doppler_end_time', ''))
    meta['radarWavelength'] = metadata.get('wavelength',
                               metadata.get('radar_wavelength', 0.05546576))
    meta['startingRange'] = metadata.get('starting_range',
                             metadata.get('slant_range_time', 800000.0))
    # Convert slant_range_time (seconds) to meters if needed
    if isinstance(meta['startingRange'], (int, float)):
        try:
            if meta['startingRange'] < 1.0:
                SPEED_OF_LIGHT = 299792458.0
                meta['startingRange'] = meta['startingRange'] * SPEED_OF_LIGHT / 2.0
        except Exception:
            meta['startingRange'] = meta.get('starting_range', 800000.0)
    meta['passDirection'] = metadata.get('orbit_direction',
                             metadata.get('orbit_pass_direction', 'ascending'))
    meta['polarization'] = metadata.get('polarization', 'VV')
    meta['trackNumber'] = metadata.get('track_number', 0)
    meta['orbitNumber'] = metadata.get('absolute_orbit_number', 0)

    # Additional attributes needed for calculations
    meta['sensingMid'] = metadata.get('sensing_mid', metadata.get('sensing_start', ''))
    meta['azimuthTimeInterval'] = metadata.get('azimuth_time_interval',
                                   metadata.get('azimuth_time_interval_', 0.002))
    meta['rangePixelSize'] = metadata.get('range_pixel_spacing',
                              metadata.get('range_pixel_spacing_', 2.3))
    meta['swathNumber'] = metadata.get('swathNumber',
                           metadata.get('swath_number', 2))
    meta['ascendingNodeTime'] = None

    # Calculate satellite speed (Vs) from orbit data (if available)
    has_orbit = all(k in metadata for k in ['time', 'position_x', 'velocity_x'])
    if has_orbit and isce3_available:
        try:
            sensingMid = _to_seconds(metadata['sensing_mid'], metadata['reference_epoch'])
            sv = _orbit_interp_hermite(metadata, sensingMid)
            velocity = np.linalg.norm(sv[1])
            position = sv[0]
            ellipsoid = isce3.core.Ellipsoid()
            llh = ellipsoid.xyz_to_lon_lat(position)
            heading = _compute_heading(sv)
            meta['satelliteSpeed'] = velocity
            meta['position'] = position
            meta['HEADING'] = heading
            meta['earthRadius'] = ellipsoid.r_dir(math.radians(heading), llh[1])
            meta['altitude'] = llh[2]
        except Exception as e:
            print(f"WARNING: orbit interpolation failed: {e}. Using default values.")
            isce3_available = False

    if not has_orbit or not isce3_available:
        meta['satelliteSpeed'] = 7545.0
        meta['position'] = [0.0, 0.0, 0.0]
        meta['HEADING'] = 0.0
        meta['earthRadius'] = 6371000.0
        meta['altitude'] = 693000.0

    # Calculate frame numbers
    if meta['ascendingNodeTime'] is not None and meta['burstStartUTC']:
        try:
            from datetime import datetime
            start_dt = datetime.strptime(meta['burstStartUTC'], '%Y-%m-%d %H:%M:%S.%f')
            if isinstance(meta['ascendingNodeTime'], str):
                node_dt = datetime.strptime(meta['ascendingNodeTime'], '%Y-%m-%d %H:%M:%S.%f')
            else:
                node_dt = meta['ascendingNodeTime']
            time_diff_start = (start_dt - node_dt).total_seconds()
            time_diff_stop = (datetime.strptime(meta['burstStopUTC'], '%Y-%m-%d %H:%M:%S.%f') - node_dt).total_seconds()
            meta['firstFrameNumber'] = int(0.2 * time_diff_start)
            meta['lastFrameNumber'] = int(0.2 * time_diff_stop)
        except Exception:
            meta['firstFrameNumber'] = 0
            meta['lastFrameNumber'] = 0
    else:
        meta['firstFrameNumber'] = 0
        meta['lastFrameNumber'] = 0

    return meta


def save_burst_attributes_to_xml(burst_attrs: Dict[str, Any], output_path: Union[str, Path]) -> bool:
    """Save burst attributes to XML file in MintPy-compatible format."""
    root = ET.Element('burst_metadata')
    info = ET.SubElement(root, 'info')
    ET.SubElement(info, 'generation_time').text = datetime.now().isoformat() + 'Z'
    ET.SubElement(info, 'source').text = 'static_layers extraction'

    attrs_elem = ET.SubElement(root, 'burst_attributes')
    attributes_to_save = [
        'prf', 'burstStartUTC', 'burstStopUTC', 'radarWavelength',
        'startingRange', 'passDirection', 'polarization', 'trackNumber',
        'orbitNumber', 'sensingMid', 'azimuthTimeInterval', 'rangePixelSize',
        'swathNumber', 'ascendingNodeTime', 'satelliteSpeed', 'position',
        'HEADING', 'earthRadius', 'altitude', 'firstFrameNumber', 'lastFrameNumber'
    ]
    for key in attributes_to_save:
        if key in burst_attrs:
            elem = ET.SubElement(attrs_elem, key)
            value = burst_attrs[key]
            if value is None:
                elem.text = 'None'
                elem.set('type', 'NoneType')
            elif isinstance(value, (int, np.integer)):
                elem.text = str(int(value))
                elem.set('type', 'int')
            elif isinstance(value, (float, np.floating)):
                elem.text = f"{value:.12f}"
                elem.set('type', 'float')
            elif isinstance(value, str):
                elem.text = value
                elem.set('type', 'str')
            elif isinstance(value, np.ndarray):
                elem.text = str(value)
                elem.set('type', 'np.ndarray')
            else:
                elem.text = str(value)
                elem.set('type', 'other')

    xml_str = ET.tostring(root, encoding='utf-8')
    dom = xml.dom.minidom.parseString(xml_str)
    pretty_xml = dom.toprettyxml(indent='  ')
    lines = [line for line in pretty_xml.split('\n') if line.strip()]
    pretty_xml = '\n'.join(lines)
    with open(output_path, 'w') as f:
        f.write(pretty_xml)
    return True


def generate_burst_xml_from_static(static_h5_file: Union[str, Path], output_xml: Union[str, Path]) -> str:
    """Generate burst XML file from a single static_layers HDF5."""
    print(f'Generating burst XML from {static_h5_file}')
    meta = prepare_mintpy_metadata(static_h5_file)
    attrs = extract_required_attributes(meta)
    save_burst_attributes_to_xml(attrs, output_xml)
    print(f'Saved burst XML to {output_xml}')
    return str(output_xml)
