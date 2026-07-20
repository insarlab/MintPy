import glob
import os
import re
from pathlib import Path

from mintpy.utils import isce3_utils, ptime, readfile, writefile


#########################################################################
def add_ifgram_metadata(metadata_in, dates=[], baseline_dict={}):
    """Add metadata unique for each interferogram.

    Parameters: metadata_in   : dict, input common metadata for the entire dataset
                dates         : list of str in YYYYMMDD format
                baseline_dict : dict, output of baseline_timeseries()
    Returns:    metadata      : dict, updated metadata
    """
    metadata = metadata_in.copy()
    metadata['DATE12'] = f'{dates[0][2:]}-{dates[1][2:]}'

    if baseline_dict:
        if dates[0] in baseline_dict and dates[1] in baseline_dict:
            bperp_top = baseline_dict[dates[1]][0] - baseline_dict[dates[0]][0]
            bperp_bottom = baseline_dict[dates[1]][1] - baseline_dict[dates[0]][1]
            metadata['P_BASELINE_TOP_HDR'] = str(bperp_top)
            metadata['P_BASELINE_BOTTOM_HDR'] = str(bperp_bottom)
    return metadata


def prepare_geometry_isce3(geom_dir, out_dir, geom_files=None, metadata=None,
                           processor='tops', update_mode=True, ref_int_file=None,
                           target_shape=None, geom_dirs=None):
    """Prepare geometry files from ISCE3/Dolphin static_layers HDF5.

    Parameters
    ----------
    target_shape : tuple of (length, width), optional
        Interferogram dimensions.
    """
    from pathlib import Path

    print('preparing geometry files from ISCE3/Dolphin static layers')
    geom_dir = os.path.abspath(geom_dir)
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    if geom_files is None:
        geom_files = ['height.tif', 'los_east.tif', 'los_north.tif',
                      'layover_shadow_mask.tif', 'local_incidence_angle.tif']

    # Step 0: Read full-resolution pixel size from first static_layers.h5
    burst_full_dx = None
    burst_full_dy = None
    num_bursts = 0
    geom_path = Path(geom_dir)
    burst_subdirs = sorted([d for d in geom_path.iterdir()
                            if d.is_dir() and list(d.glob('static_layers*.h5'))])
    h5_files_root = sorted(geom_path.glob('static_layers*.h5'))
    first_h5_path = None
    if burst_subdirs:
        first_h5_path = sorted(burst_subdirs[0].glob('static_layers*.h5'))[0]
    elif h5_files_root:
        first_h5_path = h5_files_root[0]
    num_bursts = len(burst_subdirs) if burst_subdirs else (1 if h5_files_root else 0)
    if geom_dirs and len(geom_dirs) > 1:
        for edir in geom_dirs[1:]:
            epath = Path(edir)
            if epath.is_dir():
                if list(epath.glob('static_layers*.h5')):
                    num_bursts += 1
                else:
                    for sub in epath.iterdir():
                        if sub.is_dir() and list(sub.glob('static_layers*.h5')):
                            num_bursts += 1
    if first_h5_path:
        try:
            import h5py
            with h5py.File(first_h5_path, 'r') as h5:
                x_coords = h5['/data/x_coordinates'][:]
                y_coords = h5['/data/y_coordinates'][:]
            burst_full_dx = abs(x_coords[1] - x_coords[0])
            burst_full_dy = abs(y_coords[1] - y_coords[0])
            print(f'Number of bursts: {num_bursts}')
            print(f'Full-resolution pixel size: dx={burst_full_dx}, dy={burst_full_dy}')
        except Exception as e:
            print(f'WARNING: could not read full-res pixel size: {e}')

    # Step 1: Merge and crop geometry to interferogram extent and resolution
    extra_dirs = geom_dirs[1:] if (geom_dirs and len(geom_dirs) > 1) else None
    geometry_dict = isce3_utils.extract_merge_geometry(
        geom_dir=geom_dir,
        output_dir=out_dir,
        geom_types=geom_files,
        ref_int_file=ref_int_file,
        metadata=None,
        extra_dirs=extra_dirs,
    )

    # Step 2: Compute ALOOKS/RLOOKS from full-res pixel size vs target pixel size
    if metadata is not None:
        lks_y = 1
        lks_x = 1

        if burst_full_dx is not None and burst_full_dy is not None and ref_int_file:
            from osgeo import gdal
            ref_ds = gdal.Open(str(ref_int_file))
            if ref_ds:
                ref_gt = ref_ds.GetGeoTransform()
                ref_dx = abs(ref_gt[1])
                ref_dy = abs(ref_gt[5])
                ref_ds = None
                lks_x = max(1, int(round(ref_dx / burst_full_dx)))
                lks_y = max(1, int(round(ref_dy / burst_full_dy)))

        metadata['ALOOKS'] = str(lks_y)
        metadata['RLOOKS'] = str(lks_x)
        print(f'Set ALOOKS={lks_y}, RLOOKS={lks_x} '
              f'(full-res dx={burst_full_dx}, dy={burst_full_dy})')

    # Write .rsc files
    for geom_name, geom_path in geometry_dict.items():
        if geom_path and os.path.isfile(geom_path):
            geom_path = str(geom_path)
            rsc_file = geom_path + '.rsc'
            geom_meta = metadata.copy() if metadata else {}
            geom_meta.update(readfile.read_attribute(geom_path))
            writefile.write_roipac_rsc(geom_meta, rsc_file,
                                       update_mode=update_mode,
                                       print_msg=True)

    # Update metadata with final dimensions
    height_file = geometry_dict.get('height.tif')
    if height_file and os.path.isfile(height_file):
        geom_atr = readfile.read_attribute(str(height_file))
        if metadata is not None:
            metadata['LENGTH'] = geom_atr['LENGTH']
            metadata['WIDTH'] = geom_atr['WIDTH']
            print(f"Updated metadata with LENGTH={metadata['LENGTH']}, WIDTH={metadata['WIDTH']}")

    return metadata


def prepare_stack_isce3(obs_file, metadata=None, baseline_dict=None, update_mode=True):
    print(f'preparing RSC file for: {obs_file}')
    tif_files = sorted(glob.glob(obs_file))
    if not tif_files:
        raise FileNotFoundError(f'No file found with pattern: {obs_file}')

    meta = metadata.copy() if metadata else {}
    num_file = len(tif_files)

    # Normalize ALOOKS/RLOOKS to integer strings if present
    for key in ['ALOOKS', 'RLOOKS']:
        if key in meta:
            try:
                meta[key] = str(int(float(meta[key])))
            except (ValueError, TypeError):
                pass

    prog_bar = ptime.progressBar(maxValue=num_file, print_msg=(num_file > 5))
    num_valid = 0
    for i, tif_file in enumerate(tif_files):
        try:
            fbase = os.path.basename(tif_file)

            # Extract date pair(s) from filename
            date_nums = re.findall(r'\d{8}', fbase)
            if len(date_nums) < 2:
                prog_bar.update(i+1, suffix=f'skipped {i+1}/{num_file}')
                continue
            dates = sorted(set(date_nums))[:2]

            if not all(19000101 <= int(d) <= 20991231 for d in dates):
                prog_bar.update(i+1, suffix=f'skipped {i+1}/{num_file}')
                continue

            num_valid += 1
            prog_bar.update(i+1, suffix=f'{dates[0]}_{dates[1]} {i+1}/{num_file}')

            ifg_meta = meta.copy()
            ifg_meta.update(readfile.read_attribute(tif_file))
            ifg_meta = add_ifgram_metadata(ifg_meta, dates, baseline_dict)

            rsc_file = tif_file + '.rsc'
            writefile.write_roipac_rsc(ifg_meta, rsc_file,
                                       update_mode=update_mode,
                                       print_msg=False)

        except Exception as e:
            prog_bar.update(i+1, suffix=f'error {i+1}/{num_file}')
            continue

    prog_bar.close()
    if num_valid == 0:
        print(f'WARNING: no valid files processed for pattern: {obs_file}')
    return


#########################################################################
def prep_isce3(inps):
    """Prepare ISCE3/Dolphin metadata files."""
    # If no meta file is provided or it is 'auto', generate one from static_layers
    if not inps.meta_file or inps.meta_file == 'auto':
        print('No meta file provided. Generating burst XML from static_layers.h5...')
        geom_path = Path(inps.geom_dir)
        # Find all subdirectories containing static_layers*.h5
        burst_dirs = sorted([d for d in geom_path.iterdir()
                             if d.is_dir() and list(d.glob('static_layers*.h5'))])
        # Also check for HDF5 files directly in geom_dir
        h5_files_in_root = sorted(geom_path.glob('static_layers*.h5'))
        if not burst_dirs and not h5_files_in_root:
            raise FileNotFoundError(f'No static_layers HDF5 found in {inps.geom_dir}')

        if burst_dirs:
            first_burst_dir = burst_dirs[0]
            h5_list = sorted(first_burst_dir.glob('static_layers*.h5'))
            first_h5 = h5_list[0]
            burst_id = first_burst_dir.name
        else:
            first_h5 = h5_files_in_root[0]
            burst_id = os.path.basename(geom_path)

        xml_out = Path(inps.out_dir) / f'{burst_id}.burst.xml'
        xml_out.parent.mkdir(parents=True, exist_ok=True)

        isce3_utils.generate_burst_xml_from_static(str(first_h5), str(xml_out))
        inps.meta_file = str(xml_out)
        print(f'Generated meta file: {inps.meta_file}')

    # Read common metadata from reference burst XML (now guaranteed to exist)
    metadata = {}
    if inps.meta_file:
        metadata = isce3_utils.extract_isce3_metadata(
            inps.meta_file,
            update_mode=inps.update_mode
        )

    # Determine target shape from first interferogram if available
    target_shape = None
    ref_int = None
    if inps.obs_files:
        int_list = glob.glob(inps.obs_files[0])
        if int_list:
            ref_int = int_list[0]
            int_atr = readfile.read_attribute(ref_int)
            target_shape = (int(int_atr['LENGTH']), int(int_atr['WIDTH']))
            print(f'Target shape from interferogram: {target_shape}')

    # Prepare geometry (updates metadata with LENGTH/WIDTH)
    if inps.geom_dir:
        metadata = prepare_geometry_isce3(
            geom_dir=inps.geom_dir,
            out_dir=inps.out_dir,
            geom_files=inps.geom_files,
            metadata=metadata,
            processor=inps.processor,
            update_mode=inps.update_mode,
            ref_int_file=ref_int,
            target_shape=target_shape,
            geom_dirs=getattr(inps, 'geom_dirs', None)
        )

    # Read baseline info
    baseline_dict = {}
    if inps.baseline_dir:
        baseline_dict = isce3_utils.read_baseline_timeseries_isce3(
            inps.baseline_dir,
            processor=inps.processor
        )

    # Prepare metadata for interferogram stack(s)
    if inps.obs_files:
        for obs_file in inps.obs_files:
            prepare_stack_isce3(
                obs_file,
                metadata=metadata,
                baseline_dict=baseline_dict,
                update_mode=inps.update_mode
            )

    print('Done.')
    return