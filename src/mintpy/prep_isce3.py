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
                           target_shape=None):
    """Prepare geometry files from ISCE3/Dolphin static_layers HDF5.

    Parameters
    ----------
    target_shape : tuple of (length, width), optional
        Interferogram dimensions: used for multilooking geometry and
        computing ALOOKS/RLOOKS from burst full-resolution dimensions.
    """
    from pathlib import Path

    print('preparing geometry files from ISCE3/Dolphin static layers')
    geom_dir = os.path.abspath(geom_dir)
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    if geom_files is None:
        geom_files = ['height.tif', 'los_east.tif', 'los_north.tif',
                      'layover_shadow_mask.tif', 'local_incidence_angle.tif']

    # Step 0: Read burst count and info from static_layers.h5
    num_bursts = 0
    geom_path = Path(geom_dir)
    burst_subdirs = sorted([d for d in geom_path.iterdir()
                            if d.is_dir() and list(d.glob('static_layers*.h5'))])
    num_bursts = len(burst_subdirs)
    if burst_subdirs:
        print(f'Number of bursts: {num_bursts}')

    # Step 1: Merge full-resolution geometry WITHOUT cropping to interferogram extent
    geometry_dict = isce3_utils.extract_merge_geometry(
        geom_dir=geom_dir,
        output_dir=out_dir,
        geom_types=geom_files,
        ref_int_file=None,
        metadata=metadata
    )

    # Step 2: Read merged full-resolution geometry dimensions
    height_file = geometry_dict.get('height.tif')
    if height_file and os.path.isfile(str(height_file)):
        geom_atr = readfile.read_attribute(str(height_file))
        geom_shape = (int(geom_atr['LENGTH']), int(geom_atr['WIDTH']))
    else:
        geom_shape = None

    # Step 3: Multilook geometry to match interferogram resolution
    lks_y = 1
    lks_x = 1
    if target_shape is not None and geom_shape is not None and geom_shape != target_shape:
        if (geom_shape[0] % target_shape[0] != 0 or
            geom_shape[1] % target_shape[1] != 0):
            raise ValueError(
                f'Geometry shape {geom_shape} is not an integer multiple of '
                f'target interferogram shape {target_shape}.')
        lks_y = geom_shape[0] // target_shape[0]
        lks_x = geom_shape[1] // target_shape[1]
        if lks_y > 1 or lks_x > 1:
            print(f'Multilooking geometry from {geom_shape} to {target_shape} '
                  f'(looks: {lks_y} x {lks_x}) and overwriting originals')
            geometry_dict = isce3_utils.multilook_geometry_files(
                geometry_dict, lks_y, lks_x, output_dir=None, overwrite=True)

    # Step 4: Crop geometry to match interferogram extent
    if ref_int_file and os.path.isfile(ref_int_file):
        geometry_dict = isce3_utils.crop_geometry_to_reference(
            geometry_dict, ref_int_file, output_dir=None, overwrite=True)

    # Step 5: Set ALOOKS/RLOOKS from multilook factors
    if metadata is not None:
        metadata['ALOOKS'] = str(lks_y)
        metadata['RLOOKS'] = str(lks_x)
        msg = f'Set ALOOKS={lks_y}, RLOOKS={lks_x}'
        if geom_shape is not None:
            msg += f' (geom {geom_shape} -> target {target_shape})'
        print(msg)

    # Write .rsc files (after possible multilook)
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
        if not burst_dirs:
            raise FileNotFoundError(f'No static_layers HDF5 found in {inps.geom_dir}')

        first_burst_dir = burst_dirs[0]
        # Get the first HDF5 file in the first burst directory
        h5_list = sorted(first_burst_dir.glob('static_layers*.h5'))
        if not h5_list:
            raise FileNotFoundError(f'No static_layers HDF5 found in {first_burst_dir}')
        first_h5 = h5_list[0]

        # Use the burst directory name as burst ID
        burst_id = first_burst_dir.name
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
            target_shape=target_shape
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