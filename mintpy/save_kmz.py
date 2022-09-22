############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################
# recommended import
#   from mintpy import save_kmz


import os
import shutil
from zipfile import ZipFile

import numpy as np
from lxml import etree
from matplotlib import colorbar, colors, pyplot as plt, ticker
from pykml.factory import KML_ElementMaker as KML

import mintpy
from mintpy import subset
from mintpy.utils import (
    attribute as attr,
    plot as pp,
    ptime,
    readfile,
    utils as ut,
)


############################################################
def plot_colorbar(out_file, vmin, vmax, unit='cm/year', cmap='jet', figsize=(0.18, 3.6),
                  nbins=None, label='Mean LOS velocity'):
    fig, cax = plt.subplots(figsize=figsize)
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cbar = colorbar.ColorbarBase(cax, cmap=plt.get_cmap(cmap), norm=norm, orientation='vertical')
    cbar.set_label(f'{label} [{unit}]', fontsize=12)

    # update ticks
    if nbins:
        cbar.locator = ticker.MaxNLocator(nbins=nbins)
        cbar.update_ticks()

    cbar.ax.tick_params(which='both', labelsize=12)
    fig.patch.set_facecolor('white')
    fig.patch.set_alpha(0.7)
    print('writing', out_file)
    fig.savefig(out_file, bbox_inches='tight', facecolor=fig.get_facecolor(), dpi=300)
    return out_file


def generate_cbar_element(cbar_file, cmap, vmin, vmax, unit='cm/year', loc='lower left',
                          nbins=None, label='Mean LOS velocity'):
    """Generate colorbar as an screen overlay object.

    Parameters: cbar_file - str, colorbar image file path
                cmap      - matplotlib.colors.Colormap instance
                vmin/vmax - float, min/max value to display
                unit      - str, display unit
                loc       - str, location of colorbar on the screen.
                            lower-left, lower-right, upper-left, upper-right
    Returns:    cbar_overlay - KML.ScreenOverlay object
    """
    # plot colobar and save as an image
    cbar_file = plot_colorbar(
        out_file=cbar_file,
        vmin=vmin,
        vmax=vmax,
        unit=unit,
        cmap=cmap,
        nbins=nbins,
        label=label)

    # colobar location
    if loc.split()[0] == 'lower':
        oy, sy = '0', '0'
    elif loc.split()[0] == 'upper':
        oy, sy = '1', '1'
    if loc.split()[1] == 'left':
        ox, sx = '0', '0'
    elif loc.split()[1] == 'right':
        ox, sx = '1', '1'

    # generate KML screen overlay object
    cbar_overlay = KML.ScreenOverlay(
        KML.name('colorbar'),
        KML.Icon(
            KML.href(f"{os.path.basename(cbar_file)}"),
            KML.viewBoundScale(0.75)
        ),
        KML.overlayXY(x=ox, y=oy, xunits="fraction", yunits="fraction"),
        KML.screenXY(x=sx, y=sy, xunits="fraction", yunits="fraction"),
        KML.size(x="0", y="250", xunits="pixel", yunits="pixel"),
        KML.rotation(0),
        KML.visibility(1),
        KML.open(0)
    )
    #print('add colorbar.')
    return cbar_overlay


def get_hex_color(v, colormap, norm):
    """Get color name in hex format.
    Parameters: v        : float, number of interest
                colormap : matplotlib.colors.Colormap instance
                norm     : matplotlib.colors.Normalize instance
    Returns:    c_hex    : color name in hex format
    """
    # get rgba color components for point velocity
    rgba = colormap(norm(v))
    # rgba to hex
    c_hex = colors.to_hex([rgba[3], rgba[2], rgba[1], rgba[0]], keep_alpha=True)[1:]
    return c_hex


def create_placemark_element(lat, lon, row, col, val, icon_file, inps):
    """Create an KMZ Placemark element.
    Parameters: lat/lon   - float, latitude / longitude in degrees
                row/col   - int, row / column number
                val       - float, value
                icon_file - str, path of the icon file
                inps      - Namespace object
    Returns:    placemark - KMZ.Placemark() object
    """

    placemark = KML.Placemark(
        # style
        KML.Style(
            KML.IconStyle(
                    KML.color(get_hex_color(val, inps.colormap, inps.norm)),
                    KML.scale(0.5),
                    KML.Icon(KML.href(icon_file),
                ),
            ),
        ),

        # extended data
        KML.ExtendedData(
            KML.Data(
                KML.value(f"{lat:.6f}˚"),
                name="Latitude",
            ),
            KML.Data(
                KML.value(f"{lon:.6f}˚"),
                name="Longitude",
            ),
            KML.Data(
                KML.value(f"{row:.0f}"),
                name="Row",
            ),
            KML.Data(
                KML.value(f"{col:.0f}"),
                name="Column",
            ),
            KML.Data(
                KML.value(f"{val:.2f} {inps.disp_unit}"),
                name="Value",
            ),
        ),

        # point coord
        KML.Point(
            KML.coordinates(f"{lon},{lat}"),
        ),
    )

    return placemark


def write_kmz_file(out_file_base, kml_doc, data_files=None, res_files=None, keep_kml_file=False):
    """Write KML and KMZ files.
    Parameters: out_file_base - str, output file name without extension
                kml_doc       - KML.Document() object
                data_files    - list of str, rel path of data files
                res_files     - list of str, rel path of resource files
                keep_kml_file - bool, do not remove KML files after zipping.
    Returns:    kmz_file      - str, zipped KMZ file.
    """
    # default values
    data_files = [] if data_files is None else data_files
    res_files  = [] if res_files  is None else res_files

    work_dir = os.path.dirname(out_file_base)
    kml_file = f'{out_file_base}.kml'
    kmz_file = f'{out_file_base}.kmz'

    # 1. Write KML file
    kml = KML.kml()
    kml.append(kml_doc)

    print('writing '+kml_file)
    with open(kml_file, 'w') as f:
        f.write(etree.tostring(kml, pretty_print=True).decode('utf8'))

    # 2. Copy resource files
    if res_files:
        res_dir = os.path.join(os.path.dirname(mintpy.__file__), "data")
        for fname in res_files:
            src_file = os.path.join(res_dir, os.path.basename(fname))
            shutil.copy2(src_file, work_dir)
            print(f"copy {src_file} to the local directory")

    # 3. Generate KMZ file, by
    # 1) go to the directory of kmz file
    run_dir = os.path.abspath(os.getcwd())
    os.chdir(work_dir)

    # 2) zip all data files
    with ZipFile(kmz_file, 'w') as fz:
        for fname in [kml_file] + data_files + res_files:
            fz.write(os.path.relpath(fname))
            if not keep_kml_file:
                os.remove(fname)
                print(f'remove {fname}')

    # 3) go back to the running directory
    os.chdir(run_dir)
    print(f'merged all files to {kmz_file}')

    return kmz_file


############################################################
def write_kmz_overlay(data, meta, out_file, inps):
    """Generate Google Earth Overlay KMZ file for data in GEO coordinates.
    Parameters: data     - 2D np.array in int/float, data matrix to write
                meta     - dict, containing the following attributes:
                           WIDTH/LENGTH      : required, file size
                           X/Y_FIRST/STEP    : required, for lat/lon spatial converage
                           REF_X/Y           : optional, column/row number of reference pixel
                out_file - string, output file name
                inps     - Namespace, optional, input options for display
    Returns:    kmz_file - string, output KMZ filename
    """

    south, north, west, east = ut.four_corners(meta)

    # 1. Make PNG file - Data
    print('plotting data ...')

    # Figure size
    if not inps.fig_size:
        inps.fig_size = pp.auto_figure_size(ds_shape=[north-south, east-west], scale=2.0)
    fig = plt.figure(figsize=inps.fig_size, frameon=False)
    ax = fig.add_axes([0., 0., 1., 1.])
    ax.set_axis_off()

    # Plot - data matrix
    ax.imshow(data, vmin=inps.vlim[0], vmax=inps.vlim[1], cmap=inps.colormap,
              aspect='auto', interpolation='nearest')

    # Plot - reference pixel
    rx = meta.get('REF_X', None)
    ry = meta.get('REF_Y', None)
    if inps.disp_ref_pixel and rx is not None and ry is not None:
        ax.plot(int(rx), int(ry), inps.ref_marker,
                color=inps.ref_marker_color,
                ms=inps.ref_marker_size)
        print('show reference point')
    else:
        print('no plot for reference point.')

    width = int(meta['WIDTH'])
    length = int(meta['LENGTH'])
    ax.set_xlim([0, width])
    ax.set_ylim([length, 0])

    out_file_base = os.path.splitext(out_file)[0]
    data_png_file = out_file_base + '.png'
    print(f'writing {data_png_file} with dpi={inps.fig_dpi}')
    plt.savefig(data_png_file, pad_inches=0.0, transparent=True, dpi=inps.fig_dpi)

    # 2. Generate KML file
    kml_doc = KML.Document()

    # Add data png file
    img_name = os.path.splitext(os.path.basename(data_png_file))[0]
    img_overlay = KML.GroundOverlay(
        KML.name(img_name),
        KML.Icon(
            KML.href(os.path.basename(data_png_file))
        ),
        KML.altitudeMode('clampToGround'),
        KML.LatLonBox(
            KML.north(str(north)),
            KML.east(str(east)),
            KML.south(str(south)),
            KML.west(str(west)),
        ),
    )
    kml_doc.append(img_overlay)

    # Add colorbar png file
    cbar_file = f'{out_file_base}_cbar.png'
    cbar_overlay = generate_cbar_element(
        cbar_file,
        cmap=inps.colormap,
        vmin=inps.vlim[0],
        vmax=inps.vlim[1],
        unit=inps.disp_unit,
        loc=inps.cbar_loc,
        nbins=inps.cbar_bin_num,
        label=inps.cbar_label)
    kml_doc.append(cbar_overlay)

    # Write KML file
    kmz_file = write_kmz_file(
        out_file_base,
        kml_doc,
        data_files=[data_png_file, cbar_file],
        keep_kml_file=inps.keep_kml_file)

    return kmz_file


def write_kmz_placemark(data, meta, out_file, geom_file, inps):
    """Generate Google Earth Placemark KMZ file for data in RADAR coordinates.
    Parameters: data      - 2D np.array in int/float, data matrix to write
                meta      - dict, containing the following attributes:
                            WIDTH/LENGTH      : required, file size
                            X/Y_FIRST/STEP    : required, for lat/lon spatial converage
                            REF_X/Y           : optional, column/row number of reference pixel
                geom_file - str, path of the geometry file with latitude/longitude datasets
                out_file  - string, output file name
                inps      - Namespace, optional, input options for display
    Returns:    kmz_file  - string, output KMZ filename
    """

    out_file_base = os.path.splitext(out_file)[0]
    dot_file = 'shaded_dot.png'
    star_file = 'star.png'

    # read latitude / longitude
    lats = readfile.read(geom_file, datasetName='latitude',  box=inps.pix_box)[0]
    lons = readfile.read(geom_file, datasetName='longitude', box=inps.pix_box)[0]

    ## Generate KML file
    kml_doc = KML.Document()

    # 1. colorbar png file
    print('plot and add colorbar as a ScreenOverlay element')
    cbar_file = f'{out_file_base}_cbar.png'
    cbar_overlay = generate_cbar_element(
        cbar_file,
        cmap=inps.colormap,
        vmin=inps.vlim[0],
        vmax=inps.vlim[1],
        unit=inps.disp_unit,
        loc=inps.cbar_loc,
        nbins=inps.cbar_bin_num,
        label=inps.cbar_label)
    kml_doc.append(cbar_overlay)

    # 2. reference point
    xmin = int(meta.get('SUBSET_XMIN', 0))
    ymin = int(meta.get('SUBSET_YMIN', 0))

    if 'REF_Y' in meta.keys():
        print('add reference point as a star icon')
        ry, rx = int(meta['REF_Y']), int(meta['REF_X'])
        rlat = lats[ry, rx]
        rlon = lons[ry, rx]
        ref_point = create_placemark_element(
            lat=rlat,
            lon=rlon,
            row=ry + ymin,
            col=rx + xmin,
            val=0.0,
            icon_file=star_file,
            inps=inps)
        ref_point.name = 'ReferencePoint'
        ref_point.Style.IconStyle.scale = 1.0
        kml_doc.append(ref_point)

        # do not plot reference point as data again
        data[ry, rx] = np.nan

    # 3. data folder for all points
    data_folder = KML.Folder(KML.name("Data"))

    print(f'generating point element with step size of {inps.step} pixels')
    length, width = data.shape
    prog_bar = ptime.progressBar(maxValue=length)
    for y in range(0, length, inps.step):
        for x in range(0, width, inps.step):
            value = data[y, x]
            if not np.isnan(value):
                lat = lats[y, x]
                lon = lons[y, x]

                # create KML icon element
                placemark = create_placemark_element(
                    lat=lat,
                    lon=lon,
                    row=y + ymin,
                    col=x + xmin,
                    val=value,
                    icon_file=dot_file,
                    inps=inps)
                data_folder.append(placemark)

        prog_bar.update(y+1, every=1, suffix=f'row={y+1}/{length}')
    prog_bar.close()
    kml_doc.append(data_folder)

    # Write KML file
    kmz_file = write_kmz_file(
        out_file_base,
        kml_doc,
        data_files=[cbar_file],
        res_files=[dot_file, star_file],
        keep_kml_file=inps.keep_kml_file)

    return kmz_file


############################################################
def save_kmz(inps):

    # matplotlib backend setting
    plt.switch_backend('Agg')

    ## 1. Read metadata and data
    ftype = readfile.read_attribute(inps.file)['FILE_TYPE']
    if ftype == 'timeseries' and inps.dset and '_' in inps.dset:
        inps.ref_date, inps.dset = inps.dset.split('_')
    else:
        inps.ref_date = None
    atr = readfile.read_attribute(inps.file, datasetName=inps.dset)

    # pix_box
    inps.pix_box = subset.subset_input_dict2box(vars(inps), atr)[0]
    inps.pix_box = ut.coordinate(atr).check_box_within_data_coverage(inps.pix_box)
    data_box = (0, 0, int(atr['WIDTH']), int(atr['LENGTH']))
    print(f'data   coverage in y/x: {data_box}')
    print(f'subset coverage in y/x: {inps.pix_box}')
    atr = attr.update_attribute4subset(atr, inps.pix_box)

    # read data
    data = readfile.read(inps.file, datasetName=inps.dset, box=inps.pix_box)[0]
    if ftype == 'timeseries' and inps.ref_date:
        data -= readfile.read(inps.file, datasetName=inps.ref_date, box=inps.pix_box)[0]

    # mask
    mask = pp.read_mask(
        inps.file,
        mask_file=inps.mask_file,
        datasetName=inps.dset,
        box=inps.pix_box)[0]
    if mask is not None:
        print(f'masking out pixels with zero value in file: {inps.mask_file}')
        data[mask == 0] = np.nan
    if inps.zero_mask:
        print('masking out pixels with zero value')
        data[data == 0] = np.nan
    del mask

    # Data Operation - Display Unit & Rewrapping
    data, inps.disp_unit, inps.disp_scale, inps.wrap = pp.scale_data4disp_unit_and_rewrap(
        data,
        metadata=atr,
        disp_unit=inps.disp_unit,
        wrap=inps.wrap,
        wrap_range=inps.wrap_range,
    )
    if inps.wrap:
        inps.vlim = inps.wrap_range


    ## 2. Generate Google Earth KMZ
    # disp min/max and colormap
    cmap_lut = 256
    if not inps.vlim:
        cmap_lut, inps.vlim = pp.auto_adjust_colormap_lut_and_disp_limit(data)[:2]
    inps.cmap_name = pp.auto_colormap_name(atr, inps.cmap_name)
    inps.colormap = pp.ColormapExt(inps.cmap_name, cmap_lut).colormap
    inps.norm = colors.Normalize(vmin=inps.vlim[0], vmax=inps.vlim[1])

    # output filename
    inps.fig_title = pp.auto_figure_title(inps.file, inps.dset, vars(inps))
    inps.outfile = inps.outfile if inps.outfile else f'{inps.fig_title}.kmz'
    inps.outfile = os.path.abspath(inps.outfile)

    # write KMZ file
    if 'Y_FIRST' in atr.keys():
        # create ground overlay KML for file in geo-coord
        write_kmz_overlay(
            data,
            meta=atr,
            out_file=inps.outfile,
            inps=inps)

    else:
        # create placemark KML for file in radar-coord
        write_kmz_placemark(
            data,
            meta=atr,
            out_file=inps.outfile,
            geom_file=inps.geom_file,
            inps=inps)

    return
