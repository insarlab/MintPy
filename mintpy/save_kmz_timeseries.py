############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Joshua Zahner, Zhang Yunjun, 2018                #
############################################################


import os
import shutil
from zipfile import ZipFile

import matplotlib as mpl
import numpy as np
from lxml import etree
from matplotlib import patches, pyplot as plt
from pykml.factory import KML_ElementMaker as KML

import mintpy
from mintpy import save_kmz
from mintpy.objects import deramp, timeseries
from mintpy.utils import plot as pp, readfile, utils as ut


############################################################
def get_all_file_paths(directory):
    # initializing empty file paths list
    file_paths = []

    # crawling through directory and subdirectories
    for root, _, files in os.walk(directory):
        for filename in files:
            # join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)

    # returning all file paths
    return file_paths


def flatten_lat_lon(box, ts_obj, coords=None):
    if coords is None:
        lats, lons = ut.get_lat_lon(ts_obj.metadata, box=box)
    lats = sorted(lats.flatten())
    lons = sorted(lons.flatten())
    return lats, lons


def split_into_sub_boxes(ds_shape, step=20, num_pixel=50**2, print_msg=True):
    """split input shape into multiple sub boxes
    Parameters: ds_shape  : 2-tuple of int for the shape of whole dataset in (length, width)
                step      : int, number of pixels per selection
                num_pixel : int, max number of selections per box
    Returns:    box_list  : list of 4-tuple of int, each indicating (col0, row0, col1, row1)
    """
    print('-'*30)
    print(f'step: {step} pixels')
    win_size = int(np.sqrt(num_pixel) * step)
    length, width = ds_shape

    # number of sub boxes
    nrows = np.ceil(length / win_size).astype(int)
    ncols = np.ceil(width / win_size).astype(int)
    if print_msg:
        print(f'number of output boxes: {nrows * ncols}')

    # start/end row/column number of each box
    box_list = []
    for i in range(nrows):
        r0 = i * win_size
        r1 = min([length, r0 + win_size])
        for j in range(ncols):
            c0 = j * win_size
            c1 = min([width, c0 + win_size])
            box = (c0, r0, c1, r1)
            box_list.append(box)
    return box_list


def get_boxes4deforming_area(vel_file, mask_file, step=2, num_pixel=30**2, min_perc=0.2,
                             cutoff=3, ramp_type='quadratic', display=False):
    """Get list of boxes to cover the deforming areas.
    A pixel is identified as deforming if its velocity exceeds the MAD of the whole image.
    Parameters: vel_file  - str, path of velocity file
                mask_file - str, path of mask file
                win_size  - int, length and width of the output box
                min_perc  - float between 0 and 1, minimum percentage of deforming points in the box
                ramp_type - str, type of phase ramps to be removed while evaluating the deformation
                display   - bool, plot the identification result or not
    Returns:    box_list  - list of t-tuple of int, each indicating (col0, row0, col1, row1)
    """
    win_size = int(np.sqrt(num_pixel) * step)
    print('-'*30)
    print(f'get boxes on deforming areas with step: {step} pixels')
    mask = readfile.read(mask_file)[0]
    vel, atr = readfile.read(vel_file)
    print(f'removing a {ramp_type} phase ramp from input velocity before the evaluation')
    vel = deramp(vel, mask, ramp_type=ramp_type, metadata=atr)[0]               #remove ramp before the evaluation

    # get deforming pixels
    mad = ut.median_abs_deviation_threshold(vel[mask], cutoff=cutoff)     #deformation threshold
    print(f'velocity threshold / median abs dev: {mad:.3f} cm/yr')
    vel[mask == 0] = 0
    mask_aoi = (vel >= mad) + (vel <= -1. * mad)
    print(f'number of points: {np.sum(mask_aoi)}')

    # get deforming boxes
    box_list = []
    min_num = min_perc * (win_size ** 2)
    length, width = vel.shape
    num_row = np.ceil(length / win_size).astype(int)
    num_col = np.ceil(width / win_size).astype(int)
    for i in range(num_row):
        r0 = i * win_size
        r1 = min([length, r0 + win_size])
        for j in range(num_col):
            c0 = j * win_size
            c1 = min([width, c0 + win_size])
            box = (c0, r0, c1, r1)
            if np.sum(mask_aoi[r0:r1, c0:c1]) >= min_num:
                box_list.append(box)
    print(f'number of boxes : {len(box_list)}')

    if display:
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=[8, 4], sharey=True)
        vel[mask == 0] = np.nan
        axs[0].imshow(vel, cmap='jet')
        axs[1].imshow(mask_aoi, cmap='gray_r')
        for box in box_list:
            for ax in axs:
                rect = patches.Rectangle((box[0], box[1]),
                                         width=(box[2]-box[0]),
                                         height=(box[3]-box[1]),
                                         linewidth=2, edgecolor='r', fill=False)
                ax.add_patch(rect)
        fig.tight_layout()
        out_fig = os.path.join(os.path.dirname(vel_file), 'defo_area.png')
        fig.savefig(out_fig, bbox_inches='tight', transparent=True, dpi=300)
        print(f'save figure to {out_fig}')
        plt.show()
    return box_list


def create_reference_point_element(inps, lats, lons, ts_obj):
    """Create reference point element"""
    norm = mpl.colors.Normalize(vmin=inps.vlim[0], vmax=inps.vlim[1])

    ref_yx = (int(ts_obj.metadata['REF_Y']), int(ts_obj.metadata['REF_X']))
    ref_lalo = (lats[ref_yx[0], ref_yx[1]], lons[ref_yx[0], ref_yx[1]])

    ref_point = KML.Placemark(
        KML.Style(
            KML.IconStyle(
                KML.color(save_kmz.get_hex_color(0.0, inps.colormap, norm)),
                KML.scale(1.),
                KML.Icon(
                    KML.href(f"{os.path.basename(inps.star_file)}")
                )
            )
        ),
        KML.description("Reference point <br /> \n <br /> \n" +
                        get_description_string(ref_lalo, ref_yx, 0.00, 0.00, 0.00, 1.00)
                        ),
        KML.Point(
            KML.coordinates(f"{ref_lalo[1]}, {ref_lalo[0]}")
        )
    )

    return ref_point


def get_description_string(coords, yx, v, vstd, disp, tcoh=None, font_size=4):
    """Description information of each data point."""
    des_str = f"<font size={font_size}>"
    des_str += f"Latitude: {coords[0]:.6f}˚ <br /> \n"
    des_str += f"Longitude: {coords[1]:.6f}˚ <br /> \n"
    des_str += f"Row: {yx[0]:.0f} <br /> \n"
    des_str += f"Column: {yx[1]:.0f} <br /> \n"
    des_str += " <br /> \n"
    des_str += f"Mean LOS velocity [cm/year]: {v:.2f} +/- {vstd:.2f} <br /> \n"
    des_str += f"Cumulative displacement [cm]: {disp:.2f} <br /> \n"
    if tcoh is not None:
        des_str += f"Temporal coherence: {tcoh:.2f} <br /> \n"
    des_str += "</font>"
    des_str += " <br />  <br /> "
    des_str += "*Double click to reset plot <br /> <br />\n"
    des_str += "\n\n"
    return des_str


def generate_js_datastring(dates, dygraph_file, num_date, ts):
    """String of the Java Script for interactive plot of diplacement time-series"""
    dygraph_file = f'../../../{os.path.basename(dygraph_file)}'
    js_data_string = f"<script type='text/javascript' src='{dygraph_file}'></script>"
    js_data_string += """
        <div id='graphdiv'> </div>
        <style>
            .dygraph-legend{
                left: 230px !important;
                width: 265px !important;
            }
        </style>
        <script type='text/javascript'>
            g = new Dygraph( document.getElementById('graphdiv'),
            "Date, displacement\\n" +
    """

    # append the date/displacement data
    for k in range(num_date):
        date = dates[k]
        dis = ts[k]
        date_displacement_string = f"\"{date}, {dis}\\n\" + \n"
        js_data_string += date_displacement_string

    js_data_string += """

    "",
       {
         width: 500,
         height: 300,
         axes: {
             x: {
                 axisLabelFormatter: function (d, gran) {
                     var months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
                     var date = new Date(d)
                     var dateString = months[date.getMonth()] + ' ' + date.getFullYear()
                     return dateString;
                 },
                 valueFormatter: function (d) {
                     var date = new Date(d)
                     var dateString = 'Date: ' + ('0' + date.getDate()).slice(-2) +
                                      '/' + ('0' + (date.getMonth() + 1)).slice(-2) +
                                      '/' + date.getFullYear()
                     return dateString;
                 },
                 pixelsPerLabel: 90
             },
             y: {
                 valueFormatter: function(v) {
                    if(v >= 0){
                        return (' 000' + v.toFixed(2)).slice(-8)
                    }else{
                        return '-'+(' 000' + Math.abs(v).toFixed(2)).slice(-7)
                    }
                 }
             }
         },
         ylabel: 'LOS displacement [cm]',
         yLabelWidth: 18,
         drawPoints: true,
         strokeWidth: 0,
         pointSize: 3,
         highlightCircleSize: 6,
         axisLabelFontSize: 12,
         xRangePad: 30,
         yRangePad: 30,
         hideOverlayOnMouseOut: false,
         panEdgeFraction: 0.0
       });
       </script>

    """
    return js_data_string


def create_kml_region_document(inps, box_list, ts_obj, step):
    """Create list of KML.Document() objects
    for one level of details defined by box_list and step
    """
    dot_file = f'../../{os.path.basename(inps.dot_file)}'

    ## 1. read data file into timeseries object
    region_docs = []
    num_box = len(box_list)
    for i in range(num_box):
        box = box_list[i]

        if box is None:
            box = (0, 0, ts_obj.width, ts_obj.length)

        length = box[3] - box[1]
        width = box[2] - box[0]

        # 1.1 Parse Date
        dates = np.array(ts_obj.times)  # 1D np.array of dates in datetime.datetime object in size of [num_date,]
        dates = list(map(lambda d: d.strftime("%Y-%m-%d"), dates))
        num_date = len(dates)

        # 1.2 Parse Spatial coordinates
        lats, lons = ut.get_lat_lon(ts_obj.metadata, box=box)
        rows, cols = np.mgrid[box[1]:box[3] - 1:length * 1j, box[0]:box[2] - 1:width * 1j]

        # 1.3 Read Velocity / time-series / temporal coherence data
        vel = readfile.read(inps.vel_file, datasetName='velocity', box=box)[0] * 100.
        vel_std = readfile.read(inps.vel_file, datasetName='velocityStd', box=box)[0] * 100.
        ts_data = readfile.read(inps.ts_file, box=box)[0] * 100.
        ts_data -= np.tile(ts_data[0, :, :], (ts_data.shape[0], 1, 1))  # enforce displacement starts from zero
        temp_coh = readfile.read(inps.tcoh_file, box=box)[0]
        mask = readfile.read(inps.mask_file, box=box)[0]

        vel_c = np.array(vel, dtype=np.float32)
        if inps.wrap:
            vel_c = inps.vlim[0] + np.mod(vel_c - inps.vlim[0], inps.vlim[1] - inps.vlim[0])


        ## 2. Create KML Document
        kml_doc = KML.Document()

        # 2.1 Normalize colormap to defined vlim
        norm = mpl.colors.Normalize(vmin=inps.vlim[0], vmax=inps.vlim[1])

        # 2.2 Set number of pixels to use
        num_pixel = int(length / step) * int(width / step)
        msg = f"create KML doc for box {i+1}/{num_box}: {box}"
        msg += f", step: {step} pixels, {num_pixel} pixels in total ..."
        print(msg)

        # 2.3 Create data folder for all points
        data_folder = KML.Folder(KML.name("Data"))
        for i in range(0, length, step):
            for j in range(0, width, step):
                if mask[i, j]:  # add point if it's not marked as masked out
                    lat = lats[i, j]
                    lon = lons[i, j]
                    row = rows[i, j]
                    col = cols[i, j]
                    ts = ts_data[:, i, j]
                    v = vel[i, j]
                    vc = vel_c[i, j]
                    vstd = vel_std[i, j]
                    tcoh = temp_coh[i, j]

                    # 2.3.1 Create KML icon style element
                    style = KML.Style(
                        KML.IconStyle(
                            KML.color(save_kmz.get_hex_color(vc, inps.colormap, norm)),
                            KML.scale(0.5),
                            KML.Icon(KML.href(f"{dot_file}"))
                        )
                    )

                    # 2.3.2 Create KML point element
                    point = KML.Point(KML.coordinates(f"{lon},{lat}"))

                    js_data_string = generate_js_datastring(dates, inps.dygraph_file, num_date, ts)

                    # 2.3.3 Create KML description element
                    stats_info = get_description_string((lat, lon), (row, col), v, vstd, ts[-1], tcoh=tcoh)
                    description = KML.description(stats_info, js_data_string)

                    # 2.3.4 Crate KML Placemark element to hold style, description, and point elements
                    placemark = KML.Placemark(style, description, point)

                    # 2.3.5 Append each placemark element to the KML document object
                    data_folder.append(placemark)

        # 2.4 Append data folder to KML document
        kml_doc.append(data_folder)

        # 2.5 Append KML document to list of regionalized documents
        region_docs.append(kml_doc)

    return region_docs


def write_network_link_file(region_docs, ts_obj, box_list, lod, net_link_file):
    """Write 1) the list of KML.Document() into data KML file and 2) the root kml file for the list"""

    ## 1. Create directory to store regioalized KML data files
    links_dir = os.path.splitext(net_link_file)[0]
    if not os.path.isdir(links_dir):
        os.makedirs(links_dir)
    print(f"create KML region links directory: {os.path.basename(links_dir)}")

    ## 2. Create root KML element and KML Document element
    kml = KML.kml()
    kml_doc = KML.Document()

    ## 3. Generate a new network link element for each region
    for num, (region_doc, box) in enumerate(zip(region_docs, box_list)):
        region_kml_file = os.path.join(links_dir, f"region_{num}.kml")

        ## 3.1 Write the first region_document to a file and move it to the proper subdircetory
        kml_1 = KML.kml()
        kml_1.append(region_doc)
        with open(region_kml_file, 'w') as f:
            f.write(etree.tostring(kml_1, pretty_print=True).decode('utf-8'))

        ## 3.2 Flatten lats and lons data
        lats, lons = flatten_lat_lon(box, ts_obj)

        ## 3.3 Define new NetworkLink element
        network_link = KML.NetworkLink(
            KML.name(f'Region {num}'),
            KML.visibility(1),
            KML.Region(
                KML.Lod(
                    KML.minLodPixels(lod[0]),
                    KML.maxLodPixels(lod[1])
                ),
                KML.LatLonAltBox(
                    KML.north(lats[0] + 0.5),
                    KML.south(lats[-1] - 0.5),
                    KML.east(lons[0] + 0.5),
                    KML.west(lons[-1] - 0.5)
                )
            ),
            KML.Link(
                KML.href(os.path.relpath(region_kml_file, start=os.path.dirname(links_dir))),
                KML.viewRefreshMode('onRegion')
            )
        )

        ## 3.4 Append new NetworkLink to KML document
        kml_doc.append(network_link)
    kml.append(kml_doc)

    ## 4. Write the full KML document to the output file and move it to the proper directory
    with open(net_link_file, 'w') as f:
        f.write(etree.tostring(kml, pretty_print=True).decode('utf-8'))
    return net_link_file


def create_network_link_element(net_link_file, ts_obj):
    """Create an KML.NetworkLink element for one level of details"""
    net_link_name = os.path.splitext(os.path.basename(net_link_file))[0]

    lats, lons = flatten_lat_lon(None, ts_obj)

    network_link = KML.NetworkLink(
        KML.name(net_link_name),
        KML.visibility(1),
        KML.Region(
            KML.Lod(
                KML.minLodPixels(0),
                KML.maxLodPixels(1800)
            ),
            KML.LatLonAltBox(
                KML.north(lats[-1] + 0.5),
                KML.south(lats[0] - 0.5),
                KML.east(lons[-1] + 0.5),
                KML.west(lons[0] - 0.5)
            )
        ),
        KML.Link(
            KML.href(net_link_file),
            KML.viewRefreshMode('onRegion')
        )
    )
    return network_link


def generate_network_link(inps, ts_obj, step, lod):
    """Generate the KML.NetworkLink element for one level of details, defined by step and lod"""

    net_link_file = os.path.join(inps.kml_data_dir, "{0}by{0}.kml".format(step))

    if step <= 0:
        print(f'skip step = {step}')
        return None

    elif step == inps.steps[-1]:
        box_list = get_boxes4deforming_area(
            inps.vel_file, inps.mask_file,
            step=inps.steps[-1],
            min_perc=inps.min_percentage,
            cutoff=inps.cutoff,
        )
    else:
        box_list = split_into_sub_boxes(
            ds_shape=(ts_obj.length, ts_obj.width),
            step=step,
        )

    region_docs = create_kml_region_document(inps, box_list, ts_obj, step)

    write_network_link_file(region_docs, ts_obj, box_list, lod, net_link_file)

    net_link = create_network_link_element(os.path.relpath(net_link_file, start=inps.work_dir), ts_obj)

    return net_link


######################################################################################
def save_kmz_timeseries(inps):

    # resource file paths
    inps.work_dir = os.path.abspath(os.path.dirname(inps.ts_file))
    inps.cbar_file = os.path.join(inps.work_dir, "google_earth_cbar.png")
    inps.star_file = os.path.join(inps.work_dir, "star.png")
    inps.dot_file  = os.path.join(inps.work_dir, "shaded_dot.png")
    inps.dygraph_file = os.path.join(inps.work_dir, "dygraph-combined.js")
    inps.kml_data_dir = os.path.join(inps.work_dir, "kml_data")

    if inps.outfile:
        fbase = os.path.splitext(os.path.basename(inps.outfile))[0]
    else:
        fbase = pp.auto_figure_title(inps.ts_file, inps_dict=vars(inps))
    root_file = os.path.join(inps.work_dir, f"{fbase}_root.kml")
    kmz_file = os.path.join(inps.work_dir, f"{fbase}.kmz")

    ## read data
    ts_obj = timeseries(inps.ts_file)
    ts_obj.open()
    length, width = ts_obj.length, ts_obj.width
    inps.metadata = ts_obj.metadata
    lats, lons = ut.get_lat_lon(ts_obj.metadata)
    print(f'input data shape in row/col: {length}/{width}')

    vel = readfile.read(inps.vel_file, datasetName='velocity')[0] * 100.

    # Set vmin/max and colormap
    inps.vlim = inps.vlim if inps.vlim is not None else [np.nanmin(vel), np.nanmax(vel)]
    if inps.wrap:
        print(f're-wrapping data to {inps.vlim} cm/year for color coding')
    inps.colormap = pp.ColormapExt(inps.cmap_name).colormap


    ##--------- Create root KML file with network links to data KML files --------------##
    kml_root_doc = KML.Document()

    # 1 Create Overlay element for colorbar
    cbar_overlay = save_kmz.generate_cbar_element(
        cbar_file=inps.cbar_file,
        cmap=inps.colormap,
        vmin=inps.vlim[0],
        vmax=inps.vlim[1],
    )
    kml_root_doc.append(cbar_overlay)

    # 2 Generate the placemark for the Reference Pixel
    ref_point = create_reference_point_element(inps, lats, lons, ts_obj)
    print('add reference point.')
    ref_folder = KML.Folder(KML.name("ReferencePoint"))
    ref_folder.append(ref_point)
    kml_root_doc.append(ref_folder)

    # 3 Create data folder to contain actual data elements
    data_folder = KML.Folder(KML.name("Data"))
    for i, step in enumerate(inps.steps):
        net_link = generate_network_link(
            inps,
            ts_obj,
            step=step,
            lod=(inps.lods[i], inps.lods[i+1]),
        )
        if net_link is not None:
            data_folder.append(net_link)
    kml_root_doc.append(data_folder)


    ##---------------------------- Write root KML file ------------------------------##
    print('-'*30)
    print(f'writing {root_file}')
    kml_root = KML.kml()
    kml_root.append(kml_root_doc)
    with open(root_file, 'w') as f:
        f.write(etree.tostring(kml_root, pretty_print=True).decode('utf-8'))

    ## Copy auxiliary files
    res_dir = os.path.join(os.path.dirname(mintpy.__file__), "data")
    for fname in [inps.star_file, inps.dot_file, inps.dygraph_file]:
        src_file = os.path.join(res_dir, os.path.basename(fname))
        shutil.copy2(src_file, inps.work_dir)
        print(f"copy {src_file} to the local directory")

    ## Generate KMZ file
    # 1) go to the directory of kmz file
    run_dir = os.path.abspath(os.getcwd())
    os.chdir(inps.work_dir)

    # 2) zip all data files
    with ZipFile(kmz_file, 'w') as fz:
        kml_data_files = get_all_file_paths(inps.kml_data_dir)
        for fname in [root_file,
                      inps.cbar_file,
                      inps.dygraph_file,
                      inps.dot_file,
                      inps.star_file] + kml_data_files:
            fz.write(os.path.relpath(fname))
            if not inps.keep_kml_file:
                os.remove(fname)
        if not inps.keep_kml_file:
            shutil.rmtree(inps.kml_data_dir)

    # 3) go back to the running directory
    os.chdir(run_dir)
    print(f'merged all files to {kmz_file}')
    print('Done.')
    print(f'Open {kmz_file} in Google Earth and play!')

    return
