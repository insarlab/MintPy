#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2018-2019, Joshua Zahner, Zhang Yunjun     #
# Author:  Joshua Zahner, Zhang Yunjun                    #
############################################################

import os
import sys
import argparse

try:
    from pykml.factory import KML_ElementMaker as KML
except ImportError:
    raise ImportError('Can not import pykml!')

from lxml import etree
import numpy as np
import matplotlib as mpl
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt

from pysar.objects import timeseries, deramp
from pysar.utils import readfile, plot, utils as ut


############################################################
EXAMPLE = """example:
  cd $PROJECT_NAME/PYSAR/GEOCODE
  save_kmz_timeseries.py geo_timeseries_ECMWF_ramp_demErr.h5
  save_kmz_timeseries.py geo_timeseries_ECMWF_ramp_demErr.h5 --vel geo_velocity_masked.h5 --tcoh geo_temporalCoherence.h5
"""

def create_parser():

    parser = argparse.ArgumentParser(description='Generare Google Earth KMZ file for time-series file.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    args = parser.add_argument_group('Input files', 'File/Dataset to display')

    args.add_argument('ts_file', metavar='timeseries_file', help='Timeseries file to generate KML for')
    args.add_argument('--vel', dest='vel_file', metavar='velocity_file', default='geo_velocity_masked.h5',
                      help='Velocity file')
    args.add_argument('--tcoh', dest='tcoh_file', metavar='temporal_coh_file', default='geo_temporalCoherence.h5',
                      help='temporal coherence file')
    args.add_argument('--mask', dest='mask_file', metavar='mask_file', default='geo_maskTempCoh.h5', help='Mask file')

    opts = parser.add_argument_group('Display options', 'configurations for the display')
    opts.add_argument('--step', type=int, metavar='NUM', default=3, help='output pixel step number, default: 3')
    opts.add_argument('-v','--vlim', dest='vlim', nargs=2, metavar=('VMIN', 'VMAX'), type=float,
                      help='Display limits for matrix plotting.')
    opts.add_argument('-c', '--colormap', dest='colormap', default='jet',
                      help='colormap used for display, i.e. jet, RdBu, hsv, jet_r, temperature, viridis,  etc.\n'
                           'colormaps in Matplotlib - http://matplotlib.org/users/colormaps.html\n'
                           'colormaps in GMT - http://soliton.vm.bytemark.co.uk/pub/cpt-city/')
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


def get_lat_lon(meta, box=None):
    """extract lat/lon info of all grids into 2D matrix
    Parameters: meta : dict, including X/Y_FIRST/STEP and LENGTH/WIDTH info
                box  : 4-tuple of int for (x0, y0, x1, y1)
    Returns:    lats : 2D np.array for latitude  in size of (length, width)
                lons : 2D np.array for longitude in size of (length, width)
    """
    length, width = int(meta['LENGTH']), int(meta['WIDTH'])
    if box is None:
        box = (0, 0, width, length)

    # generate 2D matrix for lat/lon
    lat_num = box[3] - box[1]
    lon_num = box[2] - box[0]
    lat_step = float(meta['Y_STEP'])
    lon_step = float(meta['X_STEP'])
    lat0 = float(meta['Y_FIRST']) + lat_step * box[1]
    lon0 = float(meta['X_FIRST']) + lon_step * box[0]
    lat1 = lat0 + lat_step * (lat_num - 1)
    lon1 = lon0 + lon_step * (lon_num - 1)
    lats, lons = np.mgrid[lat0:lat1:lat_num*1j,
                          lon0:lon1:lon_num*1j]
    return lats, lons

def flatten_lalos(box, ts_obj, coords=None):
    if coords is None:
        lats, lons = get_lat_lon(ts_obj.metadata, box=box)

    lats = sorted(lats.flatten())
    lons = sorted(lons.flatten())
    return lats, lons

def split_into_sub_boxes(ds_shape, win_size=300, print_msg=True):
    """split input shape into multiple sub boxes
    Parameters: ds_shape : 2-tuple of int for the shape of whole dataset in (length, width)
                win_size : int, size of path in one direction
    Returns:    box_list : list of 4-tuple of int, each indicating (col0, row0, col1, row1)
    """
    length, width = ds_shape
    # number of sub boxes
    nrows = np.ceil(length / win_size).astype(int)
    ncols = np.ceil(width / win_size).astype(int)
    if print_msg:
        print('input data shape in row/col: {}/{}'.format(length, width))
        print('max box size: {}'.format(win_size))
        print('number of output boxes: {}'.format(nrows * ncols))
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

def get_boxes4deforming_area(vel_file, mask_file, win_size=30, min_percentage=0.2, ramp_type='quadratic', display=False):
    """Get list of boxes to cover the deforming areas.
    A pixel is identified as deforming if its velocity exceeds the MAD of the whole image.
    Parameters: vel_file : str, path of velocity file
                mask_file : str, path of mask file
                win_size  : int, length and width of the output box
                min_percentage : float between 0 and 1, minimum percentage of deforming points in the box
                ramp_type : str, type of phase ramps to be removed while evaluating the deformation
                display   : bool, plot the identification result or not
    Returns:    box_list  : list of t-tuple of int, each indicating (col0, row0, col1, row1)
    """
    print('-'*30)
    print('get boxes on deforming areas')
    mask = readfile.read(mask_file)[0]
    vel, atr = readfile.read(vel_file)
    print('removing a {} phase ramp from input velocity before the evaluation'.format(ramp_type))
    vel = deramp(vel, mask, ramp_type=ramp_type, metadata=atr)[0]               #remove ramp before the evaluation

    # get deforming pixels
    mad = ut.median_abs_deviation_threshold(vel[mask], center=0., cutoff=3)     #deformation threshold
    print('velocity threshold / median abs dev: {:.3f} cm/yr'.format(mad))
    vel[mask == 0] = 0
    mask_aoi = (vel >= mad) + (vel <= -1. * mad)
    print('number of points: {}'.format(np.sum(mask_aoi)))

    # get deforming boxes
    box_list = []
    min_num = min_percentage * (win_size ** 2)
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
    print('number of boxes : {}'.format(len(box_list)))

    if display:
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=[12, 8], sharey=True)
        vel[mask == 0] = np.nan
        axs[0].imshow(vel, cmap='jet')
        axs[1].imshow(mask_aoi, cmap='gray')
        for box in box_list:
            for ax in axs:
                rect = Rectangle((box[0],box[1]), (box[2]-box[0]), (box[3]-box[1]), linewidth=2, edgecolor='r', fill=False)
                ax.add_patch(rect)
        plt.show()
    return box_list


def plot_colorbar(out_file, vmin, vmax, cmap='jet', figsize=(0.18, 3.6)):
    fig, cax = plt.subplots(figsize=figsize)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)  # normalize velocity colors between 0.0 and 1.0
    cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical')
    cbar.set_label('{} [{}]'.format("Mean LOS velocity", "cm/year"), fontsize=12)
    cbar.locator = mpl.ticker.MaxNLocator(nbins=7)
    cbar.update_ticks()
    cbar.ax.tick_params(which='both', labelsize=12)
    fig.patch.set_facecolor('white')
    fig.patch.set_alpha(0.7)
    print('writing', out_file)
    fig.savefig(out_file, bbox_inches='tight', facecolor=fig.get_facecolor(), dpi=300)
    return out_file

def get_color_for_velocity(v, colormap, norm):
    rgba = colormap(norm(v))  # get rgba color components for point velocity
    hex = mpl.colors.to_hex([rgba[3], rgba[2], rgba[1], rgba[0]], keep_alpha=True)[1:]
    return hex


def generate_description_string(coords, yx, v, vstd, disp, tcoh=None, font_size=4):
    des_str = "<font size={}>".format(font_size)
    des_str += "Latitude: {:.6f}˚ <br /> \n".format(coords[0])
    des_str += "Longitude: {:.6f}˚ <br /> \n".format(coords[1])
    des_str += "Row: {:.0f} <br /> \n".format(yx[0])
    des_str += "Column: {:.0f} <br /> \n".format(yx[1])
    des_str += " <br /> \n"
    des_str += "Mean LOS velocity [cm/year]: {:.2f} +/- {:.2f} <br /> \n".format(v, vstd)
    des_str += "Cumulative displacement [cm]: {:.2f} <br /> \n".format(disp)
    if tcoh is not None:
        des_str += "Temporal coherence: {:.2f} <br /> \n".format(tcoh)
    des_str += "</font>"
    des_str += " <br />  <br /> "
    des_str += "* Double click to reset plot <br /> <br />\n"
    des_str += "\n\n"
    return des_str

def generate_js_datastring(dates, dygraph_file, num_date, ts):

    js_data_string = "<script type='text/javascript' src='../../../" + dygraph_file + "'></script>"
    js_data_string += """
        <div id='graphdiv'> </div>
        <style>
            .dygraph-legend{
                left: 435px !important;
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
        date_displacement_string = "\"{}, {}\\n\" + \n".format(date, dis)
        js_data_string += date_displacement_string

    js_data_string += """
    
    "",
       {
         width: 700,
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
                     var dateString = 'Date: ' + ('0' + date.getDate()).slice(-2) + '/' + ('0' + (date.getMonth() + 1)).slice(-2) + '/' + date.getFullYear()
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


def create_kml_document(inps, box_list, ts_obj, step):

    dot_file = "shaded_dot.png"
    dygraph_file = "dygraph-combined.js"

    ## 1. read data file into timeseries object

    kml_region_documents = []

    for i in range(len(box_list)):
        print(i)
        box = box_list[i]

        if box is None:
            box = (0, 0, width, length)

        length = box[3] - box[1]
        width = box[2] - box[0]

        # 1.1 Parse Date
        dates = np.array(ts_obj.times)  # 1D np.array of dates in datetime.datetime object in size of [num_date,]
        dates = list(map(lambda d: d.strftime("%Y-%m-%d"), dates))
        num_date = len(dates)

        # 1.2 Parse Spatial coordinates
        lats, lons = get_lat_lon(ts_obj.metadata, box=box)
        rows, cols = np.mgrid[box[1]:box[3] - 1:length * 1j, box[0]:box[2] - 1:width * 1j]

        # 1.3 Read Velocity / time-series / temporal coherence data
        print('read velocity data')
        vel = readfile.read(inps.vel_file, datasetName='velocity', box=box)[0] * 100.
        vel_std = readfile.read(inps.vel_file, datasetName='velocityStd', box=box)[0] * 100.
        print('read time-series data')
        ts_data = readfile.read(inps.ts_file, box=box)[0] * 100.
        ts_data -= np.tile(ts_data[0, :, :], (ts_data.shape[0], 1, 1))  # enforce displacement starts from zero
        print('read temporal coherence data')
        temp_coh = readfile.read(inps.tcoh_file, box=box)[0]
        mask = ~np.isnan(vel)


        ## 2. Create KML Document
        print('create KML file.')
        kml_document = KML.Document()

        # 2.1 Set and normalize colormap to defined vlim
        colormap = mpl.cm.get_cmap(inps.colormap)
        norm = mpl.colors.Normalize(vmin=inps.vlim[0], vmax=inps.vlim[1])

        # 2.2 Set number of pixels to use
        num_pixel = int(length / step) * int(width / step)
        msg = "add point time-series data "
        msg += "(select 1 out of every {} by {} pixels; select {} pixels in total) ...".format(step, step, num_pixel)
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
                    vstd = vel_std[i, j]
                    tcoh = temp_coh[i, j]

                    # 2.3.1 Create KML icon style element
                    style = KML.Style(
                        KML.IconStyle(
                            KML.color(get_color_for_velocity(v, colormap, norm)),
                            KML.scale(0.5),
                            KML.Icon(KML.href("../../{}".format(dot_file)))
                        )
                    )

                    # 2.3.2 Create KML point element
                    point = KML.Point(KML.coordinates("{},{}".format(lon, lat)))

                    js_data_string = generate_js_datastring(dates, dygraph_file, num_date, ts)

                    # 2.3.3 Create KML description element
                    stats_info = generate_description_string((lat, lon), (row, col), v, vstd, ts[-1], tcoh=tcoh)
                    description = KML.description(stats_info, js_data_string)

                    # 2.3.4 Crate KML Placemark element to hold style, description, and point elements
                    placemark = KML.Placemark(style, description, point)

                    # 2.3.5 Append each placemark element to the KML document object
                    data_folder.append(placemark)

        # 2.4 Append data folder to KML document
        kml_document.append(data_folder)

        # 2.5 Append KML document to list of regionalized documents
        kml_region_documents.append(kml_document)

    return kml_region_documents, dot_file, dygraph_file

def create_regionalized_networklinks_file(regions, ts_obj, box_list, lod, step, output_file):

    ## 1. Create directory to store regioalized KML data files
    kml_data_files_directory = "kml_data_files"
    links_directory = "{0}by{0}_links".format(step)
    mkdir_kml_data_file_command = "mkdir {}".format(kml_data_files_directory)
    os.system(mkdir_kml_data_file_command)
    cmdDirectory = "cd {}; mkdir {}/".format(kml_data_files_directory, links_directory)
    print("Creating KML links directory")
    os.system(cmdDirectory)

    ## 2. Create master KML element and KML Document element
    kml = KML.kml()
    kml_document = KML.Document()

    ## 3. Define min and max levels of detail for large vs small data
    min_lod = lod[0]
    max_lod = lod[1]

    ## 4. Define list of regionalized boxes and list of number of regions
    region_nums = list(range(0, len(regions)))

    ## 5. Generate a new network link element for each region
    for ele in zip(regions, box_list, region_nums):

        region_document = ele[0]
        box = ele[1]
        num = ele[2]

        kml_file = "region_{}.kml".format(num)

        ## 5.1 Write the first region_document to a file and move it to the proper subdircetory
        kml_1 = KML.kml()
        kml_1.append(region_document)
        print('writing ' + kml_file)
        with open(kml_file, 'w') as f:
            f.write(etree.tostring(kml_1, pretty_print=True).decode('utf-8'))

        cmdMove = "mv {sm} {dir}/{sm}".format(sm=kml_file, dir="{}/{}".format(kml_data_files_directory, links_directory))
        print("Moving KML Data Files to directory")
        os.system(cmdMove)

        ## 5.2 Flatten lats and lons data
        lats, lons = flatten_lalos(box, ts_obj)

        ## 5.3 Define new NetworkLink element
        network_link = KML.NetworkLink(
            KML.name('Region {}'.format(num)),
            KML.visibility(1),
            KML.Region(
                KML.Lod(
                    KML.minLodPixels(min_lod),
                    KML.maxLodPixels(max_lod)
                ),
                KML.LatLonAltBox(
                    KML.north(lats[0] + 0.5),
                    KML.south(lats[-1] - 0.5),
                    KML.east(lons[0] + 0.5),
                    KML.west(lons[-1] - 0.5)
                )
            ),
            KML.Link(
                KML.href("{}/{}".format(links_directory, kml_file)),
                KML.viewRefreshMode('onRegion')
            )
        )

        ## 5.4 Append new NetworkLink to KML document
        kml_document.append(network_link)

    ## 6. Write the full KML document to the output file and move it to the proper directory
    kml.append(kml_document)
    with open(output_file, 'w') as f:
        f.write(etree.tostring(kml, pretty_print=True).decode('utf-8'))

    return "{}/{}".format(kml_data_files_directory, output_file), kml_data_files_directory

def create_network_link_element(name, kml_file, ts_obj):

    lats, lons = flatten_lalos(None, ts_obj)

    network_link = KML.NetworkLink(
        KML.name(name),
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
            KML.href(kml_file),
            KML.viewRefreshMode('onRegion')
        )
    )
    return network_link

def create_reference_point_element(inps, lats, lons, ts_obj):

    star_file = "star.png"

    colormap = mpl.cm.get_cmap(inps.colormap)  # set colormap
    norm = mpl.colors.Normalize(vmin=inps.vlim[0], vmax=inps.vlim[1])

    ref_yx = (int(ts_obj.metadata['REF_Y']), int(ts_obj.metadata['REF_X']))
    ref_lalo = (lats[ref_yx[0], ref_yx[1]], lons[ref_yx[0], ref_yx[1]])

    reference_point = KML.Placemark(
        KML.Style(
            KML.IconStyle(
                KML.color(get_color_for_velocity(0.0, colormap, norm)),
                KML.scale(1.),
                KML.Icon(KML.href("{}".format(star_file)))
            )
        ),
        KML.description("Reference point <br /> \n <br /> \n" +
                        generate_description_string(ref_lalo, ref_yx, 0.00, 0.00, 0.00, 1.00)
                        ),
        KML.Point(
            KML.coordinates("{}, {}".format(ref_lalo[1], ref_lalo[0]))
        )
    )

    return reference_point, star_file

def generate_cbar_element(inps):

    out_name_base = plot.auto_figure_title(inps.ts_file, inps_dict=vars(inps))
    cbar_png_file = '{}_cbar.png'.format(out_name_base)

    cbar_png_file = plot_colorbar(out_file=cbar_png_file, vmin=inps.vlim[0], vmax=inps.vlim[1], cmap=inps.colormap)
    cbar_overlay = KML.ScreenOverlay(
        KML.name('Colorbar'),
        KML.Icon(
            KML.href("{}".format(cbar_png_file)),
            KML.viewBoundScale(0.75)
        ),
        KML.overlayXY(x="0", y="0", xunits="fraction", yunits="fraction"),
        KML.screenXY(x="0", y="0", xunits="fraction", yunits="fraction"),
        KML.size(x="0", y="250", xunits="pixel", yunits="pixel"),
        KML.rotation(0),
        KML.visibility(1),
        KML.open(0)
    )
    print('add colorbar.')
    return cbar_overlay, cbar_png_file


def generate_network_link(inps, ts_obj, box_list, step, lod, output_file=None):

    out_name_base = plot.auto_figure_title(inps.ts_file, inps_dict=vars(inps))

    if output_file is None:
        output_file = "{0}_{1}by{1}.kml".format(out_name_base, step)

    nwklink_name = "{0} by {0}".format(step)

    regions, dot_file, dygraph_file = create_kml_document(inps, box_list, ts_obj, step)
    kml_file, kml_data_files_directory = create_regionalized_networklinks_file(regions, ts_obj, box_list, lod, step, output_file)

    network_link = create_network_link_element(nwklink_name, kml_file, ts_obj)

    cmdMove = "mv {0} {1}/{0}".format(output_file, kml_data_files_directory)
    print("Moving KML Data Files to directory")
    os.system(cmdMove)

    return network_link, dot_file, dygraph_file, kml_data_files_directory


def main(iargs=None):

    ## 1. Read command line variables
    inps = cmd_line_parse(iargs)

    ## 2. Define file names
    out_name_base = plot.auto_figure_title(inps.ts_file, inps_dict=vars(inps))
    kml_file_master = '{}_master.kml'.format(out_name_base)
    kmz_file = '{}.kmz'.format(out_name_base)

    ## 4. read data
    ts_obj = timeseries(inps.ts_file)
    ts_obj.open()
    length, width = ts_obj.length, ts_obj.width
    lats, lons = get_lat_lon(ts_obj.metadata)

    vel = readfile.read(inps.vel_file, datasetName='velocity')[0] * 100.
    # Set min/max velocity for colormap
    if inps.vlim is None:
        inps.vlim = [np.nanmin(vel), np.nanmax(vel)]

    ## 5. Generate large and small KML files for different view heights
    small_dset_step = 20  # Increase this for coarser resolution in small dset
    large_dset_step = 3   # Decrease this for finer resolution in large dset
    full_dset_step = 1    # Full resolution for deforming regions

    box_list = split_into_sub_boxes((length, width))  # Create list of unique regions
    deforming_box_list = get_boxes4deforming_area(inps.vel_file, inps.mask_file)

    ## 6. Create master KML file with network links to data KML files
    kml_master = KML.kml()
    kml_master_document = KML.Document()

    # 6.1 Create Overlay element for colorbar
    cbar_overlay, cbar_png_file = generate_cbar_element(inps)
    kml_master_document.append(cbar_overlay)

    # 6.2 Generate the placemark for the Reference Pixel
    reference_point, star_file = create_reference_point_element(inps, lats, lons, ts_obj)
    print('add reference point.')
    reference_folder = KML.Folder(KML.name("ReferencePoint"))
    reference_folder.append(reference_point)
    kml_master_document.append(reference_folder)

    # 6.3 Create data folder to contain actual data elements
    data_folder = KML.Folder(KML.name("Data"))

    network_link_1, dot_file, _, _ = generate_network_link(inps, ts_obj, box_list, small_dset_step, (0, 1500))
    network_link_2, _, dygraph_file, _ = generate_network_link(inps, ts_obj, box_list, large_dset_step, (1500, 4000))
    network_link_3, _, _, kml_data_files_directory = generate_network_link(inps, ts_obj, deforming_box_list, full_dset_step, (4000, -1))

    # 6.3.3 Append network links to data folder
    data_folder.append(network_link_1)
    data_folder.append(network_link_2)
    data_folder.append(network_link_3)

    kml_master_document.append(data_folder)

    kml_master.append(kml_master_document)


    ## 7 Write master KML file
    print('writing ' + kml_file_master)
    with open(kml_file_master, 'w') as f:
        f.write(etree.tostring(kml_master, pretty_print=True).decode('utf-8'))


    ## 8 Copy auxiliary files

    # 8.1 shaded_dot file
    dot_path = os.path.join(os.path.dirname(__file__), "../docs/resources", dot_file)
    cmdDot = "cp {} {}".format(dot_path, dot_file)
    print("copying {} for point.".format(dot_file))
    os.system(cmdDot)

    # 8.2 star file
    star_path = os.path.join(os.path.dirname(__file__), "../docs/resources", star_file)
    cmdStar = "cp {} {}".format(star_path, star_file)
    print("copying {} for reference point.".format(star_file))
    os.system(cmdStar)

    # 8.3 dygraph-combined.js file
    dygraph_path = os.path.join(os.path.dirname(__file__), "../docs/resources", dygraph_file)
    cmdDygraph = "cp {} {}".format(dygraph_path, dygraph_file)
    print("copying {} for interactive plotting.".format(dygraph_file))
    os.system(cmdDygraph)

    ## 9. Generate KMZ file
    cmdKMZ = 'zip -r {} {} {} {} {} {} {}'.format(kmz_file, kml_data_files_directory, kml_file_master, cbar_png_file, dygraph_file, dot_file, star_file)
    print('writing {}\n{}'.format(kmz_file, cmdKMZ))
    os.system(cmdKMZ)

    # 9.1 Remove extra files from file tree after KMZ generation
    cmdClean = 'rm -r {} {} {} {} {} {}'.format(kml_data_files_directory, kml_file_master, cbar_png_file, dygraph_file, dot_file, star_file)
    os.system(cmdClean)

    print('Done.')
    return


######################################################################################
if __name__ == '__main__':
    main()
