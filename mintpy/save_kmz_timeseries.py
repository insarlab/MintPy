#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Joshua Zahner, Zhang Yunjun, 2018                #
############################################################


import os
import sys
import argparse
from lxml import etree
from zipfile import ZipFile
import shutil
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt, patches

try:
    from pykml.factory import KML_ElementMaker as KML
except ImportError:
    raise ImportError('Can not import pykml!')

import mintpy
from mintpy.objects import timeseries, deramp
from mintpy.utils import readfile, plot, utils as ut
from mintpy import save_kmz


############################################################
EXAMPLE = """example:
  cd $PROJECT_NAME/mintpy/geo
  save_kmz_timeseries.py geo_timeseries_ERA5_ramp_demErr.h5
  save_kmz_timeseries.py geo_timeseries_ERA5_ramp_demErr.h5 -v -5 5 --wrap
  save_kmz_timeseries.py timeseries_ERA5_demErr.h5 --vel velocity.h5 --tcoh temporalCoherence.h5 --mask maskTempCoh.h5
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Generare Google Earth KMZ file for time-series file.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    args = parser.add_argument_group('Input files', 'File/Dataset to display')

    args.add_argument('ts_file', metavar='timeseries_file', help='Timeseries file to generate KML for')
    args.add_argument('--vel', dest='vel_file', metavar='FILE',
                      help='Velocity file, used for the color of dot')
    args.add_argument('--tcoh', dest='tcoh_file', metavar='FILE',
                      help='temporal coherence file, used for stat info')
    args.add_argument('--mask', dest='mask_file', metavar='FILE',
                      help='Mask file')
    args.add_argument('-o','--output', dest='outfile', help='Output KMZ file name.')

    opts = parser.add_argument_group('Display options', 'configurations for the display')
    opts.add_argument('--steps', type=int, nargs=3, default=[20, 5, 2],
                      help='list of steps for output pixel. Default: 20 5 2')
    opts.add_argument('--level-of-details','--lods', dest='lods', type=int, nargs=4, default=[0, 1500, 4000, -1],
                      help='list of level of details to determine the visible range while browering. Default: 0, 1500, 4000, -1.\n'+
                           'Ref: https://developers.google.com/kml/documentation/kml_21tutorial')
    opts.add_argument('--vlim','-v', dest='vlim', nargs=2, metavar=('VMIN', 'VMAX'), type=float,
                      help='min/max range in cm/yr for color coding.')
    opts.add_argument('--wrap', dest='wrap', action='store_true',
                      help='re-wrap data to [VMIN, VMAX) for color coding.')
    opts.add_argument('--colormap','-c', dest='colormap', default='jet',
                      help='colormap used for display, i.e. jet, RdBu, hsv, jet_r, temperature, viridis,  etc.\n'
                           'colormaps in Matplotlib - http://matplotlib.org/users/colormaps.html\n'
                           'colormaps in GMT - http://soliton.vm.bytemark.co.uk/pub/cpt-city/')

    defo = parser.add_argument_group('HD for deforming areas', 'High resolution output for deforming areas')
    defo.add_argument('--cutoff', dest='cutoff', type=int, default=3,
                      help='choose points with velocity >= cutoff * MAD. Default: 3.')
    defo.add_argument('--min-percentage','--min-perc', dest='min_percentage', type=float, default=0.2,
                      help='choose boxes with >= min percentage of pixels are deforming. Default: 0.2.')

    parser.add_argument('--kk','--keep-kml','--keep-kml-file', dest='keep_kml_file', action='store_true',
                        help='Do not remove KML and data/resource files after compressing into KMZ file.')

    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check if in geo coordinates
    atr = readfile.read_attribute(inps.ts_file)
    if "Y_FIRST" not in atr.keys():
        raise ValueError("input file {} is NOT geocoded".format(inps.ts_file))

    inps = get_aux_filename(inps)
    for fname in [inps.vel_file, inps.tcoh_file, inps.mask_file]:
        if not os.path.isfile(fname):
            raise FileNotFoundError('auxliary file {} not found.'.format(fname))
    return inps


def get_all_file_paths(directory):
    # initializing empty file paths list
    file_paths = []

    # crawling through directory and subdirectories
    for root, directories, files in os.walk(directory):
        for filename in files:
            # join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)

    # returning all file paths
    return file_paths


def get_aux_filename(inps):
    """Get auxliary files' default filename"""
    # path feature of time-series file
    ts_dir = os.path.dirname(inps.ts_file)
    ts_prefix = os.path.basename(inps.ts_file).split('timeseries')[0]

    if not inps.vel_file:
        inps.vel_file = os.path.join(ts_dir, '{}velocity.h5'.format(ts_prefix))
    if not inps.tcoh_file:
        inps.tcoh_file = os.path.join(ts_dir, '{}temporalCoherence.h5'.format(ts_prefix))
    if not inps.mask_file:
        inps.mask_file = os.path.join(ts_dir, '{}maskTempCoh.h5'.format(ts_prefix))
    return inps


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
    print('step: {} pixels'.format(step))
    win_size = int(np.sqrt(num_pixel) * step)
    length, width = ds_shape

    # number of sub boxes
    nrows = np.ceil(length / win_size).astype(int)
    ncols = np.ceil(width / win_size).astype(int)
    if print_msg:
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


def get_boxes4deforming_area(vel_file, mask_file, step=2, num_pixel=30**2, min_percentage=0.2,
                             cutoff=3, ramp_type='quadratic', display=False):
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
    win_size = int(np.sqrt(num_pixel) * step)
    print('-'*30)
    print('get boxes on deforming areas with step: {} pixels'.format(step))
    mask = readfile.read(mask_file)[0]
    vel, atr = readfile.read(vel_file)
    print('removing a {} phase ramp from input velocity before the evaluation'.format(ramp_type))
    vel = deramp(vel, mask, ramp_type=ramp_type, metadata=atr)[0]               #remove ramp before the evaluation

    # get deforming pixels
    mad = ut.median_abs_deviation_threshold(vel[mask], cutoff=cutoff)     #deformation threshold
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
        print('save figure to {}'.format(out_fig))
        plt.show()
    return box_list


def create_reference_point_element(inps, lats, lons, ts_obj):
    """Create reference point element"""
    colormap = mpl.cm.get_cmap(inps.colormap)  # set colormap
    norm = mpl.colors.Normalize(vmin=inps.vlim[0], vmax=inps.vlim[1])

    ref_yx = (int(ts_obj.metadata['REF_Y']), int(ts_obj.metadata['REF_X']))
    ref_lalo = (lats[ref_yx[0], ref_yx[1]], lons[ref_yx[0], ref_yx[1]])

    ref_point = KML.Placemark(
        KML.Style(
            KML.IconStyle(
                KML.color(save_kmz.get_hex_color(0.0, colormap, norm)),
                KML.scale(1.),
                KML.Icon(
                    KML.href("{}".format(os.path.basename(inps.star_file)))
                )
            )
        ),
        KML.description("Reference point <br /> \n <br /> \n" +
                        get_description_string(ref_lalo, ref_yx, 0.00, 0.00, 0.00, 1.00)
                        ),
        KML.Point(
            KML.coordinates("{}, {}".format(ref_lalo[1], ref_lalo[0]))
        )
    )

    return ref_point


def get_description_string(coords, yx, v, vstd, disp, tcoh=None, font_size=4):
    """Description information of each data point."""
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
    des_str += "*Double click to reset plot <br /> <br />\n"
    des_str += "\n\n"
    return des_str


def generate_js_datastring(dates, dygraph_file, num_date, ts):
    """String of the Java Script for interactive plot of diplacement time-series"""
    dygraph_file = '../../../{}'.format(os.path.basename(dygraph_file))
    js_data_string = "<script type='text/javascript' src='{}'></script>".format(dygraph_file)
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
        date_displacement_string = "\"{}, {}\\n\" + \n".format(date, dis)
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
    dot_file = '../../{}'.format(os.path.basename(inps.dot_file))

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

        # 2.1 Set and normalize colormap to defined vlim
        colormap = mpl.cm.get_cmap(inps.colormap)
        norm = mpl.colors.Normalize(vmin=inps.vlim[0], vmax=inps.vlim[1])

        # 2.2 Set number of pixels to use
        num_pixel = int(length / step) * int(width / step)
        msg = "create KML doc for box {}/{}: {}".format(i+1, num_box, box)
        msg += ", step: {} pixels, {} pixels in total ...".format(step, num_pixel)
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
                            KML.color(save_kmz.get_hex_color(vc, colormap, norm)),
                            KML.scale(0.5),
                            KML.Icon(KML.href("{}".format(dot_file)))
                        )
                    )

                    # 2.3.2 Create KML point element
                    point = KML.Point(KML.coordinates("{},{}".format(lon, lat)))

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
    print("create KML region links directory: {}".format(os.path.basename(links_dir)))

    ## 2. Create root KML element and KML Document element
    kml = KML.kml()
    kml_doc = KML.Document()

    ## 3. Generate a new network link element for each region
    for num, (region_doc, box) in enumerate(zip(region_docs, box_list)):
        region_kml_file = os.path.join(links_dir, "region_{}.kml".format(num))

        ## 3.1 Write the first region_document to a file and move it to the proper subdircetory
        kml_1 = KML.kml()
        kml_1.append(region_doc)
        with open(region_kml_file, 'w') as f:
            f.write(etree.tostring(kml_1, pretty_print=True).decode('utf-8'))

        ## 3.2 Flatten lats and lons data
        lats, lons = flatten_lat_lon(box, ts_obj)

        ## 3.3 Define new NetworkLink element
        network_link = KML.NetworkLink(
            KML.name('Region {}'.format(num)),
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
        print('skip step = {}'.format(step))
        return None

    elif step == inps.steps[-1]:
        box_list = get_boxes4deforming_area(inps.vel_file, inps.mask_file,
                                            step=inps.steps[-1],
                                            min_percentage=inps.min_percentage,
                                            cutoff=inps.cutoff)
    else:
        box_list = split_into_sub_boxes((ts_obj.length, ts_obj.width), step=step)

    region_docs = create_kml_region_document(inps, box_list, ts_obj, step)

    write_network_link_file(region_docs, ts_obj, box_list, lod, net_link_file)

    net_link = create_network_link_element(os.path.relpath(net_link_file, start=inps.work_dir), ts_obj)

    return net_link


def main(iargs=None):
    inps = cmd_line_parse(iargs)
    inps.work_dir = os.path.abspath(os.path.dirname(inps.ts_file))
    inps.cbar_file = os.path.join(inps.work_dir, 'google_earth_cbar.png')
    inps.star_file = os.path.join(inps.work_dir, "star.png")
    inps.dot_file = os.path.join(inps.work_dir, "shaded_dot.png")
    inps.dygraph_file = os.path.join(inps.work_dir, "dygraph-combined.js")
    inps.kml_data_dir = os.path.join(inps.work_dir, 'kml_data')

    ## Define file names
    if inps.outfile:
        inps.outfile_base = os.path.splitext(os.path.basename(inps.outfile))[0]
    else:
        inps.outfile_base = plot.auto_figure_title(inps.ts_file, inps_dict=vars(inps))
    kml_root_file = os.path.join(inps.work_dir, '{}_root.kml'.format(inps.outfile_base))
    kmz_file = os.path.join(inps.work_dir, '{}.kmz'.format(inps.outfile_base))

    ## read data
    ts_obj = timeseries(inps.ts_file)
    ts_obj.open()
    length, width = ts_obj.length, ts_obj.width
    inps.metadata = ts_obj.metadata
    lats, lons = ut.get_lat_lon(ts_obj.metadata)
    print('input data shape in row/col: {}/{}'.format(length, width))

    vel = readfile.read(inps.vel_file, datasetName='velocity')[0] * 100.
    # Set min/max velocity for colormap
    if inps.vlim is None:
        inps.vlim = [np.nanmin(vel), np.nanmax(vel)]
    if inps.wrap:
        print('re-wrapping data to {} cm/year for color coding'.format(inps.vlim))

    ##--------- Create root KML file with network links to data KML files --------------##
    kml_root_doc = KML.Document()

    # 1 Create Overlay element for colorbar
    cbar_overlay = save_kmz.generate_cbar_element(
        cbar_file=inps.cbar_file,
        vmin=inps.vlim[0],
        vmax=inps.vlim[1],
        cmap=inps.colormap,
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
        net_link = generate_network_link(inps, ts_obj,
                                         step=step,
                                         lod=(inps.lods[i], inps.lods[i+1]))
        if net_link is not None:
            data_folder.append(net_link)
    kml_root_doc.append(data_folder)


    ##---------------------------- Write root KML file ------------------------------##
    print('-'*30)
    print('writing ' + kml_root_file)
    kml_root = KML.kml()
    kml_root.append(kml_root_doc)
    with open(kml_root_file, 'w') as f:
        f.write(etree.tostring(kml_root, pretty_print=True).decode('utf-8'))

    ## Copy auxiliary files
    res_dir = os.path.join(os.path.dirname(mintpy.__file__), "data")
    for fname in [inps.star_file, inps.dot_file, inps.dygraph_file]:
        src_file = os.path.join(res_dir, os.path.basename(fname))
        shutil.copy2(src_file, inps.work_dir)
        print("copy {} to the local directory".format(src_file))

    ## Generate KMZ file
    # 1) go to the directory of kmz file
    run_dir = os.path.abspath(os.getcwd())
    os.chdir(inps.work_dir)

    # 2) zip all data files
    with ZipFile(kmz_file, 'w') as fz:
        kml_data_files = get_all_file_paths(inps.kml_data_dir)
        for fname in [kml_root_file, 
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
    print('merged all files to {}'.format(kmz_file))
    print('Done.')
    print('Open {} in Google Earth and play!'.format(kmz_file))
    return


######################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
