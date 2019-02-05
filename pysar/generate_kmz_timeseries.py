#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
from lxml import etree
import matplotlib as mpl
import matplotlib.pyplot as plt
from pykml.factory import KML_ElementMaker as KML
from pysar.objects import timeseries
from pysar.utils import readfile, plot


def create_parser():

    parser = argparse.ArgumentParser(description='Generare Google Earth KML for timeseries HDF5 file.',
                                     formatter_class=argparse.RawTextHelpFormatter)

    args = parser.add_argument_group('Input File', 'File/Dataset to display')

    args.add_argument('ts_file', metavar='timeseries_file', help='Timeseries file to generate KML for')
    args.add_argument('--vel', dest='vel_file', metavar='velocity_file', default='geo_velocity_masked.h5',
                      help='Velocity file')
    args.add_argument('--tcoh', dest='tcoh_file', metavar='temporal_coh_file', default='geo_temporalCoherence.h5',
                      help='temporal coherence file')
    args.add_argument('--step', type=int, metavar='NUM', default=3, help='output pixel step number, default: 3')
    args.add_argument('-v','--vlim', dest='vlim', nargs=2, metavar=('VMIN', 'VMAX'), type=float,
                      help='Display limits for matrix plotting.')
    args.add_argument('-c', '--colormap', dest='colormap', default='jet',
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
    lat1 = lat0 + lat_step * lat_num
    lon1 = lon0 + lon_step * lon_num
    lats, lons = np.mgrid[lat0:lat1:lat_num*1j,
                          lon0:lon1:lon_num*1j]
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


def generate_description_string(coords, yx, v, vstd, disp, tcoh=None):
    des_str =  "Latitude: {:.6f}˚ <br /> \n".format(coords[0])
    des_str += "Longitude: {:.6f}˚ <br /> \n".format(coords[1])
    des_str += "Row: {:.0f} <br /> \n".format(yx[0])
    des_str += "Column: {:.0f} <br /> \n".format(yx[1])
    des_str += " <br /> \n"
    des_str += "Mean LOS velocity: {:.2f} +/- {:.2f} cm/year <br /> \n".format(v, vstd)
    des_str += "Cumulative displacement: {:.2f} cm <br /> \n".format(disp)
    if tcoh is not None:
        des_str += "Temporal coherence: {:.2f} <br /> \n".format(tcoh)
    des_str += " <br />  <br /> "
    des_str += "\n\n"
    return des_str


def plot_colorbar(out_file, vmin, vmax, cmap='jet', figsize=(5, 0.2)):
    fig, cax = plt.subplots(figsize=figsize)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)  # normalize velocity colors between 0.0 and 1.0
    cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')
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


def create_kml_document(inps, step, dot_file, dygraph_file):

    ## 1. read data
    ts_obj = timeseries(inps.ts_file)
    ts_obj.open()
    length, width = ts_obj.length, ts_obj.width

    box_list = split_into_sub_boxes((length, width))
    kml_region_documents = []

    for i in range(len(box_list)):
        print(i)
        box = box_list[i]

        if box is None:
            box = (0, 0, width, length)

        length = box[3] - box[1]
        width = box[2] - box[0]

        # 1.1 Date
        dates = np.array(ts_obj.times)  # 1D np.array of dates in datetime.datetime object in size of [num_date,]
        dates = list(map(lambda d: d.strftime("%Y-%m-%d"), dates))
        num_date = len(dates)

        # 1.2 Spatial coordinates
        lats, lons = get_lat_lon(ts_obj.metadata, box=box)
        rows, cols = np.mgrid[box[1]:box[3] - 1:length * 1j, box[0]:box[2] - 1:width * 1j]

        # 1.3 Velocity / time-series
        print('read velocity data')
        vel = readfile.read(inps.vel_file, datasetName='velocity', box=box)[0] * 100.
        vel_std = readfile.read(inps.vel_file, datasetName='velocityStd', box=box)[0] * 100.
        print('read time-series data')
        ts_data = readfile.read(inps.ts_file, box=box)[0] * 100.
        ts_data -= np.tile(ts_data[0, :, :], (ts_data.shape[0], 1, 1))  # enforce displacement starts from zero
        print('read temporal coherence data')
        temp_coh = readfile.read(inps.tcoh_file, box=box)[0]
        mask = ~np.isnan(vel)

        # data stats
        ts_min = np.nanmin(ts_data)
        ts_max = np.nanmax(ts_data)

        ## 2. Create KML Document
        print('create KML file.')
        kml_document = KML.Document()

        colormap = mpl.cm.get_cmap(inps.colormap)  # set colormap
        norm = mpl.colors.Normalize(vmin=inps.vlim[0], vmax=inps.vlim[1])

        # 2.3 Data folder for all points

        num_pixel = int(length / step) * int(width / step)
        msg = "add point time-series data "
        msg += "(select 1 out of every {} by {} pixels; select {} pixels in total) ...".format(step, step, num_pixel)
        print(msg)

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

                    # Create KML icon style element
                    style = KML.Style(
                        KML.IconStyle(
                            KML.color(get_color_for_velocity(v, colormap, norm)),
                            KML.scale(0.5),
                            KML.Icon(KML.href("../../{}".format(dot_file)))
                        )
                    )

                    # Create KML point element
                    point = KML.Point(KML.coordinates("{},{}".format(lon, lat)))

                    js_data_string = generate_js_datastring(dates, dygraph_file, num_date, ts, ts_max, ts_min)

                    # Create KML description element
                    stats_info = generate_description_string((lat, lon), (row, col), v, vstd, ts[-1], tcoh=tcoh)
                    description = KML.description(stats_info, js_data_string)

                    # Crate KML Placemark element to hold style, description, and point elements
                    placemark = KML.Placemark(style, description, point)
                    # Append each placemark element to the KML document object
                    data_folder.append(placemark)

        kml_document.append(data_folder)

        kml_region_documents.append(kml_document)

    return kml_region_documents


def generate_js_datastring(dates, dygraph_file, num_date, ts, ts_max, ts_min):
    # Javascript to embed inside the description
    js_data_string = "<script type='text/javascript' src='../../../" + dygraph_file + "'></script>\n" \
                      "<div id='graphdiv'> </div>\n" \
                      "<script type='text/javascript'>\n" \
                      "g = new Dygraph( document.getElementById('graphdiv'),\n" \
                      "\"Date, displacement\\n\" + \n"

    # append the date/displacement data
    for k in range(num_date):
        date = dates[k]
        dis = ts[k]
        date_displacement_string = "\"{}, {}\\n\" + \n".format(date, dis)
        js_data_string += date_displacement_string

    js_data_string += "\"\",\n" \
                      "{" \
                        "width: 700,\n" \
                        "height: 300,\n" \
                        "axes: {\n" \
                            "x: {\n" \
                                "axisLabelFormatter: function (d, gran) {\n" \
                                    "var months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']\n" \
                                    "var date = new Date(d)\n" \
                                    "var dateString = months[date.getMonth()] + ' ' + date.getFullYear()\n" \
                                    "return dateString;\n" \
                                "},\n" \
                                "valueFormatter: function (d) {\n" \
                                    "var date = new Date(d)\n" \
                                    "var dateString = 'Date: ' + ('0' + date.getDate()).slice(-2) + '/' + ('0' + (date.getMonth() + 1)).slice(-2) + '/' + date.getFullYear()\n" \
                                    "return dateString;\n" \
                                "},\n" \
                                "pixelsPerLabel: 90\n" \
                            "},\n" \
                            "y: {\n" \
                                "valueFormatter: function(v) {\n" \
                                    "return (' ' + v.toFixed(2)).slice(-5)\n" \
                                "}\n" \
                            "}\n" \
                        "},\n" \
                        "valueRange: [" + str(ts_min) + "," + str(ts_max) + "],\n" \
                        "ylabel: 'LOS displacement [cm]',\n" \
                        "yLabelWidth: 18,\n" \
                        "drawPoints: true,\n" \
                        "strokeWidth: 0,\n" \
                        "pointSize: 3,\n" \
                        "highlightCircleSize: 6,\n" \
                        "axisLabelFontSize: 12,\n" \
                        "xRangePad: 30,\n" \
                        "yRangePad: 30,\n" \
                        "hideOverlayOnMouseOut: false,\n" \
                        "panEdgeFraction: 0.0\n" \
                      "});\n" \
                      "</script>"
    return js_data_string


def create_regionalized_networklinks_file(regions, ts_obj, kml_data_files_directory, links_directory, output_file):

    # 5.4 create directory to store data KML files
    cmdDirectory = "cd {}; mkdir {}/".format(kml_data_files_directory, links_directory)
    print("Creating KML links directory")
    os.system(cmdDirectory)

    kml = KML.kml()

    kml_document = KML.Document()

    box_list = split_into_sub_boxes((ts_obj.length, ts_obj.width))

    region_nums = list(range(0, len(regions)))

    if "large" in links_directory:
        minLod = 1500
        maxLod = -1
    else:
        minLod = 0
        maxLod = 1500

    for ele in zip(regions, box_list, region_nums):

        box = ele[1]
        region_document = ele[0]
        num = ele[2]

        kml_file = "region_{}.kml".format(num)

        kml_1 = KML.kml()
        kml_1.append(region_document)
        print('writing ' + kml_file)
        with open(kml_file, 'w') as f:
            f.write(etree.tostring(kml_1, pretty_print=True).decode('utf-8'))

        cmdMove = "mv {sm} {dir}/{sm}".format(sm=kml_file, dir="{}/{}".format(kml_data_files_directory, links_directory))
        print("Moving KML Data Files to directory")
        os.system(cmdMove)

        lats, lons = get_lat_lon(ts_obj.metadata, box=box)
        lats = sorted(lats.flatten())
        lons = sorted(lons.flatten())

        network_link = KML.NetworkLink(
            KML.name('Region {}'.format(num)),
            KML.visibility(1),
            KML.Region(
                KML.Lod(
                    KML.minLodPixels(minLod),
                    KML.maxLodPixels(maxLod)
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

        kml_document.append(network_link)

    kml.append(kml_document)
    #print('writing ' + kml_file)
    with open(output_file, 'w') as f:
        f.write(etree.tostring(kml, pretty_print=True).decode('utf-8'))

    cmdMove = "mv {sm} {dir}/{sm};".format(sm=output_file, dir=kml_data_files_directory)
    print("Moving KML Data Files to directory")
    os.system(cmdMove)

    return output_file




def main(iargs=None):

    inps = cmd_line_parse(iargs)

    out_name_base = plot.auto_figure_title(inps.ts_file, inps_dict=vars(inps))

    kml_data_files_directory = "kml_data_files"
    kml_file_master = '{}_master.kml'.format(out_name_base)
    kml_file_sm = '{}_sm.kml'.format(out_name_base)
    kml_file_lg = '{}_lg.kml'.format(out_name_base)
    kmz_file = '{}.kmz'.format(out_name_base)

    cbar_png_file = '{}_cbar.png'.format(out_name_base)
    dygraph_file = "dygraph-combined.js"
    dot_file = "shaded_dot.png"
    star_file = "star.png"

    # 5.4 create directory to store data KML files
    cmdDirectory = "mkdir {}/".format(kml_data_files_directory)
    print("Creating KML Data File directory")
    os.system(cmdDirectory)

    small_dset_step = 20    # Increase this for coarser resolution in small dset
    large_dset_step = 3     # Decrease this for finer resolution in large dset

    ## 1. read data
    ts_obj = timeseries(inps.ts_file)
    ts_obj.open()

    lats, lons = get_lat_lon(ts_obj.metadata)

    vel = readfile.read(inps.vel_file, datasetName='velocity')[0] * 100.
    # Set min/max velocity for colormap
    if inps.vlim is None:
        inps.vlim = [np.nanmin(vel), np.nanmax(vel)]

    ## 2. Generate large and small KML files for different view heights
    kml_documents_sm = create_kml_document(inps, small_dset_step, dot_file, dygraph_file)
    kml_documents_lg = create_kml_document(inps, large_dset_step, dot_file, dygraph_file)

    ## 3. Create master KML file with network links to data KML files
    kml_master = KML.kml()

    kml_master_document = KML.Document()

    # 3.1 Create Screen Overlay element for colorbar
    cbar_png_file = plot_colorbar(out_file=cbar_png_file, vmin=inps.vlim[0], vmax=inps.vlim[1], cmap=inps.colormap)

    legend_overlay = KML.ScreenOverlay(
        KML.name('Legend'),
        KML.Icon(
            KML.href("{}".format(cbar_png_file)),
            KML.viewBoundScale(0.75)
        ),
        KML.overlayXY(x="0", y="0", xunits="pixels", yunits="insetPixels"),
        KML.screenXY(x="0", y="0", xunits="pixels", yunits="insetPixels"),
        KML.size(x="300", y="0", xunits="pixel", yunits="pixel"),
        KML.rotation(0),
        KML.visibility(1),
        KML.open(0)
    )
    print('add legend.')
    kml_master_document.append(legend_overlay)

    # 3.2 Generate the placemark for the Reference Pixel
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
    print('add reference point.')
    reference_folder = KML.Folder(KML.name("ReferencePoint"))
    reference_folder.append(reference_point)
    kml_master_document.append(reference_folder)

    lats = sorted(lats.flatten())
    lons = sorted(lons.flatten())

    data_folder = KML.Folder(KML.name("Data"))

    create_regionalized_networklinks_file(kml_documents_sm, ts_obj, kml_data_files_directory, "small_links_dir", kml_file_sm)
    create_regionalized_networklinks_file(kml_documents_lg, ts_obj, kml_data_files_directory, "large_links_dir", kml_file_lg)

    # 3.4 Create network link to small KML file
    network_link_sm = KML.NetworkLink(
                            KML.name('20 by 20'),
                            KML.visibility(1),
                            KML.Region(
                                KML.Lod(
                                    KML.minLodPixels(0),
                                    KML.maxLodPixels(1800)
                                ),
                                KML.LatLonAltBox(
                                    KML.north(lats[-1]+0.5),
                                    KML.south(lats[0]-0.5),
                                    KML.east(lons[-1]+0.5),
                                    KML.west(lons[0]-0.5)
                                )
                            ),
                            KML.Link(
                                KML.href("{}/{}".format(kml_data_files_directory, kml_file_sm)),
                                KML.viewRefreshMode('onRegion')
                            )
                      )

    # 3.4 Create network link to large KML file
    network_link_lg = KML.NetworkLink(
                            KML.name('3 by 3'),
                            KML.visibility(1),
                            KML.Region(
                                KML.Lod(
                                    KML.minLodPixels(1800),
                                    KML.maxLodPixels(-1)
                                ),
                                KML.LatLonAltBox(
                                    KML.north(lats[-1] + 0.5),
                                    KML.south(lats[0] - 0.5),
                                    KML.east(lons[-1] + 0.5),
                                    KML.west(lons[0] - 0.5)
                                )
                            ),
                            KML.Link(
                                KML.href("{}/{}".format(kml_data_files_directory, kml_file_lg)),
                                KML.viewRefreshMode('onRegion')
                            )
                        )

    data_folder.append(network_link_sm)
    data_folder.append(network_link_lg)

    kml_master_document.append(data_folder)

    kml_master.append(kml_master_document)


    ## 4 Write KML files

    # 4.3 write master KML file
    print('writing ' + kml_file_master)
    with open(kml_file_master, 'w') as f:
        f.write(etree.tostring(kml_master, pretty_print=True).decode('utf-8'))


    ## 5 Copy auxiliary files
    # 5.1 shaded_dot file
    dot_path = os.path.join(os.path.dirname(__file__), "utils/resources", dot_file)
    cmdDot = "cp {} {}".format(dot_path, dot_file)
    print("copying {} for point.".format(dot_file))
    os.system(cmdDot)

    # 5.2 star file
    star_path = os.path.join(os.path.dirname(__file__), "utils/resources", star_file)
    cmdStar = "cp {} {}".format(star_path, star_file)
    print("copying {} for reference point.".format(star_file))
    os.system(cmdStar)

    # 5.3 dygraph-combined.js file
    dygraph_path = os.path.join(os.path.dirname(__file__), "utils/resources", dygraph_file)
    cmdDygraph = "cp {} {}".format(dygraph_path, dygraph_file)
    print("copying {} for interactive plotting.".format(dygraph_file))
    os.system(cmdDygraph)

    # 5.5 move data KML files into new directory
    cmdMove = "mv {sm} {dir}/{sm}; mv {lg} {dir}/{lg};".format(sm=kml_file_sm, dir=kml_data_files_directory, lg=kml_file_lg)
    print("Moving KML Data Files to directory")
    os.system(cmdMove)

    ## 6 Generate KMZ file
    cmdKMZ = 'zip -r {} {} {} {} {} {} {}'.format(kmz_file, kml_data_files_directory, kml_file_master, cbar_png_file, dygraph_file, dot_file, star_file)
    print('writing {}\n{}'.format(kmz_file, cmdKMZ))
    os.system(cmdKMZ)

    cmdClean = 'rm -r {} {} {} {} {} {}'.format(kml_data_files_directory, kml_file_master, cbar_png_file, dygraph_file, dot_file, star_file)
    os.system(cmdClean)

    print('Done.')
    return


######################################################################################
if __name__ == '__main__':
    main()
