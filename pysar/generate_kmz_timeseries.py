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


def get_lat_lon(meta, mask=None):
    """extract lat/lon info of all grids into 2D matrix
    Parameters: meta : dict, including X/Y_FIRST/STEP and LENGTH/WIDTH info
    Returns:    lats : 2D np.array for latitude  in size of (length, width)
                lons : 2D np.array for longitude in size of (length, width)
    """
    # generate 2D matrix for lat/lon
    lat_num = int(meta['LENGTH'])
    lon_num = int(meta['WIDTH'])
    lat0 = float(meta['Y_FIRST'])
    lon0 = float(meta['X_FIRST'])
    lat_step = float(meta['Y_STEP'])
    lon_step = float(meta['X_STEP'])
    lat1 = lat0 + lat_step * lat_num
    lon1 = lon0 + lon_step * lon_num
    lats, lons = np.mgrid[lat0:lat1:lat_num*1j,
                          lon0:lon1:lon_num*1j]
    return lats, lons


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
    hex = mpl.colors.to_hex([rgba[3], rgba[2],
                             rgba[1], rgba[0]],
                            keep_alpha=True)[1:]
    return hex


def main(iargs=None):

    inps = cmd_line_parse(iargs)

    out_name_base = plot.auto_figure_title(inps.ts_file, inps_dict=vars(inps))
    cbar_png_file = '{}_cbar.png'.format(out_name_base)
    kml_file = '{}.kml'.format(out_name_base)
    kmz_file = '{}.kmz'.format(out_name_base)

    dygraph_file = "dygraph-combined.js"
    dot_file = "shaded_dot.png"
    star_file = "star.png"


    ## 1. read data
    ts_obj = timeseries(inps.ts_file)
    ts_obj.open()
    length, width = ts_obj.length, ts_obj.width

    # 1.1 Date
    dates = np.array(ts_obj.times)  # 1D np.array of dates in datetime.datetime object in size of [num_date,]
    dates = list(map(lambda d: d.strftime("%Y-%m-%d"), dates))
    num_date = len(dates)

    # 1.2 Spatial coordinates
    lats, lons = get_lat_lon(ts_obj.metadata)
    rows, cols = np.mgrid[0:length-1:length*1j,
                          0:width-1:width*1j]

    # 1.3 Velocity / time-series
    print('read velocity data')
    vel = readfile.read(inps.vel_file, datasetName='velocity')[0] * 100.
    vel_std = readfile.read(inps.vel_file, datasetName='velocityStd')[0] * 100.
    print('read time-series data')
    ts_data = readfile.read(inps.ts_file)[0] * 100.
    ts_data -= np.tile(ts_data[0, :, :], (ts_data.shape[0], 1, 1))     #enforce displacement starts from zero
    print('read temporal coherence data')
    temp_coh = readfile.read(inps.tcoh_file)[0]
    mask = ~np.isnan(vel)

    # data stats
    ts_min = np.nanmin(ts_data)
    ts_max = np.nanmax(ts_data)
    # Set min/max velocity for colormap
    if inps.vlim is None:
        inps.vlim = [np.nanmin(vel), np.nanmax(vel)]


    ## 2. Create KML Document
    print('create KML file.')
    kml_document = KML.Document()

    # 2.1 Create Screen Overlay element for colorbar
    cbar_png_file = plot_colorbar(out_file=cbar_png_file, vmin=inps.vlim[0], vmax=inps.vlim[1], cmap=inps.colormap)

    legend_overlay = KML.ScreenOverlay(
        KML.name('Legend'),
        KML.Icon(
            KML.href(cbar_png_file),
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
    kml_document.append(legend_overlay)

    colormap = mpl.cm.get_cmap(inps.colormap)                    # set colormap
    norm = mpl.colors.Normalize(vmin=inps.vlim[0], vmax=inps.vlim[1])


    # 2.2 Generate the placemark for the Reference Pixel
    ref_yx = (int(ts_obj.metadata['REF_Y']), int(ts_obj.metadata['REF_X']))
    ref_lalo = (lats[ref_yx[0], ref_yx[1]], lons[ref_yx[0], ref_yx[1]])

    reference_point =   KML.Placemark(
                            KML.Style(
                                KML.IconStyle(
                                    KML.color(get_color_for_velocity(0.0, colormap, norm)),
                                    KML.scale(1.),
                                    KML.Icon(KML.href(star_file))
                                )
                            ),
                            KML.description("Reference point <br /> \n <br /> \n"+\
                                generate_description_string(ref_lalo, ref_yx, 0.00, 0.00, 0.00, 1.00)
                            ),
                            KML.Point(
                                KML.coordinates("{}, {}".format(ref_lalo[1], ref_lalo[0]))
                            )
                        )
    print('add reference point.')
    reference_folder =  KML.Folder(KML.name("ReferencePoint"))
    reference_folder.append(reference_point)
    kml_document.append(reference_folder)


    # 2.3 Data folder for all points
    step = inps.step
    num_pixel = int(length/step) * int(width/step)
    msg = "add point time-series data "
    msg += "(select 1 out of every {} by {} pixels; select {} pixels in total) ...".format(step, step, num_pixel)
    print(msg)

    data_folder =   KML.Folder(KML.name("Data"))
    for i in range(0, length, step):
        for j in range(0, width, step):
            if mask[i,j]:          # add point if it's not marked as masked out
                lat = lats[i,j]
                lon = lons[i,j]
                row = rows[i,j]
                col = cols[i,j]
                ts = ts_data[:, i, j]
                v = vel[i,j]
                vstd = vel_std[i,j]
                tcoh = temp_coh[i,j]

                # Create KML icon style element                
                style = KML.Style(
                            KML.IconStyle(
                                KML.color(get_color_for_velocity(v, colormap, norm)),
                                KML.scale(0.5),
                                KML.Icon(KML.href(dot_file))
                            )
                        )

                # Create KML point element
                point = KML.Point(KML.coordinates("{},{}".format(lon, lat)))

                # Javascript to embed inside the description
                js_data_string = "<script type='text/javascript' src='"+dygraph_file+"'></script>\n" \
                                    "<div id='graphLegend' style='padding-left: 200px; margin-bottom: 10px; font-size: 16px'></div>\n" \
                                    "<div id='graphdiv'> </div>\n" \
                                    "<script type='text/javascript'>\n" \
                                        "g = new Dygraph( document.getElementById('graphdiv'),\n" \
                                        "\"Date, Displacement\\n\" + \n"

                # append the date/displacement data
                for k in range(num_date):
                    date = dates[k]
                    dis = ts[k]
                    date_displacement_string = "\"{}, {}\\n\" + \n".format(date, dis)  
                    js_data_string += date_displacement_string

                # TODO: Fix xRangePad and yRangePad issues
                js_data_string +=   "\"\",\n" \
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
                                                    "var months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']\n" \
                                                    "var date = new Date(d)\n" \
                                                    "var dateString = 'Date: ' + date.getDate() + ' ' + months[date.getMonth()] + ' ' + date.getFullYear()\n" \
                                                    "return dateString;\n" \
                                                "},\n" \
                                                "pixelsPerLabel: 90\n" \
                                            "},\n" \
                                            "y: {\n" \
                                                "valueFormatter: function(v) {\n" \
                                                    "return v.toFixed(5)\n" \
                                                "}\n" \
                                            "}\n" \
                                        "},\n" \
                                        "valueRange: ["+str(ts_min-2)+","+str(ts_max+2)+"],\n" \
                                        "ylabel: 'LOS displacement [cm]',\n" \
                                        "yLabelWidth: 18,\n" \
                                        "drawPoints: true,\n" \
                                        "strokeWidth: 0,\n" \
                                        "pointSize: 3,\n" \
                                        "highlightCircleSize: 6,\n" \
                                        "axisLabelFontSize: 12,\n" \
                                        "xRangePad: 30,\n" \
                                        "yRangePad: 30,\n" \
                                        "labelsDiv: 'graphLegend',\n" \
                                        "hideOverlayOnMouseOut: false\n" \
                                    "});\n" \
                                  "</script>"

                # Create KML description element
                stats_info = generate_description_string((lat, lon), (row, col), v, vstd, ts[-1], tcoh=tcoh)
                description = KML.description(stats_info, js_data_string)
        
                # Crate KML Placemark element to hold style, description, and point elements
                placemark = KML.Placemark(style, description, point)
                # Append each placemark element to the KML document object
                data_folder.append(placemark)

    kml_document.append(data_folder)

    # 2.4 Write KML file
    kml = KML.kml()
    kml.append(kml_document)
    print('writing ' + kml_file)
    with open(kml_file, 'w') as f:
        f.write(etree.tostring(kml, pretty_print=True).decode('utf-8'))

    # 2.5 Copy auxliary files
    # 2.5.1 shaded_dot file
    dot_path = os.path.join(os.path.dirname(__file__), "utils/resources", dot_file)
    cmdDot = "cp {} {}".format(dot_path, dot_file)
    print("copying {} for point.".format(dot_file))
    os.system(cmdDot)

    # 2.5.2 star file
    star_path = os.path.join(os.path.dirname(__file__), "utils/resources", star_file)
    cmdStar = "cp {} {}".format(star_path, star_file)
    print("copying {} for reference point.".format(star_file))
    os.system(cmdStar)

    # 2.5.3 dygraph-combined.js file
    dygraph_path = os.path.join(os.path.dirname(__file__), "utils/resources", dygraph_file)
    cmdDygraph = "cp {} {}".format(dygraph_path, dygraph_file)
    print("copying {} for interactive plotting.".format(dygraph_file))
    os.system(cmdDygraph)

    # 2.6 Generate KMZ file
    cmdKMZ = 'zip {} {} {} {} {} {}'.format(kmz_file, kml_file, cbar_png_file, dygraph_file, dot_file, star_file)
    print('writing {}\n{}'.format(kmz_file, cmdKMZ))
    os.system(cmdKMZ)

    cmdClean = 'rm {} {} {} {} {}'.format(kml_file, cbar_png_file, dygraph_file, dot_file, star_file)
    print(cmdClean)
    os.system(cmdClean)

    print('Done.')
    return


######################################################################################
if __name__ == '__main__':
    main()
