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

    parser = argparse.ArgumentParser(description='Generare Google Earth Compatible KML for offline timeseries analysis',
                                     formatter_class=argparse.RawTextHelpFormatter)

    args = parser.add_argument_group('Input File', 'File/Dataset to display')

    args.add_argument('--ts', dest='ts_file', metavar='FILE', help='Timeseries file to generate KML for')
    args.add_argument('--vel', dest='vel_file', metavar='FILE', help='Velocity file')
    args.add_argument('-v','--vlim', dest='vlim', nargs=2, metavar=('VMIN', 'VMAX'), type=float,
                      help='Display limits for matrix plotting.')
    args.add_argument('-c', '--colormap', dest='colormap',
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


if __name__ == "__main__":

    inps = cmd_line_parse(sys.argv[1:])

    out_name_base = plot.auto_figure_title(inps.ts_file, inps_dict=vars(inps))

    ts_file = inps.ts_file
    vel_file = inps.vel_file

    ts_obj = timeseries(ts_file)
    ts_obj.open()

    ## read data

    # Date
    dates = np.array(ts_obj.times)  # 1D np.array of dates in datetime.datetime object in size of [num_date,]
    dates = list(map(lambda d: d.strftime("%Y-%m-%d"), dates))

    # Velocity / time-series
    vel = readfile.read(vel_file)[0]*100
    mask = ~np.isnan(vel)                   # matrix indicating valid pixels (True for valid, False for invalid)
    vel = readfile.read(vel_file)[0][mask]  # 1D np.array of velocity   in np.float32 in size of [num_pixel,] in meter/year
    ts = readfile.read(ts_file)[0][:, mask] # 2D np.array of time-sries in np.float32 in size of [num_date, num_pixel] in meter


    # Set min/max velocity for colormap
    if inps.vlim is None:
        min_vel = min(vel)
        max_vel = max(vel)
    else:
        min_vel = inps.vlim[0]
        max_vel = inps.vlim[1]

    norm = mpl.colors.Normalize(vmin=min_vel, vmax=max_vel)  # normalize velocity colors between 0.0 and 1.0

    if inps.colormap is None:
        cmap = "jet"
    else:
        cmap = inps.colormap

    # Spatial coordinates
    lats = get_lat_lon(ts_obj.metadata)[0][mask]  # 1D np.array of latitude  in np.float32 in size of [num_pixel,] in degree
    lons = get_lat_lon(ts_obj.metadata)[1][mask]  # 1D np.array of longitude in np.float32 in size of [num_pixel,] in degree
    coords = list(zip(lats, lons))

    # Create and save colorbar image
    pc = plt.figure(figsize=(8, 1))
    cax = pc.add_subplot(111)
    cbar = mpl.colorbar.ColorbarBase(cax, cmap=inps.colormap, norm=norm, orientation='horizontal')
    cbar.set_label('{} [{}]'.format("Velocity", "cm"))
    cbar.update_ticks()

    pc.patch.set_facecolor('white')
    pc.patch.set_alpha(0.7)

    cbar_png_file = '{}_cbar.png'.format(out_name_base)
    print('writing ' + cbar_png_file)
    pc.savefig(cbar_png_file, bbox_inches='tight', facecolor=pc.get_facecolor(), dpi=300)

    # Create KML Document
    kml = KML.kml()
    kml_document = KML.Document()

    for i in range(0, len(coords), 10):
        lat = coords[i][0]
        lon = coords[i][1]
        v = vel[i]

        colormap = mpl.cm.get_cmap(cmap)                    # set colormap
        rgba = colormap(norm(v))                            # get rgba color components for point velocity
        hex = mpl.colors.to_hex([rgba[3], rgba[2],
                                 rgba[1], rgba[0]],
                                keep_alpha=True)[1:]    # convert rgba to hex components reversed for kml color specification

        # Create KML icon style element
        style = KML.Style(
                    KML.IconStyle(
                        KML.color(hex),
                        KML.scale(0.5),
                        KML.Icon(
                            KML.href("http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png")
                        )
                    )
                )

        # Create KML point element
        point = KML.Point(
                    KML.coordinates("{},{}".format(lon, lat))
                )

        # Javascript to embed inside the description
        js_data_string = "< ![CDATA[\n" \
                            "<script type='text/javascript' src='Velocity__mm_year__time_series_dir/dygraph-combined.js'></script>\n" \
                            "<div id='graphdiv'> </div>\n" \
                            "<script type='text/javascript'>\n" \
                                "g = new Dygraph( document.getElementById('graphdiv'),\n" \
                                "\"Date, Displacement\\n\" + \n"

        for j in range(len(dates)):
            date = dates[j]
            displacement = ts[j][i]*100

            date_displacement_string = "\"{}, {}\\n\" + \n".format(date, displacement)  # append the date/displacement data

            js_data_string += date_displacement_string

        js_data_string +=       "{" \
                                    "valueRange: [-100,100]," \
                                    "ylabel: '[mm]'," \
                                "});" \
                              "</script>" \
                          "]]>"

        # Create KML description element
        description = KML.description(js_data_string)

        # Crate KML Placemark element to hold style, description, and point elements
        placemark = KML.Placemark(style, description, point)

        kml_document.append(placemark)

    kml.append(kml_document)

    # Write KML file
    kml_file = '{}.kml'.format(out_name_base)
    print('writing ' + kml_file)
    with open(kml_file, 'w') as f:
        f.write(etree.tostring(kml, pretty_print=True).decode('utf-8'))

    # 2.4 Generate KMZ file
    kmz_file = '{}.kmz'.format(out_name_base)
    cmdKMZ = 'zip {} {} {}'.format(kmz_file, kml_file, cbar_png_file)
    print('writing {}\n{}'.format(kmz_file, cmdKMZ))
    os.system(cmdKMZ)

    cmdClean = 'rm {} {}'.format(kml_file, cbar_png_file)
    print(cmdClean)
    os.system(cmdClean)

    # # Write KML structure to .kml file
    # with open("/Users/joshua/Desktop/test.kml", 'w+') as kml_file:
    #     kml_file.write(etree.tostring(kml, pretty_print=True).decode('utf-8'))
    #
    # print("Done")

