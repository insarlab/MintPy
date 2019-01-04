#!/usr/bin/env python3

import sys
import argparse
import numpy as np
from lxml import etree
from pykml.factory import KML_ElementMaker as KML
from pysar.objects import timeseries
from pysar.utils import readfile

def create_parser():

    parser = argparse.ArgumentParser(description='Generare Google Earth Compatible KML for offline timeseries analysis',
                                     formatter_class=argparse.RawTextHelpFormatter)

    args = parser.add_argument_group('Input File', 'File/Dataset to display')

    args.add_argument('--ts', dest='ts_file', metavar='FILE', help='Timeseries file to generate KML for')
    args.add_argument('--vel', dest='vel_file', metavar='FILE', help='Velocity file')
    args.add_argument('--v','--vlim', dest='vlim', nargs=2, metavar=('VMIN', 'VMAX'), type=float,
                      help='Display limits for matrix plotting.')

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

    print(inps.vlim)

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

    if inps.vlim is None:
        min_vel = min(vel)
        max_vel = max(vel)
    else:
        min_vel = inps.vlim[0]
        max_vel = inps.vlim[1]

    print("({}, {})".format(min_vel, max_vel))

    vel_range = max_vel - min_vel
    vel_step = vel_range / 51

    vel_color_map = {}

    colors = ["ff00009f", "ff0000cf", "ff0000df", "ff0000ef", "ff0000ff", "ff000fff", "ff001fff", "ff002fff", "ff003fff",
              "ff004fff", "ff005fff", "ff006fff", "ff007fff", "ff008fff", "ff009fff", "ff00afff", "ff00bfff", "ff00cfff",
              "ff00dfff", "ff00efff", "ff00ffff", "ff0fffef", "ff1fffdf", "ff2fffcf", "ff3fffbf", "ff4fffaf", "ff5fff9f",
              "ff6fff8f", "ff7fff7f", "ff8fff6f", "ff9fff5f", "ffafff4f", "ffbfff3f", "ffcfff2f", "ffdfff1f", "ffefff0f",
              "ffffff00", "ffffef00", "ffffdf00", "ffffcf00", "ffffbf00", "ffffaf00", "ffff9f00", "ffff8f00", "ffff7f00",
              "ffff6f00", "ffff5f00", "ffff4f00", "ffff3f00", "ffff0000", "ffef0000"]

    for i in range(51):
        v_range = (min_vel+(vel_step*i), min_vel+(vel_step * (i+1)))
        color = colors[i]

        vel_color_map[v_range] = color

    # Spatial coordinates
    lats = get_lat_lon(ts_obj.metadata)[0][mask]  # 1D np.array of latitude  in np.float32 in size of [num_pixel,] in degree
    lons = get_lat_lon(ts_obj.metadata)[1][mask]  # 1D np.array of longitude in np.float32 in size of [num_pixel,] in degree
    coords = list(zip(lats, lons))


    # Create KML Document
    kml = KML.kml()
    kml_document = KML.Document()

    for i in range(0, len(coords), 10):
        lat = coords[i][0]
        lon = coords[i][1]
        v = vel[i]

        color = None

        for key in vel_color_map:
            if key[0] <= v <= key[1]:
                color = vel_color_map[key]
                break

        style = KML.Style(
                    KML.IconStyle(
                        KML.color(color),
                        KML.scale(0.5),
                        KML.Icon(
                            KML.href("http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png")
                        )
                    )
                )

        point = KML.Point(
                    KML.coordinates("{},{}".format(lon, lat))
                )

        js_data_string = "< ![CDATA[\n" \
                            "<script type='text/javascript' src='Velocity__mm_year__time_series_dir/dygraph-combined.js'></script>\n" \
                            "<div id='graphdiv'> </div>\n" \
                            "<script type='text/javascript'>\n" \
                                "g = new Dygraph( document.getElementById('graphdiv'),\n" \
                                "\"Date, Displacement\\n\" + \n"
        for j in range(len(dates)):
            date = dates[j]
            displacement = ts[j][i]*100

            date_displacement_string = "\"{}, {}\\n\" + \n".format(date, displacement)

            js_data_string += date_displacement_string

        js_data_string +=   "{" \
                                "valueRange: [-100,100]," \
                                "ylabel: '[mm]'," \
                            "});" \
                          "</script>" \
                          "]]>"

        description = KML.description(js_data_string)

        placemark = KML.Placemark(style, description, point)

        kml_document.append(placemark)
    kml.append(kml_document)

    with open("/Users/joshua/Desktop/test.kml", 'w+') as kml_file:
        kml_file.write(etree.tostring(kml, pretty_print=True).decode('utf-8'))

    print("Done")

