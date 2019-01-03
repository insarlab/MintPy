#!/usr/bin/env python3

import os
import numpy as np
from lxml import etree
from pykml.factory import KML_ElementMaker as KML
from pysar.objects import timeseries
from pysar.utils import readfile

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

    ts_file = "/Users/joshua/Desktop/pysar/test_data/FernandinaSenDT128/geo_timeseries_ECMWF_ramp_demErr_masked.h5"
    vel_file = "/Users/joshua/Desktop/pysar/test_data/FernandinaSenDT128/geo_velocity_masked.h5"

    ts_obj = timeseries(ts_file)
    ts_obj.open()

    ## read data

    # Date
    dates = np.array(ts_obj.times)  # 1D np.array of dates in datetime.datetime object in size of [num_date,]
    dates = list(map(lambda d: d.strftime("%Y-%m-%d"), dates))

    # Velocity / time-series
    vel = readfile.read(vel_file)[0]
    mask = ~np.isnan(vel)  # matrix indicating valid pixels (True for valid, False for invalid)
    vel = readfile.read(vel_file)[0][mask]  # 1D np.array of velocity   in np.float32 in size of [num_pixel,] in meter/year
    ts = readfile.read(ts_file)[0][:, mask]  # 2D np.array of time-sries in np.float32 in size of [num_date, num_pixel] in meter

    # Spatial coordinates
    lats = get_lat_lon(ts_obj.metadata)[0][mask]  # 1D np.array of latitude  in np.float32 in size of [num_pixel,] in degree
    lons = get_lat_lon(ts_obj.metadata)[1][mask]  # 1D np.array of longitude in np.float32 in size of [num_pixel,] in degree
    coords = list(zip(lats, lons))


    # Create KML Document
    kml = KML.kml()
    kml_document = KML.Document()

    for i in range(0, len(coords), 100):
        lat = coords[i][0]
        lon = coords[i][1]
        v = vel[i]
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
            displacement = ts[j][i]

            date_displacement_string = "\"{}, {}\\n\" + \n".format(date, displacement)

            js_data_string += date_displacement_string

        js_data_string +=   "{" \
                                "valueRange: [-100,100]," \
                                "ylabel: '[mm]'," \
                            "});" \
                          "</script>" \
                          "]]>"

        description = KML.description(js_data_string)

        placemark = KML.Placemark(description, point)

        kml_document.append(placemark)
    kml.append(kml_document)

    with open("/Users/joshua/Desktop/test.kml", 'w+') as kml_file:
        kml_file.write(etree.tostring(kml, pretty_print=True).decode('utf-8'))

    print("Done")

