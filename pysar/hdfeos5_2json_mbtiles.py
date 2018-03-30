#!/usr/bin/env python3

import json
import h5py
import numpy as np
from datetime import date
import math
import time
import os
import sys
import psycopg2
import geocoder
from pysar.add_attribute_insarmaps import InsarDatabaseController
from pysar.mask import mask_matrix
import argparse
import pickle

# ex: python Converter_unavco.py Alos_SM_73_2980_2990_20070107_20110420.h5

# This script takes a UNAVCO format timeseries h5 file, converts to mbtiles, 
# and sends to database which allows website to make queries and display data
# ---------------------------------------------------------------------------------------
# FUNCTIONS
# ---------------------------------------------------------------------------------------
# returns a dictionary of datasets that are stored in memory to speed up h5 read process
def get_date(date_string): 
    year = int(date_string[0:4])
    month = int(date_string[4:6])
    day = int(date_string[6:8])
    return date(year, month, day)
# ---------------------------------------------------------------------------------------
# takes a date and calculates the number of days elapsed in the year of that date
# returns year + (days_elapsed / 365), a decimal representation of the date necessary
# for calculating linear regression of displacement vs time
def get_decimal_date(d):
    start = date(d.year, 1, 1)
    return abs(d-start).days / 365.0 + d.year

def region_name_from_project_name(project_name):
    track_index = project_name.find('T')

    return project_name[:track_index]

needed_attributes = {
    "prf", "first_date", "mission", "WIDTH", "X_STEP", "processing_software",
    "wavelength", "processing_type", "beam_swath", "Y_FIRST", "look_direction",
    "flight_direction", "last_frame", "post_processing_method", "min_baseline_perp"
    "unwrap_method", "relative_orbit", "beam_mode", "LENGTH", "max_baseline_perp",
    "X_FIRST", "atmos_correct_method", "last_date", "first_frame", "frame", "Y_STEP", "history",
    "scene_footprint", "data_footprint", "downloadUnavcoUrl", "referencePdfUrl", "areaName", "referenceText",
    "ref_lat", "ref_lon"
}

def serialize_dictionary(dictionary, fileName):
    with open(fileName, "w") as file:
        pickle.dump(dictionary, file)
# ---------------------------------------------------------------------------------------
# convert h5 file to json and upload it. folder_name == unavco_name
def convert_data(attributes, decimal_dates, timeseries_datasets, dates, json_path, folder_name):

    region_file = None
    project_name = attributes["PROJECT_NAME"]
    region = region_name_from_project_name(project_name)
# get the attributes for calculating latitude and longitude
    x_step = float(attributes["X_STEP"])
    y_step = float(attributes["Y_STEP"])
    x_first = float(attributes["X_FIRST"])
    y_first = float(attributes["Y_FIRST"])
    num_columns = int(attributes["WIDTH"])
    num_rows = int(attributes["LENGTH"])
    print("columns: %d" % num_columns)
    print("rows: %d" % num_rows)
    # create a siu_man array to store json point objects
    siu_man = []
    displacement_values = []
    displacements = '{'
    # np array of decimal dates, x parameter in linear regression equation
    x = decimal_dates
    A = np.vstack([x, np.ones(len(x))]).T
    y = []
    chunk_num = 1
    point_num = 0
    CHUNK_SIZE = 20000

    # iterate through h5 file timeseries
    for (row, col), value in np.ndenumerate(timeseries_datasets[dates[0]]):
        longitude = x_first + (col * x_step)
        latitude = y_first + (row * y_step) 
        displacement = float(value) 
        # if value is not equal to naN, create a new json point object and append to siu_man array
        if not math.isnan(displacement):
            # get displacement values for all the dates into array for json and string for pgsql
            for date in dates:
                displacement = timeseries_datasets[date][row][col]
                displacements += (str(displacement) + ",")
                displacement_values.append(float(displacement))
            displacements = displacements[:len(displacements) - 1] + '}'

            # np array of displacement values, y parameter in linear regression equation
            y = displacement_values

            # y = mx + c -> we want m = slope of the linear regression line 
            m, c = np.linalg.lstsq(A, y)[0]

            data = {
            "type": "Feature",
            "geometry": {"type": "Point", "coordinates": [longitude, latitude]},    
            "properties": {"d": displacement_values, "m": m, "p": point_num}
            }   

            siu_man.append(data)

            # clear displacement array for json and the other string for dictionary, for next point
            displacement_values = []
            displacements = '{'
            point_num += 1
            # break;    # for testing purposes convert only 1 point

            # if chunk_size limit is reached, write chunk into a json file
            # then increment chunk number and clear siu_man array
            if len(siu_man) == CHUNK_SIZE:
                make_json_file(chunk_num, siu_man, dates, json_path, folder_name)
                chunk_num += 1
                siu_man = []

    # write the last chunk that might be smaller than chunk_size
    make_json_file(chunk_num, siu_man, dates, json_path, folder_name)

    # dictionary to contain metadata needed by db to be written to a file
    # and then be read by json_mbtiles2insarmaps.py
    insarmapsMetadata = {}
    # calculate mid lat and long of dataset - then use google python lib to get country
    mid_long = x_first + ((num_columns/2) * x_step)
    mid_lat = y_first + ((num_rows/2) * y_step)
    country = None
    try:
        g = geocoder.google([mid_lat,mid_long], method='reverse', timeout=60.0)
        country = str(g.country_long)
    except Exception as e:
        sys.stderr.write("timeout reverse geocoding country name")

    area = folder_name

    # for some reason pgsql only takes {} not [] - format date arrays and attributes to be inserted to pgsql
    string_dates_sql = '{'
    for k in dates:
        string_dates_sql += (str(k) + ",")
    string_dates_sql = string_dates_sql[:len(string_dates_sql) - 1] + '}'

    decimal_dates_sql = '{'
    for d in decimal_dates:
        decimal_dates_sql += (str(d) + ",")
    decimal_dates_sql = decimal_dates_sql[:len(decimal_dates_sql) - 1] + '}'
    # add keys and values to area table. TODO: this will be removed eventually
    # and all attributes will be put in extra_attributes table
    attribute_keys = '{'
    attribute_values = '{'
    for k in attributes:
        v = attributes[k]
        if k in needed_attributes:
            print(str(k) + ": " + str(v))
            attribute_keys += (str(k) + ",")
            attribute_values += (str(v) + ',')
    attribute_keys = attribute_keys[:len(attribute_keys)-1] + '}'
    attribute_values = attribute_values[:len(attribute_values)-1] + '}'

    # write out metadata to json file
    insarmapsMetadata["area"] = area
    insarmapsMetadata["project_name"] = project_name
    insarmapsMetadata["mid_long"] = mid_long
    insarmapsMetadata["mid_lat"] = mid_lat
    insarmapsMetadata["country"] = country
    insarmapsMetadata["region"] = region
    insarmapsMetadata["chunk_num"] = 1
    insarmapsMetadata["attribute_keys"] = attribute_keys
    insarmapsMetadata["attribute_values"] = attribute_values
    insarmapsMetadata["string_dates_sql"] = string_dates_sql
    insarmapsMetadata["decimal_dates_sql"] = decimal_dates_sql
    insarmapsMetadata["attributes"] = attributes
    insarmapsMetadata["needed_attributes"] = needed_attributes
    metadataFilePath = json_path + "/metadata.pickle" 
    serialize_dictionary(insarmapsMetadata, metadataFilePath)
# ---------------------------------------------------------------------------------------
# create a json file out of siu man array
# then put json file into directory named after the h5 file
def make_json_file(chunk_num, points, dates, json_path, folder_name):

    data = {
    "type": "FeatureCollection",
    "dates": dates,
    "features": points
    }

    chunk = "chunk_" + str(chunk_num) + ".json"
    json_file = open(json_path + "/" + chunk, "w")
    string_json = json.dumps(data, indent=4, separators=(',',':'))
    json_file.write("%s" % string_json)
    json_file.close()

    print("converted chunk " + str(chunk_num))

# ---------------------------------------------------------------------------------------
def build_parser():
    dbHost = "insarmaps.rsmas.miami.edu"
    parser = argparse.ArgumentParser(description='Convert a Unavco format H5 file for ingestion into insarmaps.')
    required = parser.add_argument_group("required arguments")
    required.add_argument("file", help="unavco file to ingest")
    required.add_argument("outputDir", help="directory to place json files and mbtiles file")

    return parser

# ---------------------------------------------------------------------------------------
# START OF EXECUTABLE
# ---------------------------------------------------------------------------------------
def main():
    parser = build_parser()
    parseArgs = parser.parse_args()
    file_name = parseArgs.file
    output_folder = parseArgs.outputDir
    should_mask = True

    path_name_and_extension = os.path.basename(file_name).split(".")
    path_name = path_name_and_extension[0]
    extension = path_name_and_extension[1]
# ---------------------------------------------------------------------------------------
# start clock to track how long conversion process takes
    start_time = time.clock()

# use h5py to open specified group(s) in the h5 file 
# then read datasets from h5 file into memory for faster reading of data
    file = h5py.File(file_name,  "r")
    timeseries_group = file["HDFEOS"]["GRIDS"]["timeseries"]
    displacement_3d_matrix = timeseries_group["observation"]["displacement"]

# get attributes (stored at root) of UNAVCO timeseries file
    attributes = dict(file.attrs)

# in timeseries displacement_3d_matrix, there are datasets
# need to get datasets with dates - strings that can be converted to integers
    dates = displacement_3d_matrix.attrs["DATE_TIMESERIES"].split(" ")

# array that stores dates from dates that have been converted to decimal
    decimal_dates = []

# read datasets in the group into a dictionary of 2d arrays and intialize decimal dates
    timeseries_datasets = {}
    i = 0
    for displacement_2d_matrix in displacement_3d_matrix:
        dataset = displacement_2d_matrix[:]
        if should_mask:
            print("Masking " + dates[i])
            mask = timeseries_group["quality"]["mask"][:]
            dataset = mask_matrix(dataset, mask)

        timeseries_datasets[dates[i]] = dataset
        d = get_date(dates[i])
        decimal = get_decimal_date(d)
        decimal_dates.append(decimal)
        i += 1

# close h5 file
    file.close()

    path_list = path_name.split("/")
    folder_name = path_name.split("/")[len(path_list)-1]

    try: # create path for output
        os.mkdir(output_folder)
    except:
        print(output_folder + " already exists")

# read and convert the datasets, then write them into json files and insert into database
    convert_data(attributes, decimal_dates, timeseries_datasets, dates, output_folder, folder_name)

# run tippecanoe command to get mbtiles file
    os.chdir(os.path.abspath(output_folder))
    os.system("tippecanoe *.json -l chunk_1 -x d -pf -pk -Bg -d9 -D12 -g12 -r0 -o " + folder_name + ".mbtiles")

# ---------------------------------------------------------------------------------------
# check how long it took to read h5 file data and create json files
    end_time =  time.clock()
    print(("time elapsed: " + str(end_time - start_time)))
# ---------------------------------------------------------------------------------------

if __name__ == '__main__':
    main()
