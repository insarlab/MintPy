#! /usr/bin/env python

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
import getopt
from pysar.add_attributes_insarmaps import InsarDatabaseController
import argparse

# ex: python Converter_unavco.py Alos_SM_73_2980_2990_20070107_20110420.h5

# This script takes a UNAVCO format timeseries h5 file, converts to mbtiles, 
# and sends to database which allows website to make queries and display data
# ---------------------------------------------------------------------------------------
# FUNCTIONS
# ---------------------------------------------------------------------------------------
dbUsername = "INSERT"
dbPassword = "INSERT"
dbHost = "INSERT"
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
# ---------------------------------------------------------------------------------------
def convert_data(): 
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

    # outer loop increments row = longitude, inner loop increments column = latitude
    for (row, col), value in np.ndenumerate(timeseries_datasets[dataset_keys[0]]):
        latitude = x_first + (col * x_step)
        longitude = y_first + (row * y_step) 
        displacement = float(value) 
        # if value is not equal to naN, create a new json point object and append to siu_man array
        if not math.isnan(displacement):
            # get displacement values for all the dates into array for json and string for pgsql
            for key in dataset_keys:
                displacement = timeseries_datasets[key][row][col]
                displacements += (str(displacement) + ",")
                displacement_values.append(float(displacement))
            displacements = displacements[:len(displacements) - 1] + '}'

            # np array of displacement values, y parameter in linear regression equation
            y = displacement_values

            # y = mx + c -> we want m = slope of the linear regression line 
            m, c = np.linalg.lstsq(A, y)[0]

            data = {
            "type": "Feature",
            "geometry": {"type": "Point", "coordinates": [latitude, longitude]},    
            "properties": {"d": displacement_values, "m": m, "p": point_num}
            }   
            # allocate memory space for siu_man array in beginning 
            siu_man.append(data)

            # clear displacement array for json and the other string for dictionary, for next point
            displacement_values = []
            displacements = '{'
            point_num += 1
            # break;    # for testing purposes convert only 1 point

            # if chunk_size limit is reached, write chunk into a json file
            # then increment chunk number and clear siu_man array
            if len(siu_man) == chunk_size:
                make_json_file(chunk_num, siu_man)
                chunk_num += 1
                siu_man = []

    # write the last chunk that might be smaller than chunk_size
    make_json_file(chunk_num, siu_man)

    # calculate mid lat and long of dataset - then use google python lib to get country
    mid_lat = x_first + ((num_columns/2) * x_step)
    mid_long = y_first + ((num_rows/2) * y_step)
    country = None
    try:
        g = geocoder.google([mid_long,mid_lat], method='reverse', timeout=60.0)
        country = str(g.country_long)
    except Exception, e:
        print "timeout reverse geocoding country name"

    area = folder_name

    # for some reason pgsql only takes {} not [] - format date arrays and attributes to be inserted to pgsql
    string_dates_sql = '{'
    for k in dataset_keys:
        string_dates_sql += (str(k) + ",")
    string_dates_sql = string_dates_sql[:len(string_dates_sql) - 1] + '}'

    decimal_dates_sql = '{'
    for d in decimal_dates:
        decimal_dates_sql += (str(d) + ",")
    decimal_dates_sql = decimal_dates_sql[:len(decimal_dates_sql) - 1] + '}'

    # scene_footprint attribute uses a wkt geometry type with format that confuses postgresql database
    # thus we have to add "Polygon(coordinates, coordinates, coordinates, coordinates)" as string
    attribute_keys = '{'
    attribute_values = '{'
    for k, v in attributes_dictionary:
        attribute_keys += (str(k) + ",")
        if "POLYGON" in str(v):
            arr = v.split(",")
            s = "\,"
            v = s.join(arr)
        attribute_values += (str(v) + ',')
    attribute_keys = attribute_keys[:len(attribute_keys)-1] + '}'
    attribute_values = attribute_values[:len(attribute_values)-1] + '}'

    try:    # connect to databse
        con = psycopg2.connect("dbname='pgis' user='" + dbUsername + "' host='" + dbHost + "' password='" + dbPassword + "'")
        cur = con.cursor()
        # create area table if not exist - limit for number of dates is 200, limt for number of attribute keys/values is 100
        cur.execute("CREATE TABLE IF NOT EXISTS area ( unavco_name varchar, project_name varchar, latitude double precision, longitude double precision, country varchar, region varchar, numchunks integer, attributekeys varchar[100], attributevalues varchar[100], stringdates varchar[200], decimaldates double precision[200] );")
        con.commit()
        print 'created area table'
    except Exception, e:
        print "unable to connect to the database"
        print e
        sys.exit()

    # put dataset into area table
    # area_data = {"latitude": mid_lat, "longitude": mid_long, "country": country, "num_chunks": chunk_num, "dates": dataset_keys}
    try:
        con = psycopg2.connect("dbname='pgis' user='" + dbUsername + "' host='" + dbHost + "' password='" + dbPassword + "'")
        cur = con.cursor()
        query = "INSERT INTO area VALUES (" + "'" + area + "','" + project_name + "','" + str(mid_lat) + "','" + str(mid_long) + "','" + country + "','" + region + "','" + str(chunk_num) + "','" + attribute_keys + "','" + attribute_values + "','" + string_dates_sql + "','" + decimal_dates_sql + "')"
        cur.execute(query)
        con.commit()
        con.close()
    except Exception, e:
        print "error inserting into area"
        print e
        sys.exit()

    # put attributes in own table. TODO: remove old way of adding attributes
    # via array
    attributesController = InsarDatabaseController(dbUsername, dbPassword, dbHost, 'pgis')
    attributesController.connect()

    for k, v in attributes_dictionary:
        if "POLYGON" in str(v):
            arr = v.split(",")
            s = "\,"
            v = s.join(arr)
        attributesController.add_attribute(project_name, k, v)
    attributesController.close()

    # create index to speed up queries:
    print "Creating index"
    try:
        con = psycopg2.connect("dbname='pgis' user='" + dbUsername + "' host='" + dbHost + "' password='" + dbPassword + "'")
        cur = con.cursor()
        query = 'CREATE INDEX ON ' + area + ' (p)'
        cur.execute(query)
        con.commit()
        con.close()
    except Exception, e:
        print "error creating index on p"
        print e
        #sys.exit() Don't exit, as what if the index already exists
    print "Done creating index"
    
# ---------------------------------------------------------------------------------------
# create a json file out of siu man array
# then put json file into directory named after the h5 file
def make_json_file(chunk_num, points): 

    data = {
    "type": "FeatureCollection",
    "dates": dataset_keys, 
    "features": points
    }

    chunk = "chunk_" + str(chunk_num) + ".json"
    json_file = open(json_path + "/" + chunk, "w")
    string_json = json.dumps(data, indent=4, separators=(',',':'))
    json_file.write("%s" % string_json)
    json_file.close()

    # insert json file to pgsql using ogr2ogr - folder_name = area name
    command = 'ogr2ogr -append -f "PostgreSQL" PG:"dbname=pgis host=' + dbHost + ' user=' + dbUsername + ' password=' + dbPassword + '" --config PG_USE_COPY YES -nln ' + folder_name + " "
    chunk_path = './mbtiles/' + folder_name + '/' + chunk
    res = os.system(command + ' ' + chunk_path)

    if res != 0:
        print "Error inserting into the database. This is most often due to running out of Memory (RAM), or incorrect database credentials... quitting"
        sys.exit()

    print "inserted chunk " + str(chunk_num) + " to db"

# ---------------------------------------------------------------------------------------
# START OF EXECUTABLE
# ---------------------------------------------------------------------------------------
# get name of h5 file and the groupname of that file's data
# ---------------------------------------------------------------------------------------
#  BEGIN EXECUTABLE
# ---------------------------------------------------------------------------------------
def build_parser():
    dbHost = "insarmaps.rsmas.miami.edu"
    parser = argparse.ArgumentParser(description='Convert a Unavco format H5 file for ingestion into insarmaps.')
    required = parser.add_argument_group("required arguments")
    required.add_argument("-f", "--file", help="unavco file to ingest", required=True)
    required.add_argument("-u", "--user", help="username for the insarmaps database", required=True)
    required.add_argument("-p", "--password", help="password for the insarmaps database", required=True)
    required.add_argument("--host", default=dbHost, help="postgres DB URL for insarmaps database", required=True)

    return parser

def main():
    parser = build_parser()
    parseArgs = parser.parse_args()

    file_name = parseArgs.file
    dbUsername = parseArgs.user
    dbPassword = parseArgs.password
    dbHost = parseArgs.host

    path_name = file_name[:len(file_name)-3]
    region_file_name = file_name[:len(file_name)-3] + '_region.txt'
# ---------------------------------------------------------------------------------------
# start clock to track how long conversion process takes
    start_time = time.clock()

# search for region file - if exist get first line which is region name
    region_file = None
    region = "null"
    project_name = "null"
    try:
        region_file = open(region_file_name, "r")
        region = region_file.readline()
        project_name = region_file.readline()
        region_file.close()
    except:
        pass

# use h5py to open specified group(s) in the h5 file 
# then read datasets from h5 file into memory for faster reading of data
# depending on UNAVCO format, the main key to access groups might be '/GEOCODE'
    file = h5py.File(file_name,  "r")
    group = file['GEOCODE/timeseries']  # assuming there is only one main key called 'GEOCODE'

# get attributes (stored at root) of UNAVCO timeseries file
    attributes = file['/'].attrs
    attributes_dictionary = attributes.items()

# in timeseries group, there are 25 datasets 
# need to get datasets with dates - strings that can be converted to integers
    dataset_keys = []
    for k in group.keys():
        if k.isdigit():
            dataset_keys.append(k)
    dataset_keys.sort()

# get the attributes for calculating latitude and longitude
    x_step = float(attributes["X_STEP"])
    y_step = float(attributes["Y_STEP"])
    x_first = float(attributes["X_FIRST"])
    y_first = float(attributes["Y_FIRST"])
    num_columns = int(attributes["WIDTH"])
    num_rows = int(attributes["FILE_LENGTH"])
    print "columns: %d" % num_columns
    print "rows: %d" % num_rows

# array that stores dates from dataset_keys that have been converted to decimal
    decimal_dates = []

# read datasets in the group into a dictionary of 2d arrays and intialize decimal dates
    timeseries_datasets = {}
    for key in dataset_keys:
        timeseries_datasets[key] = group[key][:]
        d = get_date(key)
        decimal = get_decimal_date(d)
        decimal_dates.append(decimal)

# set number of points per json chunk - then close h5 file
    chunk_size = 20000
    file.close()

# connect to postgresql database
# also create folder named after h5 file to store json files in mbtiles folder
    con = None
    cur = None

    path_list = path_name.split("/")
    mbtiles_path = os.getcwd() + "/mbtiles"
    folder_name = path_name.split("/")[len(path_list)-1]
    json_path = mbtiles_path + "/" + folder_name

    try: # create path for folder that stores all mbtiles
        os.mkdir(mbtiles_path)
    except:
        print mbtiles_path + " already exists"

    try: # create path for json
        os.mkdir(json_path)
    except:
        print json_path + " already exists"

# read and convert the datasets, then write them into json files and insert into database
    convert_data()

# run tippecanoe command to get mbtiles file and then delete the json files to save space
    os.chdir(os.path.abspath(json_path))
    os.system("tippecanoe *.json -x d -pf -pk -Bg -d9 -D12 -g12 -r0 -o " + folder_name + ".mbtiles")
    os.system("rm -rf *.json")

# move mbtiles file from json folder to mbtiles folder and then delete json folder
    os.system("mv " + folder_name + ".mbtiles " + os.path.abspath(mbtiles_path))
    os.system("rm -rf " + os.path.abspath(json_path))

# ---------------------------------------------------------------------------------------
# check how long it took to read h5 file data and create json files
    end_time =  time.clock()
    print ("time elapsed: " + str(end_time - start_time))
# ---------------------------------------------------------------------------------------

if __name__ == '__main__':
    main()
