#! /usr/bin/env python

import sys
import argparse
from pysar.add_attribute_insarmaps import InsarDatabaseController
import os
import pycurl
from cStringIO import StringIO
import urllib
import cPickle

dbUsername = "INSERT"
dbPassword = "INSERT"
dbHost = "INSERT"

def build_parser():
    dbHost = "insarmaps.rsmas.miami.edu"
    parser = argparse.ArgumentParser(description='Convert a Unavco format     H5 file for ingestion into insarmaps.')
    parser.add_argument("-f", "--folder", help="folder containing json to upload. The folder name will be used as the table name in the db to upload, so it should be as provided by unavco2json_mbtiles.py", required=False)
    parser.add_argument("-U", "--server_user", help="username for the insarmaps server (the machine where the tileserver and http server reside)", required=False)
    parser.add_argument("-P", "--server_password", help="password for the insarmaps server (the machine where the tileserver and http server reside)", required=False)
    parser.add_argument("-m", "--mbtiles", help="mbtiles file to upload", required=False)
    required = parser.add_argument_group("required arguments")
    required.add_argument("-u", "--user", help="username for the insarmaps database", required=True)
    required.add_argument("-p", "--password", help="password for the insarmaps database", required=True)
    required.add_argument("--host", default=dbHost, help="postgres DB URL     for insarmaps database", required=True)

    return parser

def get_file_name(fullPath):
    pathComponents = fullPath.split("/")
    for name in reversed(pathComponents):
        if name != "":
            return name

    return None

def upload_insarmaps_metadata(fileName):
    insarmapsMetadata = None

    with open(fileName, "r") as file:
        insarmapsMetadata = cPickle.load(file)

    area = insarmapsMetadata["area"]
    project_name = insarmapsMetadata["project_name"]
    mid_long = insarmapsMetadata["mid_long"]
    mid_lat = insarmapsMetadata["mid_lat"]
    country = insarmapsMetadata["country"]
    region = insarmapsMetadata["region"]
    chunk_num = insarmapsMetadata["chunk_num"]
    attribute_keys = insarmapsMetadata["attribute_keys"]
    attribute_values = insarmapsMetadata["attribute_values"]
    string_dates_sql = insarmapsMetadata["string_dates_sql"]
    decimal_dates_sql = insarmapsMetadata["decimal_dates_sql"]
    attributes = insarmapsMetadata["attributes"]
    needed_attributes = insarmapsMetadata["needed_attributes"]

    attributesController = InsarDatabaseController(dbUsername, dbPassword, dbHost, 'pgis')
    attributesController.connect()
    attributesController.create_area_table_if_not_exists()
    attributesController.insert_dataset_into_area_table(area, project_name,
    mid_long, mid_lat, country, region, chunk_num, attribute_keys,
    attribute_values, string_dates_sql, decimal_dates_sql)

    # put in attributes into standalone attributes table
    for k in attributes:
        v = attributes[k]
        if k in needed_attributes:
            attributesController.add_attribute(area, k, v)
        elif k == "plotAttributes":
            attributesController.add_plot_attribute(area, k, v)

    attributesController.close()

def upload_json(folder_path):
    global dbUsername, dbPassword, dbHost
    attributesController = InsarDatabaseController(dbUsername, dbPassword,     dbHost, 'pgis')
    attributesController.connect()
    print "Clearing old dataset, if it is there"
    area_name = get_file_name(folder_path)
    attributesController.remove_dataset_if_there(area_name)
    attributesController.close()

    area_name = get_file_name(folder_path)
    for json_chunk in os.listdir(folder_path):
        # insert json file to pgsql using ogr2ogr
        if json_chunk != "metadata.pickle":
            command = 'ogr2ogr -lco LAUNDER=NO -append -f "PostgreSQL" PG:"dbname=pgis host=' + dbHost + ' user=' + dbUsername + ' password=' + dbPassword + '" --config PG_USE_COPY YES -nln "' + area_name + '" ' + folder_path + '/' + json_chunk

            res = os.system(command)

            if res != 0:
                print "Error inserting into the database. This is most often due to running out of Memory (RAM), or incorrect database credentials... quitting"
                sys.exit()

            print "Inserted " + json_chunk + " to db"
        else:
            print "PICKLE!!!!!"

    # uploading metadata for area 
    upload_insarmaps_metadata(folder_path + "/metadata.pickle")
    # create index
    print "Creating index on " + area_name
    attributesController = InsarDatabaseController(dbUsername, dbPassword, dbHost, 'pgis')
    attributesController.connect()
    attributesController.index_table_on(area_name, "p", None)
    attributesController.close()

def upload_mbtiles(fileName, username, password):
    curl = pycurl.Curl()
    curl.setopt(curl.POST, 1)
    loginParams =  urllib.urlencode([("email", username), ("password", password)])
    curl.setopt(curl.POSTFIELDS, loginParams)
    loginURL = dbHost + "/auth/login"
    curl.setopt(curl.URL, loginURL)
    bodyOutput = StringIO()
    headersOutput = StringIO()
    curl.setopt(curl.WRITEFUNCTION, bodyOutput.write)
    curl.setopt(curl.HEADERFUNCTION, headersOutput.write)
    curl.setopt(pycurl.COOKIEFILE, "")
    curl.perform()
    curl.setopt(curl.HTTPPOST, [('title', fileName), (('file', (curl.FORM_FILE, fileName)))])
    uploadURL = dbHost + "/WebServices/uploadMbtiles"
    curl.setopt(curl.URL, uploadURL)
    #curl.setopt(curl.VERBOSE, 1)
    curl.perform()
    
    responseCode = curl.getinfo(pycurl.HTTP_CODE)
    if responseCode == 200:
        print "Successfully uploaded " + fileName
    elif responseCode == 302:
        print "Server redirected us... Please check username and password, and try again"
    else:
        print "The server responded with code: " + str(responseCode)

def main():
    global dbUsername, dbPassword, dbHost
    parser = build_parser()
    parseArgs = parser.parse_args()
    dbUsername = parseArgs.user
    dbPassword = parseArgs.password
    dbHost = parseArgs.host
    folder_path = parseArgs.folder

    if parseArgs.folder:
        print "Uploading json chunks..."
        upload_json(parseArgs.folder)

    if parseArgs.mbtiles:
        if not parseArgs.server_user or not parseArgs.server_password:
            print "Error: credentials for the insarmaps server not provided"
        else:
            print "Uploading mbtiles..."
            upload_mbtiles(parseArgs.mbtiles, parseArgs.server_user, parseArgs.server_password)

if __name__ == '__main__':
    main()
