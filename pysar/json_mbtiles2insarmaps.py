#! /usr/bin/env python

import sys
import argparse
from pysar.add_attribute_insarmaps import InsarDatabaseController
import os
import pycurl
from cStringIO import StringIO
import urllib

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

def upload_json(folder_path):
    global dbUsername, dbPassword, dbHost

    table_name = get_file_name(folder_path)
    for json_chunk in os.listdir(folder_path):
        # insert json file to pgsql using ogr2ogr
        command = 'ogr2ogr -append -f "PostgreSQL" PG:"dbname=pgis host=' + dbHost + ' user=' + dbUsername + ' password=' + dbPassword + '" --config PG_USE_COPY YES -nln "' + table_name + '" ' + folder_path + '/' + json_chunk

        res = os.system(command)

        if res != 0:
            print "Error inserting into the database. This is most often due to running out of Memory (RAM), or incorrect database credentials... quitting"
            sys.exit()

        print "Inserted " + json_chunk + " to db"

    # create index
    print "Creating index on " + table_name
    attributesController = InsarDatabaseController(dbUsername, dbPassword, dbHost, 'pgis')
    attributesController.connect()
    attributesController.index_table_on(table_name, "p", None)
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
        attributesController = InsarDatabaseController(dbUsername, dbPassword, dbHost, 'pgis')
        table_name = get_file_name(folder_path)
        attributesController.connect()
        if attributesController.table_exists(table_name.lower()):
            print "Deleting old timeseries table"
            attributesController.remove_point_table_if_there(table_name.lower())
        attributesController.close()
 
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
