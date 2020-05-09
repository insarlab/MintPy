#! /usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Alfredo Terreco, 2016                            #
############################################################


import os
import sys
import argparse
import pickle
import numpy

from mintpy.insarmaps.add_attribute_insarmaps import (
    InsarDatabaseController,
    InsarDatasetController,
)


dbUsername = "INSERT"
dbPassword = "INSERT"
dbHost = "INSERT"

# use pickle file to get unavco_name of dataset
def get_unavco_name(json_path):
    insarmapsMetadata = None
    fileName = json_path + "/metadata.pickle"

    with open(fileName, "rb") as file:
        insarmapsMetadata = pickle.load(file)

    return insarmapsMetadata["area"]

def upload_insarmaps_metadata(fileName):
    insarmapsMetadata = None

    with open(fileName, "rb") as file:
        insarmapsMetadata = pickle.load(file)

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
        # convert numpy.int64 objects to native python types otherwise psycopg2 can't upload to db
        # needed because we use pickle.HIGHEST_PROTOCOL to serialize with pickle now
        if isinstance(v, numpy.int64):
            v = v.item()

        if k in needed_attributes:
            attributesController.add_attribute(area, k, v)
        elif k == "plotAttributes":
            attributesController.add_plot_attribute(area, k, v)

    attributesController.close()

def upload_json(folder_path):
    global dbUsername, dbPassword, dbHost
    attributesController = InsarDatabaseController(dbUsername, dbPassword,     dbHost, 'pgis')
    attributesController.connect()
    print("Clearing old dataset, if it is there")
    area_name = get_unavco_name(folder_path)
    attributesController.remove_dataset_if_there(area_name)
    attributesController.close()

    # uploading metadata for area. this creates entry into area table.
    # we need this entry to get the db to generate an id for area which
    # we use to name the corresponding table for the dataset
    upload_insarmaps_metadata(folder_path + "/metadata.pickle")
    # create index
    print("Creating index on " + area_name)
    attributesController.connect()
    area_id = str(attributesController.get_dataset_id(area_name))
    attributesController.close()
    firstJsonFile = True

    for file in os.listdir(folder_path):
        # insert json file to pgsql using ogr2ogr
        file_extension = file.split(".")[1]
        if file != "metadata.pickle" and file_extension != "mbtiles":
            command = 'ogr2ogr -append -f "PostgreSQL" PG:"dbname=pgis ' + \
                      ' host=' + dbHost + ' user=' + dbUsername + ' password=' + dbPassword + \
                      '" --config PG_USE_COPY YES -nln ' + area_id + ' ' + folder_path + '/' + file
            # only provide layer creation options if this is the first file
            if firstJsonFile:
                command = 'ogr2ogr -lco LAUNDER=NO -append -f "PostgreSQL" PG:"dbname=pgis ' + \
                          ' host=' + dbHost + ' user=' + dbUsername + ' password=' + dbPassword + \
                          '" --config PG_USE_COPY YES -nln ' + area_id + ' ' + folder_path + '/' + file
                firstJsonFile = False

            res = os.system(command)

            if res != 0:
                msg = "Error inserting into the database."
                msg += " This is most often due to running out of Memory (RAM)"
                msg += ", or incorrect database credentials... quitting"
                sys.stderr.write(msg)
                sys.exit()

            print("Inserted " + file + " to db")

    attributesController.connect()
    attributesController.index_table_on(area_id, "p", None)
    attributesController.cluster_table_using(area_id, area_id + "_p_idx")
    attributesController.close()

def build_parser():
    dbHost = "insarmaps.rsmas.miami.edu"
    parser = argparse.ArgumentParser(description='Convert a Unavco format HDF5 file for ingestion into insarmaps.')
    parser.add_argument("--json_folder", help="folder containing json to upload.", required=False)
    parser.add_argument("json_folder_positional", help="folder containing json to upload.", nargs="?")
    parser.add_argument("-U", "--server_user", required=False, 
        help="username for the insarmaps server (the machine where the tileserver and http server reside)")
    parser.add_argument("--remove",
        help="UNAVCO name of dataset to remove from insarmaps website", required=False)
    parser.add_argument("-P", "--server_password", required=False,
        help="password for the insarmaps server (the machine where the tileserver and http server reside)")
    parser.add_argument("--mbtiles_file", help="mbtiles file to upload", required=False)
    parser.add_argument("mbtiles_file_positional",
        help="mbtiles file to upload, as a positional argument", nargs="?")

    required = parser.add_argument_group("required arguments")
    required.add_argument("-u", "--user", help="username for the insarmaps database", required=True)
    required.add_argument("-p", "--password", help="password for the insarmaps database", required=True)
    required.add_argument("--host", default=dbHost, help="postgres DB URL for insarmaps database", required=True)

    return parser

def main():
    global dbUsername, dbPassword, dbHost
    parser = build_parser()
    parseArgs = parser.parse_args()
    dbUsername = parseArgs.user
    dbPassword = parseArgs.password
    dbHost = parseArgs.host

    if parseArgs.json_folder:
        print("Uploading json chunks...")
        upload_json(parseArgs.json_folder)
    elif parseArgs.json_folder_positional:
        print("Uploading json chunks....")
        upload_json(parseArgs.json_folder_positional)

    if parseArgs.mbtiles_file or parseArgs.mbtiles_file_positional:
        dbContoller = InsarDatasetController(dbUsername,
                                             dbPassword,
                                             dbHost,
                                             'pgis',
                                             parseArgs.server_user,
                                             parseArgs.server_password)

        if not parseArgs.server_user or not parseArgs.server_password:
            sys.stderr.write("Error: credentials for the insarmaps server not provided")
        elif parseArgs.mbtiles_file:
            print("Uploading mbtiles...")
            dbContoller.upload_mbtiles(parseArgs.mbtiles_file)
        else:
            print("Uploading mbtiles....")
            dbContoller.upload_mbtiles(parseArgs.mbtiles_file_positional)

    if parseArgs.remove:
        if not parseArgs.server_user or not parseArgs.server_password:
            sys.stderr.write("Error: credentials for the insarmaps server not provided")
        else:
            print("Removing " + parseArgs.remove)
            dbContoller = InsarDatasetController(dbUsername,
                                                 dbPassword,
                                                 dbHost,
                                                 'pgis',
                                                 parseArgs.server_user,
                                                 parseArgs.server_password)

            dbContoller.connect()
            dbContoller.remove_dataset_if_there(parseArgs.remove)
            dbContoller.close()
            dbContoller.remove_mbtiles(parseArgs.remove + ".mbtiles")

if __name__ == '__main__':
    main()
