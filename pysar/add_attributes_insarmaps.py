#! /usr/bin/env python

import psycopg2
import sys
import getopt
import os
import argparse

class InsarDatabaseController:
    def __init__(self, username, password, host, db):
        self.username = username
        self.password = password
        self.host = host
        self.db = db
        self.con = None
        self.cursor = None

    def connect(self):
        try:
            self.con = psycopg2.connect("dbname='pgis' user='" + self.username + "' host='" + self.host + "' password='" + self.password + "'")
            self.cursor = self.con.cursor()
        except Exception, e:
            print "Error While Connecting"
            print e
            sys.exit()

    def close(self):
        self.con.close()
        self.con = None
        self.cursor = None

    def get_dataset_names(self):
        sql = "SELECT * FROM area"
        self.cursor.execute(sql)

        return self.cursor.fetchall()

    def get_dataset_id(self, dataset):
        sql = "SELECT id from area WHERE area.project_name = '" + dataset + "'"
        self.cursor.execute(sql)

        return self.cursor.fetchone()[0]

    def table_exists(self, table):
        sql = "SELECT exists(SELECT * FROM information_schema.tables WHERE table_name=%s)"
        self.cursor.execute(sql, (table,))

        return self.cursor.fetchone()[0]

    def attribute_exists_for_dataset(self, dataset, attributekey):
        dataset_id = self.get_dataset_id(dataset)

        sql = "SELECT exists(SELECT attributekey FROM extra_attributes WHERE area_id = " + str(dataset_id) + " AND attributekey = '" + attributekey + "');"
        self.cursor.execute(sql)

        return self.cursor.fetchone()[0]

    def add_attribute(self, dataset, attributekey, attributevalue):
        dataset_id = self.get_dataset_id(dataset)
        sql = ""
        prepared_values = None

        if not self.table_exists("extra_attributes"):
            sql = "CREATE TABLE IF NOT EXISTS extra_attributes (area_id integer, attributekey varchar, attributevalue varchar);"
            self.cursor.execute(sql)
            self.con.commit()

        if not self.attribute_exists_for_dataset(dataset, attributekey):
            sql = "INSERT INTO extra_attributes VALUES (%s, %s, %s);"
            prepared_values = (str(dataset_id), attributekey, attributevalue)
        else:
            sql = "UPDATE extra_attributes SET attributevalue = %s WHERE area_id = %s AND attributekey = %s"
            prepared_values = (attributevalue, str(dataset_id), attributekey)

        self.cursor.execute(sql, prepared_values)
        self.con.commit()

    def index_table_on(self, table, on, index_name):
        # can't remove single quotes from table name, so we do it manually
        sql = None
        if index_name:
            sql = "CREATE INDEX " + index_name + " ON " + table + " (" + on + ");"
        else:
            sql = "CREATE INDEX ON " + table + " (" + on + ");"

        try:
            self.cursor.execute(sql)
            self.con.commit()
        # index exists most probably if exception thrown
        except Exception, e:
            pass

    def remove_dataset(self, unavco_name, project_name):
        # try to drop table first in case extra_attributes or area table isn't populated
        sql = "DROP TABLE " + unavco_name.lower()
        self.cursor.execute(sql)
        self.con.commit()

        # then try to delete from area and extra_attributes
        dataset_id = self.get_dataset_id(project_name)
        sql = "DELETE from area WHERE id = " + str(dataset_id)
        self.cursor.execute(sql)
        sql = "DELETE from extra_attributes WHERE area_id = " + str(dataset_id)
        self.cursor.execute(sql)
        self.con.commit()
            

def usage():
    print "add_atributes.py -u USERNAME -p PASSWORD -h HOST -d DB -f FILE"

def parse_file_for_attributes(file):
    attributes_dict = {}

    with open(file) as f:
        file_contents = f.readlines()

        for line in file_contents:
            if line != '\n':
                key_value = line.split("=") # index 0 is key, 1 is value

                if key_value[1][-1] == '\n':
                    key_value[1] = key_value[1][:-1]
                    key_value[1] = key_value[1].lstrip() # take out leading whitespace

                attributes_dict[key_value[0]] = key_value[1]

    return attributes_dict

def build_parser():
    dbHost = "insarmaps.rsmas.miami.edu"
    parser = argparse.ArgumentParser(description='Edit attributes of an insarmaps dataset')
    required = parser.add_argument_group("required arguments")
    required.add_argument("-f", "--folder", help="folder of the dataset to look for add_Attribute.txt", required=True)
    required.add_argument("-u", "--user", help="username for the insarmaps database", required=True)
    required.add_argument("-p", "--password", help="password for the insarmaps database", required=True)
    required.add_argument("--host", default=dbHost, help="postgres DB URL for insarmaps database", required=True)
    required.add_argument("-d", "--db", help="postgres database", required=True)

    return parser

def main(argv):
    parser = build_parser()
    parseArgs = parser.parse_args()

    username = parseArgs.user
    password = parseArgs.password
    host = parseArgs.host
    db = parseArgs.db
    working_dir = parseArgs.folder


    # make sure we have a final / so the below code doesn't break
    if working_dir[-1] != "/":
        working_dir += "/"

    project_name = working_dir.split("/")[-2]
    attributes_file = working_dir + "add_Attribute.txt"
    attributes = parse_file_for_attributes(attributes_file)
    dbController = InsarDatabaseController(username, password, host, db)    
    dbController.connect()

    for key in attributes:
        print "Setting attribute " + key + " to " + attributes[key]
        dbController.add_attribute(project_name, key, attributes[key])

    dbController.index_table_on("extra_attributes", "area_id", "area_id_idx")
    dbController.close()

if __name__ == '__main__':
    main(sys.argv)
