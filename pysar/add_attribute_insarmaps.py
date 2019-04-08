#!/usr/bin/env python3

import psycopg2
import sys
import argparse
import pysar.utils.readfile as readfile
import pycurl
from io import BytesIO
import urllib.request, urllib.parse, urllib.error
import requests
from requests.packages.urllib3.exceptions import InsecureRequestWarning

requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

# TODO: fix these classes. apparantly, need to call commit() method after execute even if you fetch something

class InsarDatabaseController(object):
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
        except Exception as e:
            print("Error While Connecting")
            print(e)
            sys.exit()

    def close(self):
        self.con.close()
        self.con = None
        self.cursor = None

    def run_raw_query(self, query):
        self.cursor.execute(query)
        self.con.commit()

        if self.cursor.rowcount > 0:
            return self.cursor.fetchall()

        return None

    def get_dataset_names(self):
        sql = "SELECT * FROM area"
        self.cursor.execute(sql)

        return self.cursor.fetchall()

    def get_dataset_id(self, dataset):
        sql = "SELECT id from area WHERE area.unavco_name = '" + dataset + "'"
        self.cursor.execute(sql)
        id = self.cursor.fetchone()

        if id:
            return id[0]

        return -1

    def table_exists(self, table):
        sql = "SELECT exists(SELECT * FROM information_schema.tables WHERE table_name=%s)"
        self.cursor.execute(sql, (table,))

        return self.cursor.fetchone()[0]

    # TODO refactor below two functions
    def attribute_exists_for_dataset(self, dataset, attributekey):
        dataset_id = self.get_dataset_id(dataset)

        sql = "SELECT exists(SELECT attributekey FROM extra_attributes WHERE area_id = " + str(dataset_id) + " AND attributekey = '" + attributekey + "');"
        self.cursor.execute(sql)

        return self.cursor.fetchone()[0]
    
    def plot_attribute_exists_for_dataset(self, dataset, attributekey):
        dataset_id = self.get_dataset_id(dataset)

        sql = "SELECT exists(SELECT attributekey FROM plot_attributes WHERE area_id = " + str(dataset_id) + " AND attributekey = '" + attributekey + "');"
        self.cursor.execute(sql)

        return self.cursor.fetchone()[0]

    def add_attribute(self, dataset, attributekey, attributevalue):
        dataset_id = self.get_dataset_id(dataset)
        sql = ""
        prepared_values = None

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

    def add_plot_attribute(self, dataset, attributekey, plotAttributeJSON):
        dataset_id = self.get_dataset_id(dataset)
        sql = ""
        prepared_values = None

        sql = "CREATE TABLE IF NOT EXISTS plot_attributes (area_id integer, attributekey varchar, attributevalue json);"
        self.cursor.execute(sql)
        self.con.commit()

        if not self.plot_attribute_exists_for_dataset(dataset, attributekey):
            sql = "INSERT INTO plot_attributes VALUES (%s, %s, %s);"
            prepared_values = (str(dataset_id), attributekey, plotAttributeJSON)
        else:
            sql = "UPDATE plot_attributes SET attributevalue = %s WHERE area_id = %s AND attributekey = %s"
            prepared_values = (plotAttributeJSON, str(dataset_id), attributekey)

        self.cursor.execute(sql, prepared_values)
        self.con.commit()
        
    def index_table_on(self, table, on, index_name):
        # can't remove single quotes from table name, so we do it manually
        sql = None
        if index_name:
            sql = 'CREATE INDEX "' + index_name + '" ON "' + table + '" (' + on + ');'
        else:
            sql = 'CREATE INDEX ON "' + table + '" (' + on + ');'

        try:
            self.cursor.execute(sql)
            self.con.commit()
        # index exists most probably if exception thrown
        except Exception as e:
            print(str(e))

    def cluster_table_using(self, table, index_name):
        sql = None
        sql = 'CLUSTER "' + table + '" USING "' + index_name + '";'

        try:
            self.cursor.execute(sql)
            self.con.commit()
        # index exists most probably if exception thrown
        except Exception as e:
            print(str(e))

    def remove_point_table_if_there(self, table_name):
        sql = 'DROP TABLE IF EXISTS "' + table_name + '"'
        self.cursor.execute(sql)
        self.con.commit()
        
    def create_area_table_if_not_exists(self):
        # create area table if not exist - limit for number of dates is 2    00, limt for number of attribute keys/values is 100
        self.cursor.execute("CREATE TABLE IF NOT EXISTS area ( unavco_name varchar, project_name varchar, longitude double precision, latitude double precision, country varchar, region varchar, numchunks integer, attributekeys varchar[100], attributevalues varchar[100], stringdates varchar[200], decimaldates double precision[200] );")
        self.con.commit()

    def insert_dataset_into_area_table(self, area, project_name,
    mid_long, mid_lat, country, region, chunk_num, attribute_keys,
    attribute_values, string_dates_sql, decimal_dates_sql):
        # put dataset into area table
        query = "INSERT INTO area VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"

        preparedValues = (area, project_name,
    mid_long, mid_lat, country, region, chunk_num, attribute_keys,
    attribute_values, string_dates_sql, decimal_dates_sql)
        self.cursor.execute(query, preparedValues)
        self.con.commit()

    def remove_dataset_if_there(self, unavco_name):
        dataset_id = self.get_dataset_id(unavco_name)

        if dataset_id == -1:
            return

        dataset_id_str = str(dataset_id)
        table_name = dataset_id_str
        self.remove_point_table_if_there(table_name)
        # then try to delete from area and extra_attributes
        try:
            dataset_id = self.get_dataset_id(unavco_name)
            sql = "DELETE from area WHERE id = " + dataset_id_str
            self.cursor.execute(sql)
            sql = "DELETE from extra_attributes WHERE area_id = " + dataset_id_str
            self.cursor.execute(sql)
            sql = "DELETE from plot_attributes WHERE area_id = " + dataset_id_str
            self.cursor.execute(sql)
            self.con.commit()
        except Exception:
            pass

class InsarDatasetController(InsarDatabaseController):
    def __init__(self, username, password, host, db, serverUsername, serverPassword):
        super(InsarDatasetController, self).__init__(username, password, host, db)
        self.bodyOutput = BytesIO()
        self.headersOutput = BytesIO()
        self.serverUsername = serverUsername
        self.serverPassword = serverPassword

    def setup_curl(self):
        curl = pycurl.Curl()
        curl.setopt(curl.WRITEFUNCTION, self.bodyOutput.write)
        curl.setopt(curl.HEADERFUNCTION, self.headersOutput.write)
        curl.setopt(pycurl.COOKIEFILE, "")
        curl.setopt(pycurl.SSL_VERIFYPEER, 0)   
        curl.setopt(pycurl.SSL_VERIFYHOST, 0)

        # hackish way of finding if the host url gets redirected to https version of site
        # TODO: find more elegant solution when time permits, and use requests instead of curl
        # for all of these HTTP requests...
        self.host = requests.get("https://" + self.host, verify=False).url

        return curl

    def curl_login(self, username, password):
        curl = self.setup_curl()
        curl.setopt(curl.POST, 1)
        loginParams =  urllib.parse.urlencode([("email", username), ("password", password)])
        curl.setopt(curl.POSTFIELDS, loginParams)
        loginURL = self.host + "/auth/login"
        curl.setopt(curl.URL, loginURL)
        curl.perform()

        return curl

    def upload_mbtiles(self, fileName):
        curl = self.curl_login(self.serverUsername, self.serverPassword)

        curl.setopt(curl.HTTPPOST, [('title', fileName), (('file', (curl.FORM_FILE, fileName)))])
        uploadURL = self.host + "/WebServices/uploadMbtiles"
        curl.setopt(curl.URL, uploadURL)
        #curl.setopt(curl.VERBOSE, 1)
        curl.perform()

        responseCode = curl.getinfo(pycurl.HTTP_CODE)
        if responseCode == 200:
            print("Successfully uploaded " + fileName)
        elif responseCode == 302:
            sys.stderr.write("Server redirected us... Please check username and password, and try again")
        elif responseCode == 301:
            sys.stderr.write("Server redirected us... Please check host address, and try again")
        else:
            sys.stderr.write("The server responded with code: " + str(responseCode))

    def remove_mbtiles(self, fileName):
        curl = self.curl_login(self.serverUsername, self.serverPassword)

        curl.setopt(curl.HTTPPOST, [('fileName', fileName)])
        deleteURL = self.host + "/WebServices/deleteMbtiles"
        curl.setopt(curl.URL, deleteURL)
        #curl.setopt(curl.VERBOSE, 1)
        curl.perform()

        responseCode = curl.getinfo(pycurl.HTTP_CODE)

        if responseCode == 302:
            sys.stderr.write("Server redirected us... Please check username and password, and try again")
        elif responseCode == 301:
            sys.stderr.write("Server redirected us... Please check host address, and try again")
            
def build_parser():
    dbHost = "insarmaps.rsmas.miami.edu"
    parser = argparse.ArgumentParser(description='Edit attributes of an insarmaps dataset')
    required = parser.add_argument_group("required arguments")
    required.add_argument("-f", "--folder", help="folder of the dataset to look for add_Attribute.txt", required=True)
    required.add_argument("-u", "--user", help="username for the insarmaps database", required=True)
    required.add_argument("-p", "--password", help="password for the insarmaps database", required=True)
    required.add_argument("--host", default=dbHost, help="postgres DB URL for insarmaps database", required=True)
    required.add_argument("-d", "--db", help="postgres database", required=True)
    required.add_argument("-U", "--unavco_name", help="UNAVCO name of this dataset", required=True)

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

    unavco_name = parseArgs.unavco_name
    attributes_file = working_dir + "add_Attribute.txt"
    attributes = readfile.read_template(attributes_file)
    dbController = InsarDatabaseController(username, password, host, db)
    dbController.connect()

    for key in attributes:
        print("Setting attribute " + key + " to " + attributes[key])
        if key == "plotAttributes":
            dbController.add_plot_attribute(unavco_name, key, attributes[key])
        else:
            dbController.add_attribute(unavco_name, key, attributes[key])

    dbController.index_table_on("extra_attributes", "area_id", "area_id_idx")
    dbController.close()

if __name__ == '__main__':
    main(sys.argv)
