#! /usr/bin/env python

import sys
import os
import getopt

# extract project name from command line supplied path
def project_name_from_path(path):
	tokens = path.split("/")
	last_token = tokens[len(tokens) - 1]

	if last_token != '':
		return last_token

	return tokens[len(tokens) - 2]

# get folders in path sorted in ascending time order
def sorted_ls(path):
    mtime = lambda f: os.stat(os.path.join(path, f)).st_mtime
    return list(sorted(os.listdir(path), key=mtime))

# same as sorted_ls but reversed
def rev_sorted_ls(path):
        files = sorted_ls(path)
        files.reverse()

        return files

# figure out what the name of the h5 file to put on site is in the given path
def get_H5_filename(path):
        files = rev_sorted_ls(path)
	region_file = None

	for file in files:
		if "_region.txt" in file:
			region_file = file
			break

	h5_file = region_file.split("_region.txt")[0] + ".h5"
        
	return h5_file # h5 file to go on site should be second newest file

def usage():
	print "pysar2insarmapsBatch.py -f path/dir -u DB_USERNAME -p DBPASSWORD -h DB_PASSWORD"

dbUsername = "INSERT"
dbPassword = "INSERT"
dbHost = "INSERT"

path = None

try:
	opts, extraArgs = getopt.getopt(sys.argv[1:],'f:u:p:h:')
except getopt.GetoptError:
	print 'Error while retrieving operations - exit'
	usage()
	sys.exit()

for o, a in opts:
	if o == '-f':
		path = a
	elif o == '-u':
		dbUsername = a
	elif o == '-p':
		dbPassword = a
	elif o == '-h':
		dbHost = a
	else:
		assert False, "unhandled option - exit"
		sys.exit()

bjob_script_filename = "run_pysar2insarmaps.py"
path_absolute = os.path.abspath(path)

h5_file = get_H5_filename(path)
h5_file_partial_name = h5_file.split(".")[0]

cur_proj_name = project_name_from_path(path)

# create working directory in scratch and copy relevant files over
scratch_dir = os.environ["SCRATCHDIR"] + "/" + cur_proj_name
print "making directory " + scratch_dir
os.system("mkdir " + scratch_dir)
command = "cp " + cur_proj_name + "/" + h5_file_partial_name + "*" + " " + scratch_dir
print "copying files to scratch with command " + command
os.system(command)

# go to scratch dir, and run the bjob command
command = "echo unavco2insarmaps.py -f " + h5_file + " -u " + dbUsername + " -p " + dbPassword + " -h " + dbHost + " > " + bjob_script_filename

mbtiles_filename = h5_file.split(".")[0] + ".mbtiles"
os.chdir(scratch_dir)
os.system(command)

os.system("createBatch.pl " + bjob_script_filename)
command = "cp -r mbtiles/" + h5_file.split(".")[0] + ".mbtiles " + path_absolute
print "bjob finished, trying to execute " + command
os.system(command)

