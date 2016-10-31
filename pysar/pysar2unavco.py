#! /usr/bin/env python

import h5py
import numpy as np
import datetime
from datetime import date
import sys
import os
from osgeo import ogr
import getopt
import re

# COMMAND I USE TO RUN SCRIPT: 
# python ~/code/insar_map_mvc/storage/json/pysar2unavco_back.py -t timeseries.h5 -i incidence_angle.h5 -d DEM_error.h5 -c temporal_coherence.h5 -m mask.h5 
# ---------------------------------------------------------------------------------------
# FUNCTIONS
# ---------------------------------------------------------------------------------------
# takes string formatted (YYYY-MM-DD) and returns date object
def get_date(date_string): 
	year = int(date_string[0:4])
	month = int(date_string[4:6])
	day = int(date_string[6:8])
	return date(year, month, day)

def mask_timeseries(timeseries_file, mask_file, out_file):
	os.system("masking.py -f " + timeseries_file + " -m " + mask_file + " -o " + out_file)
	
def usage():
	print 'Correct format: python pysar2unavco.py -t timeseries.h5 -i incidence_angle.h5 -d DEM_error.h5 -c temporal_coherence.h5 -m mask.h5'
	print 'Optional: --add_option OPTION_NAME=OPTION_VALUE, --mask_data (mask data using mask file)'

# ---------------------------------------------------------------------------------------
#  BEGIN EXECUTABLE
# ---------------------------------------------------------------------------------------
try:
	opts, extraArgs = getopt.getopt(sys.argv[1:],'t:i:d:c:m:', ['add_option=', 'mask_data']) 
except getopt.GetoptError:
	print 'Error while retrieving operations - exit'
	usage()
	sys.exit()

# read operations and arguments(file names):
operations = {}
added_options = {}
mask_data = False

for o, a in opts:
	if o == '--add_option':
		option_and_value = a.split('=')
		added_options[option_and_value[0].upper()] = option_and_value[1]
	elif o == "-t":
		timeseries = a
	elif o == "-i":
		incidence_angle = a
	elif o == "-d":
		dem = a
	elif o == "-c":
		temporal_coherence = a
	elif o == "-m":
		mask = a
	elif o == '--mask_data':
		mask_data = True
	else:
		assert False, "unhandled option - exit"
		sys.exit()

# ---------------------------------------------------------------------------------------
#  GET TIMESERIES FILE INTO UNAVCO and encode attributes to timeseries group in unavco file
# ---------------------------------------------------------------------------------------
if mask_data:
	print "Masking timeseries file named " + timeseries
	out_name = timeseries.split('.')[0] + "_masked.h5"
	mask_timeseries(timeseries, mask, out_name)
	timeseries = out_name 

print 'trying to open ' + timeseries
timeseries_file = h5py.File(timeseries, "r")
try: 
	key = 'timeseries'
	group = timeseries_file[key]
except: 
	print "failed using key '" + key + "' on " + dem
	sys.exit()

# get attributes and dates from timeseries file, sort dates to make life easier
attributes = group.attrs
dates = group.keys()
dates.sort()

# create unavco timeseries file - will change name of file to proper format later in script
unavco = "unavco_test.h5"
unavco_file = h5py.File(unavco, "w")
print "created test file"

# copy datasets from timeseries file into unavco file - main group key is 'GEOCODE'
# similar format key to the interferogram files found in unavco portal
unavco_file.create_group('/GEOCODE/timeseries')
for d in dates:
	old_key = '/timeseries/' + d
	timeseries_file.copy(old_key, unavco_file['/GEOCODE/timeseries'], d)

# ---------------------------------------------------------------------------------------
#  ENCODE REQURIED ATTRIBUTES FROM TIMESERIES FILE INTO UNAVCO
# ---------------------------------------------------------------------------------------
x_step = float(attributes["X_STEP"])	# longitude = width = x = columns
y_step = float(attributes["Y_STEP"])	# latitude = length = y = rows
x_first = float(attributes["X_FIRST"])
y_first = float(attributes["Y_FIRST"])
num_columns = int(attributes["WIDTH"])
num_rows = int(attributes["FILE_LENGTH"])

# extract attributes from project name
# KyushuT73F2980_2990AlosD - 73 = track number, 2980_2990 = frames, Alos = mission
try:
	project_name = attributes["PROJECT_NAME"]
except Exception, e:
	print "Project name is not in the h5 file, trying to find supplied name on the command line"
	key = "project_name"

	try:
		project_name = added_options[key.upper()]
	except Exception, e:
		print "Project name not supplied on the command line... quitting"
		sys.exit()

print project_name
track_index = project_name.find('T')
frame_index = project_name.find('F')
track_number = project_name[track_index+1:frame_index]
region_name = project_name[:track_index]

# sometimes there is only one frame number instead of framenumber_framenumber - look for "_"
multipleFrames = False
try:
	underscore = re.search("_", project_name).group(0)
	multipleFrames = True
except:
	pass

if multipleFrames:
	frames = re.search("\d+_\d+", project_name).group(0)
	first_frame = frames.split("_")[0]
	last_frame = frames.split("_")[1]
else:
	frames = re.search("\d+", project_name[frame_index+1:]).group(0)
	first_frame = frames
	last_frame = frames

mission_index = project_name.find(frames) + len(frames)
mission = project_name[mission_index:len(project_name)-1]

group = unavco_file['/']
# 1) mission = last chars of a folder name - SinabungT495F40_50AlosA -> mission = Alos
group.attrs['mission'] = mission

# 2) beam_mode = Yunjun said SM (stripmap) for now, later read from encoded timeseries.h5
group.attrs['beam_mode'] = 'SM'

# 3) beam_swath = Yunjun said SM (stripmap) for now
group.attrs['beam_swath'] = 'SM'

# 4) relative_orbit = number after T -> in case of Sinabung example above, 495
group.attrs['relative_orbit'] = int(track_number)

# 5 + 6) first_date and last date of timeseries datasets
group.attrs['first_date'] = get_date(dates[0]).isoformat()
group.attrs['last_date'] = get_date(dates[len(dates)-1]).isoformat()
print "type of last_date: " + str(type(get_date(dates[len(dates)-1]).isoformat()))

# 7) scene footprint = wkt polygon
# Create polygon square with boundaries of datasets
ring = ogr.Geometry(ogr.wkbLinearRing) # points added as (x,y) or (long,lat)
ring.AddPoint( x_first, y_first )
ring.AddPoint( x_first, (y_step * num_rows + y_first) )
ring.AddPoint( (x_first + x_step * num_columns),  y_first)
ring.AddPoint( (x_first + x_step * num_columns), (y_step * num_rows + y_first) )
geom_poly = ogr.Geometry(ogr.wkbPolygon)
geom_poly.AddGeometry(ring)
wkt = geom_poly.ExportToWkt()
group.attrs['scene_footprint'] = wkt

# 8) processing type = LOS_TIMESERIES
group.attrs['processing_type'] = 'LOS_TIMESERIES'

# 9) processing software = PySAR or StaMPS - PySAR for now
group.attrs['processing_software'] = 'PySAR'

# 10) history = (not important) the day I processed h5 files into UNAVCO format
group.attrs['history'] = datetime.datetime.now().date().isoformat()

# ---------------------------------------------------------------------------------------
#  ENCODE RECOMMENDED ATTRIBUTES FROM TIMESERIES FILE INTO UNAVCO
# ---------------------------------------------------------------------------------------
# UNAVCO wants this to be an int but we have multiple frames so we have two frame attributes
group.attrs['first_frame'] = int(first_frame)
group.attrs['last_frame'] = int(last_frame)

# flight_direction = A or D (ascending or descending)
# tried to encode as char but python seems to only know string
group.attrs['flight_direction'] = chr(ord(attributes['ORBIT_DIRECTION'][0].upper()))

# for now, all falks files are Right (R)
# tried to encode as char but python seems to only know string
group.attrs['look_direction'] = chr(ord('R'))

# polarization = naN (Yunjun doesn't care for this attribute)
# group.attrs['polarization'] = np.nan

# attribute in yunjun/falk's h5 files: prf, wavelength
group.attrs['prf'] = float(attributes['PRF'])
group.attrs['wavelength'] = float(attributes['WAVELENGTH'])

# atmos_correct_method = 'ECMWF' harcode for now, later read from timeseries.h5 after Yunjun encodes attribute
group.attrs['atmos_correct_method'] = 'ECMWF'

# processing_dem = yunjun is adding attribute to geo_timeseries.h5, wait for now

# unwrap method = SNAPHU for all files right now
group.attrs['unwrap_method'] = 'SNAPHU'

# post_processing_method = PySAR or StaMPS - PySAR for now
group.attrs['post_processing_method'] = 'PySAR'

# master_platform, percent_atmos, and other interferogram attributes - dont need

# baseline_perp = attribute called P_BASELINE_TIMESERIES
# UNAVCO format pdf says: This really only applies to interferograms, but for 
# time series or velocity products this could be the min and max bperp.
group.attrs['min_baseline_perp'] = float(attributes['P_BASELINE_BOTTOM_HDR'])
group.attrs['max_baseline_perp'] = float(attributes['P_BASELINE_TOP_HDR'])

# ---------------------------------------------------------------------------------------
#  ENCODE ATTRIBUTES FROM FALK'S FILES TO CALCULATE LATITUDE AND LONGITUDE - NOT IN UNAVCO DOC
# ---------------------------------------------------------------------------------------
group.attrs['X_STEP'] = x_step
group.attrs['Y_STEP'] = y_step
group.attrs['X_FIRST'] = x_first
group.attrs['Y_FIRST'] = y_first
group.attrs['WIDTH'] = num_columns
group.attrs['FILE_LENGTH'] = num_rows

timeseries_file.close()
# ---------------------------------------------------------------------------------------
#  GET INCIDENCE_ANGLE FILE INTO UNAVCO
# ---------------------------------------------------------------------------------------
# Yunjun encoded incidence_angle.h5 with key 'mask' but it should be 'incidence_angle'
# try both keys
print 'trying to open ' + incidence_angle
incidence_angle_file = h5py.File(incidence_angle, "r")
is_open = True
try:
	key = 'incidence_angle'
	group = incidence_angle_file[key]
except:
	print "failed using key '" + key + "' on " + incidence_angle
	is_open = False

if not is_open:
	try:
		key = 'mask'
		group = incidence_angle_file[key]
	except:
		print "failed using key '" + key + "' on " + incidence_angle
		print incidence_angle + ' could not be read'
		sys.exit()

# in current format incidence_angle.h5 requires you to use key 'incidence_angle/incidence_angle'
# in order to get the datasets from the h5 file - perhaps in later version we fix this to make it
# less verbose
for k in group.keys():
	old_key = key + "/" + k
	incidence_angle_file.copy(old_key, unavco_file['/GEOCODE'], 'incidence_angle')

incidence_angle_file.close()

# ---------------------------------------------------------------------------------------
#  GET DEM FILE INTO UNAVCO
# ---------------------------------------------------------------------------------------
dem_file = h5py.File(dem, "r")
print 'trying to open ' + dem
try:
	key = 'dem'
	group = dem_file[key]
except:
	print "failed using key '" + key + "' on " + dem

for k in group.keys():
	old_key = key + "/" + k
	dem_file.copy(old_key, unavco_file['/GEOCODE'], key)

dem_file.close()

# ---------------------------------------------------------------------------------------
#  GET TEMPORAL COHERENCE FILE INTO UNAVCO
# ---------------------------------------------------------------------------------------
temporal_coherence_file = h5py.File(temporal_coherence, "r")
print 'trying to open ' + temporal_coherence
try:
	key = 'temporal_coherence'
	group = temporal_coherence_file[key]
except:
	print "failed using key '" + key + "' on " + temporal_coherence

for k in group.keys():
	old_key = key + "/" + k
	temporal_coherence_file.copy(old_key, unavco_file['/GEOCODE'], key)

temporal_coherence_file.close()

# ---------------------------------------------------------------------------------------
#  GET MASK FILE INTO UNAVCO
# ---------------------------------------------------------------------------------------
mask_file = h5py.File(mask, "r")
print 'trying to open ' + mask
try:
	key = 'mask'
	group = mask_file[key]
except:
	print "failed using key '" + key + "' on " + mask

for k in group.keys():
	old_key = key + "/" + k
	mask_file.copy(old_key, unavco_file['/GEOCODE'], key)

mask_file.close()

# CLOSE UNAVCO FILE
unavco_file.close()

# IMPORTANT: RENAME UNAVCO file to proper format based on file attributes
# example - pysar file is called Kyushu T 80 F 245_246 JersD.h5
# UNAVCO timeseries file is called JERS_SM_80_245_246_<first date>_<last date>.h5 since we dont need TBASE or BPERP for timeseries
unavco_name = mission + '_SM_' + track_number + '_' + frames + '_' + dates[0] + '_' + dates[len(dates)-1] + '.h5'
os.rename(unavco, "./" + unavco_name)

# create a text file to store region (ex: Kyushu) needed for database but not unavco format
attributes_file_name = unavco_name[:len(unavco_name)-3] + '_region.txt'
attributes_file = open(attributes_file_name, "w")
attributes_file.write(region_name)
attributes_file.close()

# ---------------------------------------------------------------------------------------
#  TEST THE UNAVCO FILE FOR DATASETS
# ---------------------------------------------------------------------------------------
# open unavco h5 file that was created earlier for testing
f = h5py.File(unavco_name, "r")
print "opened "  + unavco_name

print 'trying to print the attribute types'
attributes = f['/'].attrs.items()
for k, v in attributes:
	print "key: " + k + ", value: " + str(v) + ", type: " + str(type(v)) + "\n"

print 'trying to print main group GEOCODE'
geocode = f["GEOCODE"]
print geocode

print 'trying to print group timeseries'
group = f['GEOCODE/timeseries']
print group

print 'trying to print last dataset of timeseries'
dataset = group[dates[len(dates)-1]]
print dataset

print 'trying to print dataset incidence_angle'
dataset = f['GEOCODE/incidence_angle']
print dataset

print 'trying to print dataset dem'
dataset = f['GEOCODE/dem']
print dataset

print 'trying to print dataset temporal_coherence'
dataset = f['GEOCODE/temporal_coherence']
print dataset

print 'trying to print dataset mask'
dataset = f['GEOCODE/mask']
print dataset

f.close()
