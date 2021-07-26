#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Andre Theron, Zhang Yunjun, Jun 2019             #
# Email: andretheronsa@gmail.com                           #
############################################################


import os
import sys
import argparse
import datetime
import math
from mintpy.objects.sensor import standardize_sensor_name
from mintpy.utils import readfile, writefile, utils as ut


SPEED_OF_LIGHT = 299792458  # m / s


##################################################################################################
DESCRIPTION = """
  For each interferogram, coherence or unwrapped .dim product this script will prepare.rsc 
  metadata files for for mintpy based on .dim metadata file.

  The SNAP .dim file should contain all the required sensor / baseline metadata needed.
  The baseline metadata gets written during snap back-geocoding (co-registration).
  prep_snap is run seperately for unw/ifg/cor files so neeeds seperate .dim/.data products
  with only the relevant band in each product. Use Band Subset > save BEAM-DIMAP file.

  The file name should be yyyymmdd_yyyymmdd_type_tc.dim where type can be filt/unw/coh.

  The DEM should be prepared by adding an elevation file to a coregestered product - 
  then extract the elevation band only. Use Band Subset > save BEAM-DIMAP file

  Currently only works for geocoded (terrain correction step in SNAP) interferograms.
"""

EXAMPLE = """example:
  prep_snap.py  ../interferograms/*/*/Unw_*.img
  prep_snap.py  ../dem_tc.data/dem*.img
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Prepare attributes file for SNAP products.\n'+DESCRIPTION,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs='+', help='SNAP data file(s) in *.img format.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    inps.file = ut.get_file_list(inps.file, abspath=True)
    return inps


##################################################################################################
def get_ellpsoid_local_radius(xyz):
    """Calculate satellite height and ellpsoid local radius from orbital state vector
    Parameters: xyz          - tuple of 3 float, orbital state vector
    Reference:  height       - float, satellite altitude in m
                earth_radius - float, Earth radius in m
    """

    # Code simplified from from ISCE2.isceobj.Planet.xyz_to_llh())
    a = 6378137.000 # Semimajor of WGS84 - from ISCE2.isceobj.Planet.AstronomicalHandbook
    e2 = 0.0066943799901 # Eccentricity squared of WGS84 from ISCE2.isceobj.Planet.AstronomicalHandbook
    a2 = a**2
    e4 = e2**2
    r_llh = [0]*3
    d_llh = [[0]*3 for i in range(len(xyz))]
    xy2 = xyz[0]**2+xyz[1]**2
    p = (xy2)/a2
    q = (1.-e2)*xyz[2]**2/a2
    r = (p+q-e4)/6.
    s = e4*p*q/(4.*r**3)
    t = (1.+s+math.sqrt(s*(2.+s)))**(1./3.)
    u = r*(1.+t+1./t)
    v = math.sqrt(u**2+e4*q)
    w = e2*(u+v-q)/(2.*v)
    k = math.sqrt(u+v+w**2)-w
    D = k*math.sqrt(xy2)/(k+e2)
    r_llh[0] = math.atan2(xyz[2],D)
    r_llh[1] = math.atan2(xyz[1],xyz[0])
    r_llh[2] = (k+e2-1.)*math.sqrt(D**2+xyz[2]**2)/k
    d_llh[0] = math.degrees(r_llh[0])
    d_llh[1] = math.degrees(r_llh[1])
    d_llh[2] = r_llh[2]
    height = r_llh[2]

    # Calculate local radius
    # code from ISCE2.isceobj.Ellipsoid.localRadius()
    b = a * (1.0-e2)**0.5
    latg = math.atan(math.tan(math.radians(d_llh[0]))*a**2/b**2)
    arg = math.cos(latg)**2/a**2 + math.sin(latg)**2/b**2
    earth_radius = 1.0/math.sqrt(arg)

    return height, earth_radius


def extract_snap_metadata(fname):
    ''' read and extract attributes from a SNAP .dim file
    Parameters: fname - str, path of SNAP .dim file
    Returns:    atr   - dict, metadata
    '''

    # Read XML and extract required vals - using basic file reader 
    # xmltree or minidom might be better but this works    
    with open(fname, 'r') as f:
        lines = f.readlines()
        # use lists so that elements where there are more than one, only the first mention can be extracted - 
        # Usually the first mention (list item 0) is the final subsetted/geocoded product metadata
        bp, azimuth_looks, range_looks, az_pixel, rg_pixel, dates, x, y, z = [], [], [], [], [], [], [], [], []
        for line in lines:
            if "Perp Baseline" in line:
                bp.append(line.split(">")[1].split("<")[0])
            if "Master: " in line:
                dates.append(line.split(": ")[1].split('">')[0])
            if "radar_frequency" in line:
                freq = float(line.split(">")[1].split("<")[0]) * 1e6
                wvl = SPEED_OF_LIGHT / freq
            if "NROWS" in line:
                length = int(line.split(">")[1].split("<")[0])
            if "NCOLS" in line:
                width = int(line.split(">")[1].split("<")[0])
            if "DATA_TYPE" in line:
                data_type = line.split(">")[1].split("<")[0]
            if "IMAGE_TO_MODEL_TRANSFORM" in line:
                transform = str(line.split(">")[1].split("<")[0]).split(",")
            if "antenna_pointing" in line:
                looking = line.split(">")[1].split("<")[0]
                if looking == "right":
                    antenna_side = "-1"
            if "PRODUCT_SCENE_RASTER_START_TIME" in line:
                first_time = line.split(">")[1].split("<")[0]
            if "PRODUCT_SCENE_RASTER_STOP_TIME" in line:
                last_time = line.split(">")[1].split("<")[0]
            if "first_near_lat" in line:
                first_near_lat = line.split(">")[1].split("<")[0]
            if "first_near_long" in line:
                first_near_long = line.split(">")[1].split("<")[0]
            if "first_far_lat" in line:
                first_far_lat = line.split(">")[1].split("<")[0]
            if "first_far_long" in line:
                first_far_long = line.split(">")[1].split("<")[0]
            if "last_near_lat" in line:
                last_near_lat = line.split(">")[1].split("<")[0]
            if "last_near_long" in line:
                last_near_long = line.split(">")[1].split("<")[0]
            if "last_far_lat" in line:
                last_far_lat = line.split(">")[1].split("<")[0]
            if "last_far_long" in line:
                last_far_long = line.split(">")[1].split("<")[0]
            if "ASCENDING or DESCENDING" in line:
                direction = line.split(">")[1].split("<")[0]            
            if "azimuth_looks" in line:
                azimuth_looks.append(line.split(">")[1].split("<")[0])
            if "range_looks" in line:
                range_looks.append(line.split(">")[1].split("<")[0])
            if "pulse_repetition_frequency" in line:
                prf = line.split(">")[1].split("<")[0]
            if "slant_range_to_first_pixel" in line:
                starting_range = line.split(">")[1].split("<")[0]
            if "Satellite mission" in line:
                platform = line.split(">")[1].split("<")[0]
            if "platformHeading" in line:
                heading = line.split(">")[1].split("<")[0]
            if "Azimuth sample spacing" in line:
                az_pixel.append(line.split(">")[1].split("<")[0])
            if "Range sample spacing" in line:
                rg_pixel.append(line.split(">")[1].split("<")[0])
            if "incidenceAngleMidSwath" in line:
                inc_angle_mid = line.split(">")[1].split("<")[0]
            if "x_pos" in line:
                x.append(line.split(">")[1].split("<")[0])
            if "y_pos" in line:
                y.append(line.split(">")[1].split("<")[0])
            if "z_pos" in line:
                z.append(line.split(">")[1].split("<")[0])

    atr = {}
    # Calculate the center of the scene as floating point seconds
    start_utc = datetime.datetime.strptime(first_time, '%d-%b-%Y %H:%M:%S.%f')
    end_utc   = datetime.datetime.strptime(last_time,  '%d-%b-%Y %H:%M:%S.%f')
    center_utc = start_utc + ((end_utc - start_utc) / 2)
    center_seconds = (center_utc.hour * 3600.0 + 
                      center_utc.minute * 60. + 
                      center_utc.second)
    atr["CENTER_LINE_UTC"] = center_seconds

    # Extract reference / secondary date in yymmdd format
    date1 = datetime.datetime.strptime(dates[0], "%d%b%Y").strftime("%Y%m%d")[2:8]
    date2 = datetime.datetime.strptime(dates[1], "%d%b%Y").strftime("%Y%m%d")[2:8]
    atr["DATE12"] = "{}-{}".format(date1, date2)

    # calculate ellipsoid radius and satellite height from state vector
    # using first state vector for simplicity
    xyz = (float(x[0]), float(y[0]), float(z[0]))
    height, earth_radius = get_ellpsoid_local_radius(xyz)

    # Calculate range/azimuth pixel sizes
    range_pixel_size = float(rg_pixel[0]) * math.sin(float(inc_angle_mid))
    azimuth_pixel_size = float(az_pixel[0])*((earth_radius + height)/earth_radius)

    # Add values to dict
    atr['PROCESSOR'] = 'snap'
    atr['FILE_TYPE'] = os.path.basename(fname).split('_')[-2] # yyyymmdd_yyyymmdd_type_tc.dim
    atr["WIDTH"] = width
    atr["LENGTH"] = length
    atr["DATA_TYPE"] = data_type
    atr["WAVELENGTH"] = wvl
    atr["P_BASELINE_TOP_HDR"] = bp[1]
    atr["P_BASELINE_BOTTOM_HDR"] = bp[1]
    atr["ANTENNA_SIDE"] = antenna_side
    atr["LAT_REF1"], atr["LONG_REF1"] = first_near_lat, first_near_long
    atr["LAT_REF2"], atr["LONG_REF2"] = first_far_lat, first_far_long
    atr["LAT_REF3"], atr["LONG_REF3"] = last_near_lat, last_near_long
    atr["LAT_REF4"], atr["LONG_REF4"] = last_far_lat, last_far_long
    atr["ORBIT_DIRECTION"] = direction
    atr["ALOOKS"] = int(float(azimuth_looks[0]))
    atr["RLOOKS"] = int(float(range_looks[0]))
    atr["PRF"] = prf
    atr["PLATFORM"] = standardize_sensor_name(platform)
    atr["HEADING"] = heading
    atr["EARTH_RADIUS"] =  earth_radius
    atr["HEIGHT"] = height
    atr["RANGE_PIXEL_SIZE"] = range_pixel_size
    atr["AZIMUTH_PIXEL_SIZE"] = azimuth_pixel_size

    atr['INCIDENCE_ANGLE'] = inc_angle_mid
    atr["STARTING_RANGE"] = starting_range
    atr["SLANT_RANGE_DISTANCE"] = ut.incidence_angle2slant_range_distance(atr, inc_angle_mid)

    # Convert 3.333e-4 to 0.0003333
    transform = [str(float(i)) for i in transform]
    atr["X_STEP"], atr["Y_STEP"] = transform[0], transform[3]
    atr["X_FIRST"], atr["Y_FIRST"] = transform[4], transform[5]
    atr["X_UNIT"] = "degrees"
    atr["Y_UNIT"] = "degrees"

    # convert all key value in string format to ensure the --update checking in write_rsc()
    for key, value in atr.items():
        atr[key] = str(value)
    return atr


def write_rsc(atr, out_file):
    """Write to rsc file"""
    # grab atr from existing rsc file
    try:
        atr_orig = readfile.read_roipac_rsc(out_file)
    except:
        atr_orig = dict()

    # (over)write to rsc file if input atr has more items
    if not set(atr.items()).issubset(set(atr_orig.items())):
        atr_out = {**atr_orig, **atr}
        print('write metadata to {} '.format(os.path.basename(out_file)))
        writefile.write_roipac_rsc(atr_out, out_file=out_file)
    return out_file


##################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    for img_file in inps.file:
        if not img_file.endswith('.img'):
            raise ValueError('input data file does not end with .img: {}'.format(img_file))

        # *.dim metadata file for the input *.img data file
        # the map info from *.img.hdr file is NOT right, thus, not used.
        dim_file = os.path.dirname(img_file)[:-4]+'dim'

        # get metadata dict from *.dim file
        atr = extract_snap_metadata(dim_file)

        # write metadata dict to *.rsc file
        rsc_file = img_file + '.rsc'
        write_rsc(atr, out_file=rsc_file)

    return


##################################################################################################
if __name__ == "__main__":
    main(sys.argv[1:])
