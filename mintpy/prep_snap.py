#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright(c) 2019, Andre Theron, Zhang Yunjun            #
# Author: Andre Theron, Zhang Yunjun, Jun 2019             #
# Email: andretheronsa@gmail.com                           #
############################################################


import os
import argparse
import datetime
import math
from mintpy.objects.sensor import standardedSensorNames
from mintpy.utils import readfile, writefile, utils as ut


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

  Metadata extraction from .dim was built around products generated via the following
  workflow:
    - Read
    - Slice assembly (if required)
    - Apply orbit file
    - Split product (extract relevant polarisation and subswath swaths)
    + Following is done per subswath [IW1, IW2, IW3]
        - Back geocoding
        - Enhanced spectral diversity (if more than one burst)
        - Interferogram generation
        - TOPSAR Deburst
        - Topophase removal
        - Goldstein phase filtering
    - Merge subswaths (if more than one swath was done)
    - Add elevation
    - Coherence
    - Subset geometry and all seperate bands (elevation, ifg, coh) (by geometry only possible after working with bursts)
    - Snaphu export for ifg
    - SNAPHU phase unwrapping
    - Snaphu import
    - Terrain correct (all products)
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
    Parameters: xyz : tuple of 3 float, orbital state vector
    Reference: isce2.isceobj.Planet
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


# Parse .dim file
def utm_dim_to_rsc(fname):
    ''' read and extract attributes from a SNAP .dim file and create mintpy .rsc file
    Inputs: 
        fname: SNAP .dim interferogram/coherence/unwrapped filepath
    Outputs:
        atr : dict, Attributes dictionary
    '''

    # Read XML and extract required vals - using basic file reader 
    # xmltree or minidom might be better but this works    
    with open(fname) as f:
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
                rf = 299.792458 / float(line.split(">")[1].split("<")[0])  
            if "NROWS" in line:
                lenght = int(line.split(">")[1].split("<")[0])
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
            if '"x"' in line:
                x.append(line.split(">")[1].split("<")[0])
            if '"y"' in line:
                y.append(line.split(">")[1].split("<")[0])
            if '"z"' in line:
                z.append(line.split(">")[1].split("<")[0])
            
    # Do datetime things needed after the loop
    # Calculate the center of the scene as floating point seconds
    first_time_dt = datetime.datetime.strptime(first_time, '%d-%b-%Y %H:%M:%S.%f')
    last_time_dt = datetime.datetime.strptime(last_time, '%d-%b-%Y %H:%M:%S.%f')
    center_time_dt = first_time_dt + ((last_time_dt - first_time_dt) / 2)
    center_time_s = float(datetime.timedelta(
        hours=center_time_dt.hour, 
        minutes=center_time_dt.minute,
        seconds=center_time_dt.second,
        microseconds=center_time_dt.microsecond).total_seconds())
    # Extract master slave date in yymmdd format
    master_dt = datetime.datetime.strptime(dates[0], "%d%b%Y")
    master = master_dt.strftime("%Y%m%d")[2:8]
    slave_dt = datetime.datetime.strptime(dates[1], "%d%b%Y")
    slave = slave_dt.strftime("%Y%m%d")[2:8]

    xyz = (float(x[0]), float(y[0]), float(z[0])) # Using first state vector for simplicity
    height, earth_radius = get_ellpsoid_local_radius(xyz)

    # Calculate range/azimuth pixel sizes
    range_pixel_size = float(rg_pixel[0]) * math.sin(float(inc_angle_mid))
    azimuth_pixel_size = float(az_pixel[0])*((earth_radius + height)/earth_radius)

    # Add values to dict
    atr = {}
    atr['PROCESSOR'] = 'snap'
    atr['FILE_TYPE'] = os.path.basename(fname).split('_')[-2] # yyyymmdd_yyyymmdd_type_tc.dim
    atr["WAVELENGTH"] = rf
    atr["P_BASELINE_TOP_HDR"], atr["P_BASELINE_BOTTOM_HDR"] = bp[1], bp[1]
    atr["DATE12"] = str(master) + "-" + str(slave)
    atr["WIDTH"] = width
    atr["LENGTH"], atr["FILE_LENGTH"]  = lenght, lenght
    atr["DATA_TYPE"] = data_type
    atr["ANTENNA_SIDE"] = antenna_side
    atr["CENTER_LINE_UTC"] = center_time_s
    atr["LAT_REF1"], atr["LONG_REF1"] = first_near_lat, first_near_long
    atr["LAT_REF2"], atr["LONG_REF2"] = first_far_lat, first_far_long
    atr["LAT_REF3"], atr["LONG_REF3"] = last_near_lat, last_near_long
    atr["LAT_REF4"], atr["LONG_REF4"] = last_far_lat, last_far_long
    atr["ORBIT_DIRECTION"] = direction
    atr["ALOOKS"], atr["RLOOKS"] = int(float(azimuth_looks[0])), int(float(range_looks[0]))
    atr["PRF"] = prf
    atr["STARTING_RANGE"] = starting_range
    atr["PLATFORM"] = standardedSensorNames[platform.replace('-','').lower()]
    atr["HEADING"] = heading
    atr["EARTH_RADIUS"] =  earth_radius
    atr["HEIGHT"] = height
    atr["RANGE_PIXEL_SIZE"] = range_pixel_size
    atr["AZIMUTH_PIXEL_SIZE"] = azimuth_pixel_size

    # Convert 3.333e-4 to 0.0003333
    transform = [str(float(i)) for i in transform]
    atr["X_STEP"], atr["Y_STEP"] = transform[0], transform[3]
    atr["X_FIRST"], atr["Y_FIRST"] = transform[4], transform[5]

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
        # *.dim metadata file for the input *.img data file
        dim_file = os.path.dirname(img_file)[:-4]+'dim'

        # get metadata dict from *.dim file
        atr = utm_dim_to_rsc(dim_file)

        # write metadata dict to *.rsc file
        rsc_file = os.path.splitext(img_file)[0]+'.rsc'
        rsc_file = write_rsc(atr, out_file=rsc_file)
    return


##################################################################################################
if __name__ == "__main__":
    main()
