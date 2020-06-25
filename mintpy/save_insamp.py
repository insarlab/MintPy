#!/usr/bin/env python3
#################################################################
# Program is used for preparing data for InSAMP software        #
# Author: Lv Xiaoran                                            #
# Created: June 2020                                            #
#################################################################

import os
import argparse
import shutil
import numpy as np
import pyproj
import scipy.io as sio

import mintpy
from mintpy.utils import readfile, writefile, utils as ut
######################################################################################
EXAMPLE = """example:
  
  save_insamp.py geo_20171117_20200205.unw -g ./inputs/geometryRadar.h5 -o geo_20171117_20200205_insamp
  
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Prepare data for insamp software based on ROIPAC format',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs=1, type=str, help='geocoded unw or h5 files to be converted\n')

    parser.add_argument('-g', '--geometryRadar', dest='geometry', type=str, nargs=1,
                        help='geometry file')
    parser.add_argument('-o', '--outfile', dest='outfile', nargs=1, type=str,
                        help='outfile name')
    
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

def process_roi_pac(inps):
    """prepare data for kite using roi_pac format and calculate look anlge for 4 corners from incidence angle"""
    # geometry data 
    geometryRadar = inps.geometry[0]
    print('processing geometry data {}'.format(geometryRadar))
    inc_angle = readfile.read(geometryRadar, datasetName='incidenceAngle')[0]
    azi_angle = readfile.read(geometryRadar, datasetName='azimuthAngle')[0]

    # heading angle
    head_angle = ut.azimuth2heading_angle(azi_angle)
    
    # displacement data
    disp = inps.file[0]
    print('processing displacement data {}'.format(disp))
    disp_data, atr = readfile.read(disp)
    
    # write displacement data as InSAMP required. Because we don't have the amplitude data we use phase as amplitude 
    outfile_unw = inps.outfile[0] + '.unw'
    dsDict = dict()
    dsDict['amplitude'] = disp_data
    dsDict['phase'] = disp_data
    writefile.write(datasetDict=dsDict, out_file=outfile_unw, metadata=atr)   
    
    # write angle data as InSAMP required.
    outfile_los = inps.outfile[0] + '.los' 
    los_angle= np.vstack((inc_angle, head_angle)).flatten()
    los_angle = np.array(los_angle, dtype=np.float32)
    los_angle.tofile(outfile_los)

    return
######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    
    # generate dem.jpeg
    process_roi_pac(inps)
######################################################################################
if __name__ == '__main__':
    main()
