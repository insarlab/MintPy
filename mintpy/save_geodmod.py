#!/usr/bin/env python3
############################################################
# Program is used for extract accumulated displacement of period
# Author: Lv Xiaoran                                       #
# Created: August 2019                                     #
############################################################

import os
import re
import time
import datetime
import shutil
import argparse
import subprocess
import numpy as np

import mintpy
import mintpy.workflow  #dynamic import for modules used by pysarApp workflow
from mintpy.objects import sensor, RAMP_LIST
from mintpy.utils import readfile, writefile, utils as ut
from mintpy.defaults.auto_path import autoPath

######################################################################################
EXAMPLE = """example:
  save_geodmod.py $SCRATCHDIR/DarbandikhanSenAT72/PYSAR/timeseries_ECMWF_demErr.h5 20171117_20180603 -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001 -t unw -o geo_20171117_20180603.unw -outdir $SCRATCHDIR/preGeodmod/ 
  save_geodmod.py $SCRATCHDIR/DarbandikhanSenAT72/PYSAR/INPUTS/ifgramStack.h5 20180603_20180615 -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001 -t cor -o geo_20180603_20180615.cor -outdir $SCRATCHDIR/preGeodmod/
  save_geodmod.py $SCRATCHDIR/DarbandikhanSenAT72/PYSAR/temporalCoherence.h5 20171117_20180603 -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001 -t cor -o geo_20171117_20180603.cor -outdir $SCRATCHDIR/preGeodmod/
  save_geodmod.py $SCRATCHDIR/DarbandikhanSenAT72/PYSAR/INPUTS/geometryRadar.h5 20171117_20180603 -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001 -t dem -o srtm.dem -outdir $SCRATCHDIR/preGeodmod/

"""

def create_parser():
    parser = argparse.ArgumentParser(description='Get horizontal and vertical displacements based on ascending and descending data',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs=1, help='ascending and descending timeseries files\n')
    parser.add_argument('dset', nargs=1, help='date12 of timeseries to be converted')
    parser.add_argument('-b', '--bbox', dest='SNWE', type=float, nargs=4, metavar=('S', 'N', 'W', 'E'),
                        help='Bounding box of area to be geocoded.\n' +
                        'Include the uppler left corner of the first pixel' +
                        '    and the lower right corner of the last pixel')
    parser.add_argument('-y', '--latstep', dest='latStep', type=float,
                        help='output pixel size in degree in latitude.')
    parser.add_argument('-x', '--lonstep', dest='lonStep', type=float,
                        help='output pixel size in degree in longitude.')
    parser.add_argument('-t','--type',dest='Type',type=str,
                        help='The type of data with choice [unw, cor, dem]')
    parser.add_argument('-o', '--output', dest='outfile', nargs=1,
                        help='output file name')
    parser.add_argument('-outdir','--outdir',dest='outdir',nargs=1,
                        help='output directory')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps
    
def format_args(arr, dst=""):
    """ parse list array(list item or values) to string """
    for k in arr:
        if type(k) == list:
            dst = dst + " " + format_args(k)
        else:
            dst = dst + " " + str(k)
    return dst.strip()
   
def get_path_com(path):
    """ return directory(filepath), filename(shotname), extension """
    (filepath, tempfilename) = os.path.split(os.path.abspath(path))
    (filename, extension) = os.path.splitext(tempfilename)
    return filepath, filename, extension

def write_rsc_file(inps,in_file,out_file):
    """ write rsc file for Geodmod just estract several properities from rsc file"""
    # read file
    meta = readfile.read_roipac_rsc(in_file)
    # initiate dict
    rsc = dict()
    rsc['FILE_DIR'] = "".join(inps.outdir)
    rsc['FILE_LENGTH'] = meta["FILE_LENGTH"]
    rsc['WIDTH'] = meta["WIDTH"]
    rsc['XMIN'] = 0
    rsc['XMAX'] = int(meta["WIDTH"]) - 1
    rsc['YMIN'] = 0
    rsc['YMAX'] = int(meta["FILE_LENGTH"]) - 1
    rsc['X_FIRST'] = float(meta["X_FIRST"])
    rsc['Y_FIRST'] = float(meta["Y_FIRST"])
    rsc['X_STEP'] = float(meta["X_STEP"])
    rsc['Y_STEP'] = float(meta["Y_STEP"])
    rsc['X_UNIT'] = 'degrees'
    rsc['Y_UNIT'] = 'degrees'
    rsc['RLOOKS'] = meta["RLOOKS"]
    rsc['ALOOKS'] = meta["ALOOKS"]
    rsc['Z_OFFSET'] = 0
    rsc['Z_SCALE'] = 1
    rsc['PROJECTION'] = 'LATLON'
    rsc['DATE12'] = '111111-222222'
    # write rsc file
    writefile.write_roipac_rsc(rsc, out_file, print_msg=True)
    return out_file 

def processdata(inps):
    #use geocode.py and save_roipac.py to process data"
    atr_asc = inps.file[0]
    os.chdir(os.path.dirname(atr_asc))
    
    if os.path.exists("".join(inps.outdir))=='False':
        os.mkdir("".join(inps.outdir))
    asc_args = [atr_asc, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", asc_args)
    args_str = format_args(asc_args)
    mintpy.geocode.main(args_str.split())

    os.chdir("".join(inps.outdir))
    filename, extension = get_path_com(atr_asc)[1:3]
    if filename=='ifgramStack':
        asct_args = ['geo_'+filename+extension, 'coherence-'+inps.dset[0], '-o', inps.outfile]
    elif filename=='temporalCoherence':     
        asct_args = ['geo_'+filename+extension, '-o', inps.outfile]    
    elif filename=='geometryRadar':
        asct_args = ['geo_'+filename+extension, 'height', '-o', inps.outfile]        
    else:
        asct_args = ['geo_'+filename+extension, inps.dset[0], '-o', inps.outfile]
    print("save_roipac.py", asct_args)
    asct_str = format_args(asct_args)
    mintpy.save_roipac.main(asct_str.split())

######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    
   
    if inps.Type=='unw':
        #  geocode timeseries**.h5 file and get the deformation field of two time periods
        processdata(inps)

    elif inps.Type=='cor':
        #  geocode ifgramStack.h5 file and get the coherence of two time periods
        processdata(inps)
    else:
        #  geocode geometryRadar.h5 file and get the dem
        
        processdata(inps)
        outfile = format_args([format_args(inps.outfile) + '.rsc'])
        #rscdir = os.path.dirname(os.path.abspath(outfile))
        #print(rscdir)        
        #os.chdir(rscdir)
        write_rsc_file(inps,outfile,format_args([format_args(inps.outfile) +'.rsc1']))
        os.remove(outfile)
        print('rename *.rsc1 to *.rsc')
        os.rename(format_args([format_args(inps.outfile) +'.rsc1']),outfile)
        
        

######################################################################################
if __name__ == '__main__':
    main()
