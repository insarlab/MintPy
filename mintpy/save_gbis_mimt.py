#!/usr/bin/env python3
#################################################################
# Program is used for extract accumulated displacement of period#
# Author: Lv Xiaoran                                            #
# Created: August 2019                                          #
#################################################################

import os
import argparse
import string
import shutil
import numpy as np
import re
import scipy.io as sio

import mintpy
import mintpy.workflow  #dynamic import for modules used by pysarApp workflow
from mintpy.objects import sensor, timeseries
from mintpy.utils import ptime, readfile, writefile,utils as ut
from mimtpy.utils import  multitrack_utilities
######################################################################################
EXAMPLE = """example:
  for singletrack:
  Note: startDate, endDate and outdir have default values, any of them can be not given:
    startDate default value = atr['START_DATE']
    endDate default value = atr['END_DATE']
    outdir default value = '$MODELDIR/project/gbis_startDate_endDate/'
  save_gbis_mimt.py timeseries_ECMWF_demErr.h5 -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001 -m maskTempCoh.h5
  save_gbis_mimt.py ifgramStack.h5 -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001 -s 20171117 -e 20171129 -outdir $MODELDIR/DarbandikhanSenAT73
  save_gbis_mimt.py velocity.h5 -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001 -s 20171117 -e 20171129 -outdir $MODELDIR/DarbandikhanSenAT73
  save_gbis_mimt.py S1_IW23_026_0108_0113_20171117_XXXXXXXX.he5 -s 20171128 -e 20181210 
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Prepare data for GBIS software',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs='?', help='ascending or descending files\n')

    parser.add_argument('-b', '--bbox', dest='SNWE', type=float, nargs=4, metavar=('S', 'N', 'W', 'E'),
                        help='Bounding box of area to be geocoded.\n' +
                        'Include the uppler left corner of the first pixel' +
                        '    and the lower right corner of the last pixel')
    parser.add_argument('-y', '--latstep', dest='latStep', type=float,
                        help='output pixel size in degree in latitude.')
    parser.add_argument('-x', '--lonstep', dest='lonStep', type=float,
                        help='output pixel size in degree in longitude.')
    parser.add_argument('-s','--startDate',dest='startDate',nargs='?',
                        help='date1 of timeseires to be converted.The default is the StartDate')
    parser.add_argument('-e','--endDate',dest='endDate',nargs='?',
                        help='date2 of timeseries to be converted.The default is the EndDate')
    parser.add_argument('-m', '--mask', dest='mask_file', help='mask file.')
    parser.add_argument('--ref-lalo', dest='ref_lalo', type=float, nargs='*',
                        help='custom reference pixel in lat/lon')
    
    
    parser.add_argument('-outdir','--outdir',dest='outdir',nargs='?',default=os.getenv('PWD'),
                        help='output directory')

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    # default startDate and endDate
    if not os.path.isfile("".join(inps.file)):
        file = multitrack_utilities.find_timeseries(os.getcwd())
    else:
        file = "".join(inps.file)    
    atr = readfile.read_attribute(file)
    if not inps.startDate or inps.startDate=='None':
        inps.startDate = atr['START_DATE']
    if not inps.endDate or inps.endDate=='None':
        inps.endDate = atr['END_DATE']
    # set output dir
    inps.outdir =  multitrack_utilities.set_outdir(inps,'gbis')
    return inps    

def process_geocode(inps):
    """process temporalCoherence.h5 and geometryRadar.h5 file"""
    # process cor and dem dataset
    if os.path.exists("".join(inps.outdir))==False:
        os.makedirs("".join(inps.outdir))        
    
    demname = 'geometryRadar.h5'
    if not os.path.isfile(demname):
        demname_f = './inputs/geometryRadar.h5'
    else:
        demname_f = 'geometryRadar.h5'
    cmd_args = [demname_f, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", cmd_args)
    args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    mintpy.geocode.main(args_str.split())
    
    if inps.mask_file:
        cmd_args = [inps.mask_file, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
        print("geocode.py", cmd_args)
        args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
        mintpy.geocode.main(args_str.split())

def process_timeseries(inps):
    """geocode timeseries**.h5 file and get the deformation field of two time periods"""
    atr_asc = inps.file
   
    #unw file
    cmd_args = [atr_asc, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", cmd_args)
    args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    mintpy.geocode.main(args_str.split())
        
    #save dataset of unw
    os.chdir("".join(inps.outdir))
    filename, extension = multitrack_utilities.seprate_filename_extension("".join(atr_asc))[1:3]
    
    cmd_args = ['geo_'+filename+extension, "".join([inps.startDate,'_',inps.endDate])]
    print("save_roipac.py", cmd_args)
    asct_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    os.system(multitrack_utilities.seperate_str_byspace(['save_roipac.py', asct_str.split()]))


def process_ifgramStack(inps):
    """process ifgramStack.h5 file"""
    if os.path.exists("".join(inps.outdir))==False:
        os.makedirs("".join(inps.outdir))    

    #ifgramStack file
    atr_asc = ['./inputs/'+inps.file]
    cmd_args = [atr_asc, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", cmd_args)
    args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    mintpy.geocode.main(args_str.split())
        
    #save dataset of unw
    os.chdir("".join(inps.outdir))
    filename, extension = multitrack_utilities.seprate_filename_extension("".join(atr_asc))[1:3]
    
    cmd_args = ['geo_'+filename+extension, "".join(['unwrapPhase-',inps.startDate,'_',inps.endDate])]
    print("save_roipac.py", cmd_args)
    asct_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    os.system(multitrack_utilities.seperate_str_byspace(['save_roipac.py', asct_str.split()]))   

def process_velocity(inps):
    """process velocity.h5 file"""
    atr_asc = inps.file
   
    #velocity file
    cmd_args = [atr_asc, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", cmd_args)
    args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    mintpy.geocode.main(args_str.split())
    
    os.chdir("".join(inps.outdir))
    print('save unw file')
    multitrack_utilities.velo_disp(inps)

def process_HDFEOS(inps):
    """process S1*.h5 file"""

    atr_asc = inps.file
        
    #save dataset of unw and geometry
    filename, extension = multitrack_utilities.seprate_filename_extension("".join(atr_asc))[1:3]
    cmd_args = [filename+extension, "".join(['displacement-',inps.startDate,'_',inps.endDate]), '-o', "".join(['geo_',inps.startDate,'_',inps.endDate,'.unw'])]
    print("save_roipac.py", cmd_args)
    asct_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    os.system(multitrack_utilities.seperate_str_byspace(['save_roipac.py', asct_str.split()]))   
    
     #mv 
    if os.path.exists(inps.outdir)==False:
        os.makedirs(inps.outdir)
    else:
        shutil.rmtree(inps.outdir)
        os.makedirs(inps.outdir)
    key1 = 'geo_'
    for file in os.listdir(os.getcwd()):
        if str.find(file,key1) != -1:
            shutil.move(file,inps.outdir)        

def prep_gbis(inps):
    """prepare data that has to be written in *.mat file"""

    if str.find(str(inps.file),'S1')!= -1:       
        key1 = 'S1'
        key2 = '.he5'
        for file in os.listdir(os.getcwd()):
            if str.find(file,key1) != -1 and str.find(file,key2) != -1:
                shutil.copy(file,inps.outdir) 
        os.chdir("".join(inps.outdir))
        geom_file = multitrack_utilities.find_HDFEOS_fullname(os.getcwd())
    else:

        geom_file = 'geo_geometryRadar.h5'
    # metadata
    unw_file = 'geo_'+ inps.startDate + '_' + inps.endDate +'.unw'
    inps.metadata = readfile.read_attribute(unw_file)    
    inps.phase,atr = readfile.read(unw_file)
    
    # mask
    if not inps.mask_file or inps.mask_file == 'None':
        inps.mask = np.ones((int(inps.metadata['LENGTH']),
                             int(inps.metadata['WIDTH'])), dtype=np.bool_)
    else:
        inps.mask = readfile.read('geo_'+inps.mask_file)[0]
        
    # update mask to exclude pixel with NaN value
    inps.mask *= ~np.isnan(inps.phase)
    # set all masked out pixel to NaN
    inps.phase[inps.mask==0] = np.nan

    # change reference point
    if inps.ref_lalo:
        coord = ut.coordinate(inps.metadata)
        ref_lat, ref_lon = inps.ref_lalo
        ref_y, ref_x = coord.geo2radar(ref_lat, ref_lon)[0:2]
        # update data
        inps.phase -= inps.phase[ref_y, ref_x]
        # update metadata
        inps.metadata['REF_LAT'] = ref_lat
        inps.metadata['REF_LON'] = ref_lon
        inps.metadata['REF_Y'] = ref_y
        inps.metadata['REF_X'] = ref_x

    # read geometry
    inps.lat, inps.lon = ut.get_lat_lon(inps.metadata)
    inps.inc_angle = readfile.read(geom_file, datasetName='incidenceAngle')[0]
    inps.head_angle = np.ones(inps.inc_angle.shape, dtype=np.float32) * float(inps.metadata['HEADING'])
    inps.height = readfile.read(geom_file, datasetName='height')[0]
    inps.lat[inps.mask==0] = np.nan
    inps.lon[inps.mask==0] = np.nan
    inps.inc_angle[inps.mask==0] = np.nan
    inps.head_angle[inps.mask==0] = np.nan
    inps.height[inps.mask==0] = np.nan

    # output filename
    proj_name = atr['PROJECT_NAME']
    if not proj_name:
        raise ValueError('No custom/auto output filename found.')
    inps.outfile = '{}_{}_{}.mat'.format(proj_name, inps.startDate, inps.endDate)
    inps.outfile = os.path.join(inps.outdir, inps.outfile)
    inps.outfile = os.path.abspath(inps.outfile)
    
    #delete geo_*.h5 files
    multitrack_utilities.delete_tmpgeo(inps.outdir, 'S1_', '.he5')
    multitrack_utilities.delete_tmpgeo(inps.outdir, 'geo_', '.h5')
    return
    
def save2mat(inps):
    """write mat file"""
    mdict = {}
    # required by GBIS
    mdict['Heading'] = inps.head_angle[inps.mask].reshape(-1,1)
    mdict['Inc'] = inps.inc_angle[inps.mask].reshape(-1,1)
    mdict['Lat'] = inps.lat[inps.mask].reshape(-1,1)
    mdict['Lon'] = inps.lon[inps.mask].reshape(-1,1)
    mdict['Phase'] = inps.phase[inps.mask].reshape(-1,1)
    # optional
    mdict['Height'] = inps.height[inps.mask].reshape(-1,1)
    mdict['Mask'] = inps.mask
    mdict['Metadata'] = inps.metadata
    # save to mat file
    sio.savemat(inps.outfile, mdict, long_field_names=True)
    print('save to file: {}.mat'.format(os.path.abspath(inps.outfile)))
    return    
    
######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    print('single track!')
    inps.startDate,inps.endDate =  multitrack_utilities.find_start_end_date(os.getcwd(),inps)
    print(inps)
    
    if str.find(inps.file,'ifgramStack') != -1:
        process_geocode(inps)
        process_ifgramStack(inps)
    elif str.find(inps.file,'velocity.') != -1:
        process_geocode(inps)
        process_velocity(inps)
    elif str.find(inps.file,'timeseries') != -1:
        process_geocode(inps)
        process_timeseries(inps)
    else:
        process_HDFEOS(inps)
    prep_gbis(inps)
    save2mat(inps)


    
######################################################################################
if __name__ == '__main__':
    main()
