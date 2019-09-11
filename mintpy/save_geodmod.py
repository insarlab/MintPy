#!/usr/bin/env python3
#################################################################
# Program is used for extract accumulated displacement of period#
# Author: Lv Xiaoran                                            #
# Created: August 2019                                          #
#################################################################

import os
import sys
import argparse
import string
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib.colors import LightSource
import shutil
import numpy as np
import re
from PIL import Image

import mintpy
import mintpy.workflow  #dynamic import for modules used by pysarApp workflow
from mintpy.objects import sensor
from mintpy.utils import ptime, readfile, writefile,utils as ut
from mintpy.objects import timeseries
from MimtPy.utils import multitrack_utilities
######################################################################################
EXAMPLE = """example:
  for singletrack:
  Note: startDate, endDate and outdir have default values, any of them can be not given:
    startDate default value = atr['START_DATE']
    endDate default value = atr['END_DATE']
    outdir default value = '$PWD/geodmod_startDate_endDate/'
  save_geodmod.py timeseries_ECMWF_demErr.h5 -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001  
  save_geodmod.py ifgramStack.h5  -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001 -s 20171117 -e 20171129 -outdir $MODELDIR/Darbandikhan/SenAT73
  save_geodmod.py velocity.h5  -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001 -s 20171117 -e 20171129 -outdir $MODELDIR/Darbandikhan/SenAT73
  save_geodmod.py S1_IW23_026_0108_0113_20171117_XXXXXXXX.he5 -s 20171128 -e 20181210 
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Prepare data for Geodmod software',
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
    
    parser.add_argument('-outdir','--outdir',dest='outdir',nargs='?',default=os.getenv('PWD'),
                        help='output directory')

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
  
    # default startDate and endDate
    if not os.path.isfile("".join(inps.file)):
        file = multitrack_utilities.find_timeseries(os.getcwd()+'/')
    else:
        file = "".join(inps.file)  
    atr = readfile.read_attribute(file)
    if not inps.startDate or inps.startDate=='None':
        inps.startDate = atr['START_DATE']
    if not inps.endDate or inps.endDate=='None':
        inps.endDate = atr['END_DATE']
    # set output dir
    inps.outdir =  multitrack_utilities.set_outdir(inps,'geodmod')
    
    return inps
    
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

def dem_jpeg(dem_file):
    """generate dem.jepg file based on Yunjun's code"""
    out_file = dem_file+'.jpeg'
    rsc_file = out_file+'.rsc'
    shutil.copy2(dem_file+'.rsc', rsc_file)
    # read data
    dem = readfile.read(dem_file)[0]
    print('dem.shape:',dem.shape)
    # figure size
    ds_shape = tuple(reversed(dem.shape))
    fig_dpi = 300
    fig_size = [i / fig_dpi for i in ds_shape]
    print('fig_size:',fig_size)
    # color range
    disp_min = np.nanmin(dem) - 4000
    disp_max = np.nanmax(dem) + 2000
    # prepare shaded relief
    ls = LightSource(azdeg=315, altdeg=45)
    dem_shade = ls.shade(dem, vert_exag=0.3, cmap=plt.get_cmap('gray'), vmin=disp_min, vmax=disp_max)
    dem_shade[np.isnan(dem_shade[:, :, 0])] = np.nan
    print('dem_shade.shape:',dem_shade.shape)
    # plot
    fig, ax = plt.subplots(figsize=fig_size)
    ax.imshow(dem_shade, interpolation='spline16', origin='upper')
    # get rid of whitespace on the side
    ax.axis('off')
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
    # output
    print('save figure to file {}'.format(out_file))
    plt.savefig(out_file, transparent=True, dpi=300, pad_inches=0.0)
    
    #resize to desired size(FA 8/19, unclear why size is wrong)
    im = Image.open(out_file)
    im_out = im.resize(dem.shape, Image.NEAREST)
    im_out.save(out_file)
    
    plt.show()

def process_geocode(inps):
    """process temporalCoherence.h5 and geometryRadar.h5 file"""
    # process cor and dem dataset
    if os.path.exists("".join(inps.outdir))==False:
        os.makedirs("".join(inps.outdir))
        
    # geocode
    corname = 'temporalCoherence.h5'
    cmd_args = [corname, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", cmd_args)
    args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    mintpy.geocode.main(args_str.split())
    
    demname = 'geometryRadar.h5'
    if not os.path.isfile(demname):
        demname_f = './inputs/geometryRadar.h5'
    else:
        demname_f = 'geometryRadar.h5'
    cmd_args = [demname_f, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", cmd_args)
    args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    mintpy.geocode.main(args_str.split())

def process_saveroi(inps):    
    #save_roipac
    cmd_args = ['geo_temporalCoherence.h5', '-o', "".join(['geo_',inps.startDate,'_',inps.endDate,'.cor'])]    
    print("save_roipac.py", cmd_args)
    asct_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    os.system(multitrack_utilities.seperate_str_byspace(['save_roipac.py', asct_str.split()]))
    
    cmd_args = ['geo_geometryRadar.h5', 'height', '-o', 'srtm.dem']
    print("save_roipac.py", cmd_args)
    asct_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    os.system(multitrack_utilities.seperate_str_byspace(['save_roipac.py', asct_str.split()]))

def process_timeseries(inps):
    """geocode timeseries**.h5 file and get the deformation field of two time periods"""
    atr_asc = inps.file
   
    #unw file
    cmd_args = [atr_asc, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", cmd_args)
    args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    mintpy.geocode.main(args_str.split())
        
    #save dataset of unw cor and dem
    os.chdir("".join(inps.outdir))
    filename, extension = multitrack_utilities.seprate_filename_extension("".join(atr_asc))[1:3]
    
    cmd_args = ['geo_'+filename+extension, "".join([inps.startDate,'_',inps.endDate])]
    print("save_roipac.py", cmd_args)
    asct_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    os.system(multitrack_utilities.seperate_str_byspace(['save_roipac.py', asct_str.split()]))

    process_saveroi(inps)
    multitrack_utilities.delete_tmpgeo(inps.outdir,'geo_','.h5')

def process_ifgramStack(inps):
    """process ifgramStack.h5 file"""
    if os.path.exists("".join(inps.outdir))==False:
        os.makedirs("".join(inps.outdir))
    
    # dem file
    demname='geometryRadar.h5'
    if not os.path.isfile(demname):
        demname_f = './inputs/geometryRadar.h5'
    else:
        demname_f = 'geometryRadar.h5'
    cmd_args = [demname_f, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", cmd_args)
    args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    mintpy.geocode.main(args_str.split())
    
    #ifgramStack file
    atr_asc = ['./inputs/'+inps.file]
    cmd_args = [atr_asc, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", cmd_args)
    args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    mintpy.geocode.main(args_str.split())
        
    #save dataset of unw cor and dem
    os.chdir("".join(inps.outdir))
    filename, extension = multitrack_utilities.seprate_filename_extension("".join(atr_asc))[1:3]
    
    cmd_args = ['geo_'+filename+extension, "".join(['unwrapPhase-',inps.startDate,'_',inps.endDate])]
    print("save_roipac.py", cmd_args)
    asct_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    os.system(multitrack_utilities.seperate_str_byspace(['save_roipac.py', asct_str.split()]))   

    cmd_args = ['geo_'+filename+extension, "".join(['coherence-',inps.startDate,'_',inps.endDate])]
    print("save_roipac.py", cmd_args)
    asct_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    completion_status = os.system(multitrack_utilities.seperate_str_byspace(['save_roipac.py', asct_str.split()])) 
    
    cmd_args = ['geo_geometryRadar.h5', 'height', '-o', 'srtm.dem']
    print("save_roipac.py", cmd_args)
    asct_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    os.system(multitrack_utilities.seperate_str_byspace(['save_roipac.py', asct_str.split()]))
    
    multitrack_utilities.delete_tmpgeo(inps.outdir,'geo_','.h5')

def process_velocity(inps):
    """process velocity.h5 file"""
    atr_asc = inps.file
   
    #velocity file
    cmd_args = [atr_asc, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", cmd_args)
    args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    mintpy.geocode.main(args_str.split())
    
    os.chdir("".join(inps.outdir))
    process_saveroi(inps)
    print('save unw file')
    multitrack_utilities.velo_disp(inps)
    
    multitrack_utilities.delete_tmpgeo(inps.outdir,'geo_','.h5')

def process_HDFEOS(inps):
    """process *.he5 file"""
    atr_asc = inps.file
        
    #save dataset of unw cor and dem
    filename, extension = multitrack_utilities.seprate_filename_extension("".join(atr_asc))[1:3]
    cmd_args = [filename+extension, "".join(['displacement-',inps.startDate,'_',inps.endDate]), '-o', "".join(['geo_',inps.startDate,'_',inps.endDate,'.unw'])]
    print("save_roipac.py", cmd_args)
    asct_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    os.system(multitrack_utilities.seperate_str_byspace(['save_roipac.py', asct_str.split()]))   

    cmd_args = [filename+extension, 'temporalCoherence', '-o', "".join(['geo_',inps.startDate,'_',inps.endDate,'.cor'])]
    print("save_roipac.py", cmd_args)
    asct_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    completion_status=os.system(multitrack_utilities.seperate_str_byspace(['save_roipac.py', asct_str.split()])) 
    
    cmd_args = [filename+extension, 'height', '-o', 'srtm.dem']
    print("save_roipac.py", cmd_args)
    asct_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    os.system(multitrack_utilities.seperate_str_byspace(['save_roipac.py', asct_str.split()]))
    
    # mv 
    if os.path.exists(inps.outdir)==False:
        os.makedirs(inps.outdir)
    else:
        shutil.rmtree(inps.outdir)
        os.makedirs(inps.outdir)
    key1 = 'geo_'
    key2 = 'srtm'
    for file in os.listdir(os.getcwd()):
        if str.find(file,key1) != -1 or str.find(file,key2) != -1:
            shutil.move(file,inps.outdir) 
    os.chdir("".join(inps.outdir))  
    
    multitrack_utilities.delete_tmpgeo(inps.outdir,'geo_','.h5')    
    
######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   
    
    print('single track!')
    inps.startDate,inps.endDate =  multitrack_utilities.find_start_end_date(os.getcwd(),inps)
    
    print(inps)
    
    if str.find(inps.file,'ifgramStack') != -1:
        process_ifgramStack(inps)
    elif str.find(inps.file,'velocity') != -1:
        process_geocode(inps)
        process_velocity(inps)
    elif str.find(inps.file,'timeseries') != -1:
        process_geocode(inps)
        process_timeseries(inps)
    else:
        process_HDFEOS(inps)
    # rename *.rsc1 to *.rsc
    outfile = multitrack_utilities.seperate_str_byspace(['srtm.dem' + '.rsc'])
    write_rsc_file(inps,outfile,multitrack_utilities.seperate_str_byspace(['srtm.dem' +'.rsc1']))
    os.remove(outfile)
    print('rename *.rsc1 to *.rsc')
    os.rename(multitrack_utilities.seperate_str_byspace(['srtm.dem' +'.rsc1']),outfile)
    # generate dem.jpeg
    dem_jpeg('srtm.dem')


######################################################################################
if __name__ == '__main__':
    main()
