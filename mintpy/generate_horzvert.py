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
  Without template file
  Note: startDate, endDate and outdir have default values, any of them can be not given:
    startDate default value = atr['START_DATE']
    endDate default value = atr['END_DATE']
    outdir default value = '$MODELDIR/projectSen/horzvert/'
  generate_horzvert.py DarbandikhanSenAT73 DarbandikhanSenDT80 -dt velocity -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001 
  generate_horzvert.py DarbandikhanSenAT73 DarbandikhanSenDT80 -dt timeseries -b 34.2 35.2 45.0 46.3 -y 0.001 -x  0.001 -m maskTempCoh.h5 -az 45 --outname horizontal.h5 vertical.h5
  generate_horzvert.py DarbandikhanSenAT73 DarbandikhanSenDT80 -dt ifgramStack -b 34.2 35.2 45.0 46.3 -y 0.001  -x 0.001 -e 20171130 --outname horizontal.h5 vertical.h5
  With template file:
  Note: startDate, endDate and DataSet can be not given in template:
  generate_horzvert.py -t $MIMTFILES/Darbandikhan.txt
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Generate horizontal and vertical datafile',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('DataSet', nargs='*', help='ascending or descending files\n')
    parser.add_argument('-t','--template', dest='templateFile',
                        help="Template file with geocoding options.")
    parser.add_argument('--model-software', dest='model_software',
                        help="Goephysical model's name.")                    
                        
    #parser.add_argument('-ds', '--dataset', dest='DataSet',nargs=2,metavar=('AT_FILE', 'DT_FILE'),
    #                    help='name of dataset.')
    parser.add_argument('-dt', '--datatype',dest='DataType',nargs='?',
                        help='clarify the type of data.[velocity, ifgramStack, timeseries, S1*.he5].Only used in template file')
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
    parser.add_argument('-m', '--mask', dest='maskfile', help='mask file.')
    parser.add_argument('-az', '--azimuth', dest='azimuth', type=float, default=90.0,
                        help='azimuth angle in degree (clockwise) of the direction of the horizontal movement\n' +
                             'default is 90.0 for E-W component, assuming no N-S displacement.\n' +
                             'i.e. azimuth angle of strike-slip fault\n\n')
    
    parser.add_argument('--outname',dest='outname',nargs=2,metavar=('HZ_FILE','UP_FILE'),default=['hz.h5', 'up.h5'],
                        help='output file name for vertical adn horizontal components')
    parser.add_argument('-outdir','--outdir',dest='outdir',nargs='?',default=os.getenv('MODELDIR'),
                        help='output directory')

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    # read templatefile
    if inps.templateFile:
        inps = read_template2inps("".join(inps.templateFile), inps)
 
    # default startDate and endDate
    tempath = os.getenv('SCRATCHDIR')+'/'+inps.DataSet[0]+'/mintpy/'
    file = multitrack_utilities.find_timeseries_horzvert(tempath)  
    
    atr = readfile.read_attribute(tempath+file)
    if not inps.startDate or inps.startDate=='None':
        inps.startDate = atr['START_DATE']
    if not inps.endDate or inps.endDate=='None':
        inps.endDate = atr['END_DATE']
    return inps    


def read_template2inps(templatefile, inps):
    """Read input template options into Namespace inps"""
    print('read input option from template file: ' + templatefile)
    if not inps:
        inps = cmd_line_parse()
    inps_dict = vars(inps)

    template = readfile.read_template(templatefile)    
    template = ut.check_template_auto_value(template)
    
    prefix = 'horzvert.'
    key_list = [i for i in list(inps_dict.keys()) if prefix + i in template.keys()]
    for key in key_list:
        value = template[prefix + key]
        if value:
            if key == 'DataSet':
                inps_dict[key] = list(tuple([i for i in value.split(',')]))
            elif key == 'DataType':
                inps_dict[key] = value
            elif key == 'SNWE':
                inps_dict[key] = list(tuple([float(i) for i in value.split(',')]))
            elif key in ['latStep', 'lonStep']:
                inps_dict[key] = float(value)
            elif key in ['startDate','endDate']:
                inps_dict[key] = value
            elif key in ['maskfile']:
                inps_dict[key] = value
            elif key in ['azimuth']:
                inps_dict[key] = value
            elif key in ['outname']:
                inps_dict[key] = list(tuple([i for i in value.split(',')]))

    inps.laloStep = [inps.latStep, inps.lonStep]
    if None in inps.laloStep:
        inps.laloStep = None
    return inps

def generate_outdir_name(inps,project_name):
    """set output directory"""
    dirname = inps.outdir+'/'+project_name+'/horzvert/'+inps.startDate+'_'+inps.endDate+'/'
    return dirname
    
def copy_file(key1,key2,outdir):
    """mv data file"""
    for file in os.listdir(os.getcwd()):
        if str.find(file,key1) != -1 and str.find(file,key2) != -1:
            print('move file')
            if os.path.isfile(outdir+file) == True:
                os.remove(outdir+file)
            shutil.move(file,outdir)
    
def velocity_displacement(atr_asc, startDate, endDate):
    """calculated displacement during startDate_endDate period based on linear assumption and velocity.h5"""
    data, atr = readfile.read(atr_asc)
    # calculate disp
    dt1, dt2 = ptime.date_list2vector([startDate, endDate])[0]
    data *= (dt2 - dt1).days / 365.25
    # displacement to phase
    range2phase =  -4. * np.pi / float(atr['WAVELENGTH'])
    data *= range2phase
    # write atr
    atr['PROCESSOR'] = 'roipac'
    atr['FILE_TYPE'] = '.unw'
    atr['UNIT'] = 'radian'
    out_file = 'geo_'+'{}_{}.unw'.format(startDate, endDate)
    writefile.write(data, out_file=out_file, metadata=atr)

def process_timeseries(inps, track, startDate, endDate, outdir, tmpdir):
    """geocode timeseries**.h5 file and get the deformation field of two time periods"""
    if os.path.exists("".join(tmpdir))==False:
        os.makedirs("".join(tmpdir))
        
    atr_asc = multitrack_utilities.find_timeseries_horzvert("".join([os.getenv('SCRATCHDIR')+'/'+track+'/mintpy/']))
    
    # mask
    if inps.maskfile and inps.maskfile != 'None':
        cmd_args = [atr_asc, '-m', inps.maskfile]
        print('mask.py', cmd_args)
        args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
        os.system(multitrack_utilities.seperate_str_byspace(['mask.py', args_str.split()]))
        atr_asc = multitrack_utilities.seprate_filename_extension(atr_asc)[1]+'_msk.h5'
    
    #geocode timeseries file
    cmd_args = [atr_asc, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir', tmpdir]
    print("geocode.py", cmd_args)
    args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    mintpy.geocode.main(args_str.split())
        
    #save dataset of unw
    os.chdir("".join(tmpdir))
    filename, extension = multitrack_utilities.seprate_filename_extension("".join(atr_asc))[1:3]    
    cmd_args = ['geo_'+filename+extension, "".join([startDate,'_',endDate])]
    print("save_roipac.py", cmd_args)
    args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    os.system(multitrack_utilities.seperate_str_byspace(['save_roipac.py', args_str.split()]))
    
    copy_file('geo','unw',outdir)
    shutil.rmtree(tmpdir)

def process_ifgramStack(inps, track, startDate, endDate, outdir, tmpdir):
    """process ifgramStack.h5 file"""
    if os.path.exists("".join(tmpdir))==False:
        os.makedirs("".join(tmpdir)) 

    atr_asc = ['./inputs/'+str(inps.DataType)+'.h5']
    # mask
    if inps.maskfile and inps.maskfile != 'None':
        date_msk1 = re.split(r'[_.]', str(inps.maskfile))[0]
        date_msk1f = multitrack_utilities.find_nearest_date(os.getcwd(),str(date_msk1))
        date_msk2 = re.split(r'[_.]', str(inps.maskfile))[1]
        date_msk2f = multitrack_utilities.find_nearest_date(os.getcwd(),str(date_msk2))
        maskfile = date_msk1f+'_'+date_msk2f +'.cor'
        cmd_args = [atr_asc, '-m', maskfile]
        print('mask.py', cmd_args)
        args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
        os.system(multitrack_utilities.seperate_str_byspace(['mask.py', args_str.split()]))
        atr_asc = './inputs/'+multitrack_utilities.seprate_filename_extension("".join(atr_asc))[1]+'_msk.h5'
    
    #geocode ifgramStack file
   
    cmd_args = [atr_asc, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(tmpdir)]
    print("geocode.py", cmd_args)
    args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    mintpy.geocode.main(args_str.split())
        
    #save dataset of unw
    os.chdir("".join(tmpdir))
    filename, extension = multitrack_utilities.seprate_filename_extension("".join(atr_asc))[1:3]    
    cmd_args = ['geo_'+filename+extension, "".join(['unwrapPhase-',startDate,'_',endDate])]
    print("save_roipac.py", cmd_args)
    args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    os.system(multitrack_utilities.seperate_str_byspace(['save_roipac.py', args_str.split()]))  

    copy_file('geo','unw',outdir)
    shutil.rmtree(tmpdir)

def process_HDFEOS(inps, track, startDate, endDate, outdir, tmpdir):
    """process S1*.h5 file"""

    atr_asc = multitrack_utilities.find_HDFEOS_fullname("".join([os.getenv('SCRATCHDIR')+'/'+track+'/mintpy/']))
        
    #save dataset of unw 
    filename, extension = multitrack_utilities.seprate_filename_extension("".join(atr_asc))[1:3]
    cmd_args = [filename+extension, "".join(['displacement-',startDate,'_',endDate]), '-o', "".join(['geo_',startDate,'_',endDate,'.unw'])]
    print("save_roipac.py", cmd_args)
    args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    os.system(multitrack_utilities.seperate_str_byspace(['save_roipac.py', args_str.split()]))   
    
    # mv 
    if os.path.exists(tmpdir)==False:
        os.makedirs(tmpdir)
    else:
        shutil.rmtree(tmpdir)
        os.makedirs(tmpdir)
    key1 = 'geo_'
    for file in os.listdir(os.getcwd()):
        if str.find(file,key1) != -1:
            shutil.move(file,tmpdir) 
    
    copy_file('geo','unw',outdir) 
    shutil.rmtree(tmpdir)    

def process_velocity(inps, track, startDate, endDate, outdir, tmpdir):
    """process velocity.h5 file"""
    if os.path.exists("".join(tmpdir))==False:
        os.makedirs("".join(tmpdir))
        
    atr_asc = "".join([str(inps.DataType)+'.h5'])

    # mask
    if inps.maskfile and inps.maskfile != 'None':
        cmd_args = [atr_asc, '-m', inps.maskfile]
        print('mask.py', cmd_args)
        args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
        os.system(multitrack_utilities.seperate_str_byspace(['mask.py', args_str.split()]))
        atr_asc = multitrack_utilities.seprate_filename_extension(atr_asc)[1]+'_msk.h5'
        print(atr_asc)
    #geocode velocity file
    cmd_args = [atr_asc, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(tmpdir)]
    print("geocode.py", cmd_args)
    args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    mintpy.geocode.main(args_str.split())
    
    # save dataset of unw
    os.chdir("".join(tmpdir))
    print('save unw file')
    atr_asc = 'geo_' + atr_asc
    velocity_displacement(atr_asc,startDate, endDate)
    
    copy_file('geo','unw',outdir)
    return atr_asc

def calculate_hozrvert_velocity(atr_asc,outdir,inps):
    """generate horizontal and vertical velocity data"""
    os.chdir("".join(outdir))
    
    folders = multitrack_utilities.find_folder_horzvert('Sen',outdir)
    print(folders)
    file1 = outdir+folders[0]+'/'+atr_asc
    print(file1)
    file2 = outdir+folders[1]+'/'+atr_asc
    print(file2)
    
    outname = ['vel_horizontal.h5','vel_vertical.h5']
    if inps.azimuth != 'None':
        cmd_args = [file1, file2, '--az', inps.azimuth, '-o', outname]
    else:
        cmd_args = [file1, file2, '-o', outname]
    print("asc_desc2horz_vert.py", cmd_args)
    args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    os.system(multitrack_utilities.seperate_str_byspace(['asc_desc2horz_vert.py', args_str.split()])) 

    for folder in folders:
        delete_dir = outdir+folder+'/'
        shutil.rmtree(delete_dir)
    
def calculate_horzvert_displacement(inps,outdir):
    """generate horizontal and vertical displacment data"""
    os.chdir(outdir)
    
    datafiles = []
    key1 = 'geo_'
    for file in os.listdir(outdir):
        if os.path.splitext(file)[1] =='.unw':
            if str.find(file,key1) != -1:
                datafiles.append(file)
    if inps.azimuth != 'None':
        if inps.outname[0] != 'None':
            cmd_args = [datafiles[0], datafiles[1], '--az', inps.azimuth, '-o', inps.outname]
        else:
            cmd_args = [datafiles[0], datafiles[1], '--az', inps.azimuth]
    else:
        if inps.outname[0] != 'None':
            cmd_args = [datafiles[0], datafiles[1], '-o', inps.outname]
        else:
            cmd_args = [datafiles[0], datafiles[1]]
    print("asc_desc2horz_vert.py", cmd_args)
    args_str = multitrack_utilities.seperate_str_byspace(cmd_args)
    os.system(multitrack_utilities.seperate_str_byspace(['asc_desc2horz_vert.py', args_str.split()])) 
    return
      
    
######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    print(inps)    

    # if datatype is S1,first judge whether they have same lat_step and lon_step
    if inps.DataType == 'HDFEOS':
        flag = multitrack_utilities.check_X_Y_step(inps.DataSet)
    # get project name
    project_name,project_length = multitrack_utilities.find_intersection_part(inps.DataSet[0],inps.DataSet[1])
    # generate outdir
    outdir = generate_outdir_name(inps,project_name)
    # run each track
    for project in inps.DataSet:
        print(project)
        # get project and track name
        ret = re.findall(r"^(.+?)(\d+)$", project)
        satellite = project_name[-3:]
        track_No = ret[0][1]
        
        # set parameter for each track
        os.chdir(os.getenv('SCRATCHDIR')+'/'+project+'/mintpy/')
        startDate,endDate = multitrack_utilities.find_start_end_date(os.getcwd(),inps)
        tmpdir = outdir + satellite + track_No + '/'

        # process data
        if inps.DataType == 'ifgramStack':   
            process_ifgramStack(inps, project,startDate,endDate, outdir, tmpdir)
        elif inps.DataType == 'velocity':    
            atr_asc = process_velocity(inps, project,startDate,endDate, outdir, tmpdir)
        elif inps.DataType == 'timeseries':    
            process_timeseries(inps, project,startDate,endDate, outdir, tmpdir)
        else:
            if flag:
                process_HDFEOS(inps, project,startDate,endDate, outdir, tmpdir)
    # generate horizontal and vertical deformation           
    if inps.DataType == 'velocity':
        calculate_hozrvert_velocity(atr_asc,outdir,inps)
    calculate_horzvert_displacement(inps,outdir)
    
######################################################################################
if __name__ == '__main__':
    main()
