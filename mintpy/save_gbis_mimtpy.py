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

######################################################################################
EXAMPLE = """example:
  for singletrack:
  Note: startDate, endDate and outdir have default values, any of them can be not given:
    startDate default value = atr['START_DATE']
    endDate default value = atr['END_DATE']
    outdir default value = '$MODELDIR/project/Sen**/gbis_startDate_endDate/'
  save_gbis.py timeseries_ECMWF_demErr.h5 -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001 -m maskTempCoh.h5
  save_gbis.py ifgramStack.h5 -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001 -s 20171117 -e 20171129 -outdir $MODELDIR/Darbandikhan/SenAT73/
  save_gbis.py velocity.h5 -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001 -s 20171117 -e 20171129 -outdir $MODELDIR/Darbandikhan/SenAT73/
  save_gbis.py S1_IW23_026_0108_0113_20171117_XXXXXXXX.he5 -s 20171128 -e 20181210 
  for multitrack:
  Note: startDate, endDate and DataSet can be not given in template:
  save_gbis.py -t $MIMTFILES/Darbandikhan.txt
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Prepare data for GBIS software',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs='?', help='ascending or descending files\n')
    parser.add_argument('-t', '--template', dest='templateFile',
                        help="Template file with geocoding options.")                   
                        
    parser.add_argument('-ds', '--dataset', dest='DataSet',nargs='?',
                        help="name of dataset.Seperating by ','. ")
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
    parser.add_argument('-m', '--mask', dest='mask_file', help='mask file.')
    parser.add_argument('--ref-lalo', dest='ref_lalo', type=float, nargs='?',
                        help='custom reference pixel in lat/lon')
    
    
    parser.add_argument('-outdir','--outdir',dest='outdir',nargs=1,
                        help='output directory')

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    # read templatefile
    if inps.templateFile:
        inps = read_template2inps("".join(inps.templateFile), inps)
    else:    
        # default startDate and endDate
        if not os.path.isfile("".join(inps.file)):
            file = find_timeseries(os.getcwd())
        else:
            file = "".join(inps.file)    
        atr = readfile.read_attribute(file)
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
    
    prefix = 'gbis.'
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
            elif key in ['mask_file']:
                inps_dict[key] = value
            elif key in ['ref_lalo']:
                inps_dict[key] = list(tuple([float(i) for i in value.split(',')]))

    inps.laloStep = [inps.latStep, inps.lonStep]
    if None in inps.laloStep:
        inps.laloStep = None
    print(inps)
    return inps

def generate_outdir_name(inps,pwdDir):
    """set output directory"""
    if not inps.outdir or inps.outdir[0] == 'None':
        projectTrack = pwdDir.split('/')[-2] #get projectSenAT***
        ret = re.findall(r"^(.+)(Sen[AD]T\d+)$", projectTrack)
        project = ret[0][0]
        Track = ret[0][1]
        dirname = "".join([os.getenv('MODELDIR')+'/'+project+'/'+Track+'/gbis'+'_'+inps.startDate+'_'+inps.endDate])
    else:
        dirname = inps.outdir[0]+'/gbis'+'_'+inps.startDate+'_'+inps.endDate
    return dirname
    
def find_folder(tempfilename):
    """find the project folder and sort the folder in [*AT *DT]"""
    dir = os.getenv('SCRATCHDIR')
    folders = ["".join([dir +'/'])]
    project_folder = []
    for folder in folders:
        for x in os.listdir(folder):
            if os.path.isdir(os.path.join(folder,x)) and str.find(x,tempfilename)!= -1:
                project_folder.append(x)
    project_folder_sort = sorted(project_folder)
    return project_folder_sort
 
def find_timeseries(datadir):
    """find timeseries***.h5 file. The best results is timeseries_***_demErr.h5"""
    datafiles = []
    key1 = 'timeseries'
    key2 = 'Residual'
    for file in os.listdir(datadir):
        if os.path.splitext(file)[1] =='.h5':
            if str.find(file,key1) != -1 and str.find(file,key2) == -1:
                datafiles.append(file)
    datafile = []
    for file in datafiles:
        if len(file)>len(datafile):
            datafile = file
    return datafile

def find_S1_fullname(datadir):
    """find full name for S1 datatype"""
    datafiles = []
    key1 = 'S1'
    for file in os.listdir(datadir):
        if os.path.splitext(file)[1] =='.he5':
            if str.find(file,key1) != -1 :
                datafiles.append(file)
    return datafiles[0]
    
def track_date(datadir,date):
    """get the date close to the given date"""
    datafile = find_timeseries(datadir)
    completion_status = os.system(format_args(['info.py', datafile, '--date', '>', 'date_list.txt']))
    if completion_status == 1:
        raise Exception('error when runing info.py')

    if os.path.isfile('date_list.txt'):
        f = open('date_list.txt', 'r')
        lines = f.readlines() 
        f.close()        
    date_part = "".join(date)[0:6]
    date2 = []
    sub = 31
    for dates in lines:
        if str.find(dates,date_part) != -1:
            if abs(int(dates[6:8])-int(date[6:8]))<sub:
                sub = abs(int(dates[6:8])-int(date[6:8]))
                date2 = dates
    return date2.strip()

def find_date(datadir,inps):
    """find the startdate and enddate of each track"""   
    if not inps.startDate:
        startdate2 = inps.startDate
    if not inps.endDate:
        enddate2 = inps.endDate
    if inps.startDate:
        startdate2 = track_date(datadir,inps.startDate)
    if inps.endDate:
        enddate2 = track_date(datadir,inps.endDate)
    return startdate2,enddate2

def check_step(folders):
    """for S1*h5 file, check whether the lat_step and lon_step are same for different projects"""
    x_step = []
    y_step = []
    for project in folders:
       os.chdir("".join([os.getenv('SCRATCHDIR')+'/'+project+'/PYSARTEST/']))
       # find S1*.h5 whole name
       datafile = find_S1_fullname("".join([os.getenv('SCRATCHDIR')+'/'+project+'/PYSARTEST/']))
       atr = readfile.read_attribute(datafile)
       x_step.append(atr['X_STEP'])
       y_step.append(atr['Y_STEP'])
    if len(set(x_step))!= 1:
        raise Exception("error! lon_Step between different tracks is not same!")
    elif len(set(y_step))!= 1:
        raise Exception("error! lat_Step between different tracks is not same!")
    else:
        return True

def multitrack_run_save_gbis(inps,folders):
    """run save_gbis for each track"""
    for project in folders:        
        os.chdir("".join([os.getenv('SCRATCHDIR')+'/'+project+'/PYSARTEST/']))
        if inps.DataType == 'S1':
            datafile = find_S1_fullname("".join([os.getenv('SCRATCHDIR')+'/'+project+'/PYSARTEST/']))
        elif inps.DataType == 'timeseries':
            datafile = find_timeseries("".join([os.getenv('SCRATCHDIR')+'/'+project+'/PYSARTEST/']))
        elif inps.DataType == 'ifgramStack':
            datafile = "".join([str(inps.DataType)+'.h5'])
        elif inps.DataType == 'velocity':
            datafile = "".join([str(inps.DataType)+'.h5'])                    
        
        if not inps.ref_lalo:
            print(format_args(['save_gbis_mimtpy.py', datafile, '-b', inps.SNWE, '-y', inps.latStep, '-x', inps.lonStep, '-s', inps.startDate, '-e', inps.endDate, '-m', inps.mask_file, '-outdir', inps.outdir]))
            completion_status = os.system(format_args(['save_gbis_mimtpy.py', datafile, '-b', inps.SNWE, '-y', inps.latStep, '-x', inps.lonStep, '-s', inps.startDate, '-e', inps.endDate, '-m', inps.mask_file, '-outdir', inps.outdir]))
            if completion_status == 1:
                raise Exception('error when runing save_gbis.py')        
        else:
            print(format_args(['save_gbis_mimtpy.py', datafile, '-b', inps.SNWE, '-y', inps.latStep, '-x', inps.lonStep, '-s', inps.startDate, '-e', inps.endDate, '-m', inps.mask_file, '--ref-lalo', inps.ref_lalo, '-outdir', inps.outdir]))
            completion_status = os.system(format_args(['save_gbis_mimtpy.py', datafile, '-b', inps.SNWE, '-y', inps.latStep, '-x', inps.lonStep, '-s', inps.startDate, '-e', inps.endDate, '-m', inps.mask_file, '--ref-lalo', inps.ref_lalo, '-outdir', inps.outdir]))
            if completion_status == 1:
                raise Exception('error when runing save_gbis.py')

def run_save_gbis(inps):
    """run save_gbis.py in proper directory"""
    if not inps.DataSet:
        tempfilename = inps.templateFile
        folders = find_folder(seprate_filename_exten(tempfilename)[1])
        print(folders)
    else:
        folders = inps.DataSet
        print(folders)
    # if datatype is S1,first judge whether they have same lat_step and lon_step
    if inps.DataType == 'S1':
        flag = check_step(folders)
        if flag:
            multitrack_run_save_gbis(inps,folders)
    else:
        multitrack_run_save_gbis(inps,folders)
    
def format_args(arr, dst=""):
    """ parse list array(list item or values) to string """
    for k in arr:
        if isinstance(k,list):
            dst = dst + " " + format_args(k)
        else:
            dst = dst + " " + str(k)
    return dst.strip()
   
def seprate_filename_exten(path):
    """ return directory(filepath), filename(shotname), extension """
    (filepath, tempfilename) = os.path.split(os.path.abspath(path))
    (filename, extension) = os.path.splitext(tempfilename)
    return filepath, filename, extension

def delete_tmpgeo(datadir, key1, key2):
    """delete all geo_*.h5 files in $MODLEDIR/project/SenAT(DT)/gbis_startdate_enddate/"""
    for file in os.listdir(datadir):
        if os.path.splitext(file)[1] ==key2:
            if str.find(file,key1) != -1:
                os.remove(datadir+'/'+file)

def velo_disp(inps):
    """calculated displacement during startDate_endDate period based on linear assumption and velocity.h5"""
    data, atr = readfile.read('geo_velocity.h5')
    # calculate disp
    dt1, dt2 = ptime.date_list2vector([inps.startDate, inps.endDate])[0]
    data *= (dt2 - dt1).days / 365.25
    # displacement to phase
    range2phase =  -4. * np.pi / float(atr['WAVELENGTH'])
    data *= range2phase
    # write atr
    atr['PROCESSOR'] = 'roipac'
    atr['FILE_TYPE'] = '.unw'
    atr['UNIT'] = 'radian'
    out_file = 'geo_'+'{}_{}.unw'.format(inps.startDate, inps.endDate)
    writefile.write(data, out_file=out_file, metadata=atr)

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
    args_str = format_args(cmd_args)
    mintpy.geocode.main(args_str.split())
    
    if inps.mask_file:
        cmd_args = [inps.mask_file, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
        print("geocode.py", cmd_args)
        args_str = format_args(cmd_args)
        mintpy.geocode.main(args_str.split())

def process_time(inps):
    """geocode timeseries**.h5 file and get the deformation field of two time periods"""
    atr_asc = inps.file
   
    #unw file
    cmd_args = [atr_asc, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", cmd_args)
    args_str = format_args(cmd_args)
    mintpy.geocode.main(args_str.split())
        
    #save dataset of unw
    os.chdir("".join(inps.outdir))
    filename, extension = seprate_filename_exten("".join(atr_asc))[1:3]
    
    cmd_args = ['geo_'+filename+extension, "".join([inps.startDate,'_',inps.endDate])]
    print("save_roipac.py", cmd_args)
    asct_str = format_args(cmd_args)
    os.system(format_args(['save_roipac.py', asct_str.split()]))


def process_ifgS(inps):
    """process ifgramStack.h5 file"""
    if os.path.exists("".join(inps.outdir))==False:
        os.makedirs("".join(inps.outdir))    

    #ifgramStack file
    atr_asc = ['./inputs/'+inps.file]
    cmd_args = [atr_asc, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", cmd_args)
    args_str = format_args(cmd_args)
    mintpy.geocode.main(args_str.split())
        
    #save dataset of unw
    os.chdir("".join(inps.outdir))
    filename, extension = seprate_filename_exten("".join(atr_asc))[1:3]
    
    cmd_args = ['geo_'+filename+extension, "".join(['unwrapPhase-',inps.startDate,'_',inps.endDate])]
    print("save_roipac.py", cmd_args)
    asct_str = format_args(cmd_args)
    os.system(format_args(['save_roipac.py', asct_str.split()]))   

def process_vel(inps):
    """process velocity.h5 file"""
    atr_asc = inps.file
   
    #velocity file
    cmd_args = [atr_asc, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", cmd_args)
    args_str = format_args(cmd_args)
    mintpy.geocode.main(args_str.split())
    
    os.chdir("".join(inps.outdir))
    print('save unw file')
    velo_disp(inps)

def process_S1(inps):
    """process S1*.h5 file"""

    atr_asc = inps.file
        
    #save dataset of unw and geometry
    filename, extension = seprate_filename_exten("".join(atr_asc))[1:3]
    cmd_args = [filename+extension, "".join(['displacement-',inps.startDate,'_',inps.endDate]), '-o', "".join(['geo_',inps.startDate,'_',inps.endDate,'.unw'])]
    print("save_roipac.py", cmd_args)
    asct_str = format_args(cmd_args)
    os.system(format_args(['save_roipac.py', asct_str.split()]))   
    
    # mv 
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
        geom_file = find_S1_fullname(os.getcwd())
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
    delete_tmpgeo(inps.outdir, 'S1_', '.he5')
    delete_tmpgeo(inps.outdir, 'geo_', '.h5')
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
    print(inps)    
    if not inps.templateFile:
        print('single track!')
        inps.startDate,inps.endDate = find_date(os.getcwd(),inps)
        inps.outdir = generate_outdir_name(inps,os.getcwd())
        print(inps)
        if str.find(inps.file,'ifgramStack.h5') != -1:
            process_geocode(inps)
            process_ifgS(inps)
        elif str.find(inps.file,'velocity.h5') != -1:
            process_geocode(inps)
            process_vel(inps)
        elif str.find(inps.file,'timeseries') != -1:
            process_geocode(inps)
            process_time(inps)
        else:
            process_S1(inps)
        prep_gbis(inps)
        save2mat(inps)
    else:
        print('multi track!')
        run_save_gbis(inps)

    
######################################################################################
if __name__ == '__main__':
    main()
