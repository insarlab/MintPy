#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# Yunjun, Jul 2015: Add check_num/check_file_size to .int/.cor file
# Yunjun, Jan 2017: Add auto_path_miami(), copy_roipac_file()
#                   Add roipac2pysar_multi_group_hdf5()
#                   Add r+ mode loading of multi_group hdf5 file


import os
import sys
import glob
import argparse
import warnings

import h5py
import numpy as np
import shutil

import pysar
import pysar._readfile as readfile
import pysar._pysar_utilities as ut
from pysar._readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file


############################ Sub Functions ###################################
##################################################################
def auto_path_miami(inps, template_dict={}):
    '''Auto File Path Setting for Geodesy Lab - University of Miami'''
    print 'Use auto path setting in University of Miami.'+\
          '(To turn it off, change miami_path value to False in pysar/__init__.py)'
    # PYSAR working directory
    if not inps.timeseries_dir:
        inps.timeseries_dir = os.getenv('SCRATCHDIR')+'/'+inps.project_name+'/PYSAR'

    # .unw/.cor/.int files
    process_dir = os.getenv('SCRATCHDIR')+'/'+inps.project_name+'/PROCESS'
    print "PROCESS directory: "+process_dir
    if not inps.unw:   inps.unw = process_dir+'/DONE/IFGRAM*/filt_*.unw'
    if not inps.cor:   inps.cor = process_dir+'/DONE/IFGRAM*/filt_*rlks.cor'
    if not inps.int:   inps.int = process_dir+'/DONE/IFGRAM*/filt_*rlks.int'

    # master interferogram for geomap*.trans and DEM in radar coord
    try:    m_date12 = np.loadtxt(process_dir+'/master_ifgram.txt', dtype=str).tolist()
    except: m_date12 = os.walk(process_dir+'/GEO').next()[1][0].split('geo_')[1]

    if not inps.geomap:
        try: inps.geomap = process_dir+'/GEO/*'+m_date12+'*/geomap*.trans'
        except: warnings.warn('Can not locate geomap*.trans file for geocoding!')
    
    if not inps.dem_radar:
        try: inps.dem_radar = process_dir+'/DONE/*'+m_date12+'*/radar*.hgt'
        except: warnings.warn('Can not locate DEM in radar coord!')

    # Use DEMg/DEM option if dem_geo is not specified in pysar option
    dem_dir = os.getenv('SCRATCHDIR')+'/'+inps.project_name+'/DEM'
    if not inps.dem_geo:
        if os.path.isdir(dem_dir):            inps.dem_geo = dem_dir+'*.dem'
        elif 'DEMg' in template_dict.keys():  inps.dem_geo = template_dict['DEMg']
        elif 'DEM'  in template_dict.keys():  inps.dem_geo = template_dict['DEM']
        else:  warnings.warn('Can not locate DEM in geo coord!')

    return inps


def mode (thelist):
    '''Find Mode (most common) item in the list'''
    if not thelist:
        return None
    if len(thelist) == 1:
        return thelist[0]

    counts = {}
    for item in thelist:
        counts[item] = counts.get(item, 0) + 1
    maxcount = 0
    maxitem  = None
    for k, v in counts.items():
        if v > maxcount:
            maxitem  = k
            maxcount = v
    
    if maxcount == 1:
        print "All values only appear once"
        return None
    elif counts.values().count(maxcount) > 1:
        print "List has multiple modes"
        return None
    else:
        return maxitem


##################################################################
def check_file_size(fileList, mode_width=None, mode_length=None):
    '''Update file list and drop those not in the same size with majority.'''
    # If input file list is empty
    if not fileList:
        return fileList, None, None

    # Read Width/Length list
    widthList =[]
    lengthList=[]
    for file in fileList:
        rsc = readfile.read_attribute(file)
        widthList.append(rsc['WIDTH'])
        lengthList.append(rsc['FILE_LENGTH'])
    # Mode of Width and Length
    if not mode_width:   mode_width  = mode(widthList)
    if not mode_length:  mode_length = mode(lengthList)
    
    # Update Input List
    ext = os.path.splitext(fileList[0])[1]
    fileListOut = list(fileList)
    if widthList.count(mode_width)!=len(widthList) or lengthList.count(mode_length)!=len(lengthList):
        print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        print 'WARNING: Some '+ext+' may have the wrong dimensions!'
        print 'All '+ext+' should have the same size.'
        print 'The width and length of the majority of '+ext+' are: '+str(mode_width)+', '+str(mode_length)
        print 'But the following '+ext+' have different dimensions and thus will not be loaded:'
        for i in range(len(fileList)):
            if widthList[i] != mode_width or lengthList[i] != mode_length:
                print fileList[i]+'    width: '+widthList[i]+'  length: '+lengthList[i]
                fileListOut.remove(fileList[i])
        print '\nNumber of '+ext+' left: '+str(len(fileListOut))
        print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    return fileListOut, mode_width, mode_length


def check_existed_hdf5_file(roipacFileList, hdf5File):
    '''Check file list with existed hdf5 file'''
    # If input file list is empty
    outFileList = list(roipacFileList)
    if not outFileList:
        return outFileList
    
    # if previous hdf5 file existed
    if os.path.isfile(hdf5File):
        print os.path.basename(hdf5File)+'  already exists.'
        atr = readfile.read_attribute(hdf5File)
        k = atr['FILE_TYPE']
        h5 = h5py.File(hdf5File, 'r')
        epochList = sorted(h5[k].keys())
        h5.close()
        
        # Remove file/epoch that already existed
        for epoch in epochList:
            for file in roipacFileList:
                if epoch in file:
                    outFileList.remove(file)

        # Check mode length/width with existed hdf5 file
        if outFileList:
            ext = os.path.splitext(outFileList[0])[1]
            outFileList, mode_width, mode_length = check_file_size(outFileList)
            if mode_width != atr['WIDTH'] or mode_length != atr['FILE_LENGTH']:
                print 'WARNING: input ROI_PAC files have different size than existed hdf5 file:'
                print 'ROI_PAC file size: '+mode_length+', '+mode_width
                print 'HDF5    file size: '+atr['FILE_LENGTH']+', '+atr['WIDTH']
                print 'Continue WITHOUT loading '+ext+' file'
                print 'To enforse loading, change/remove existed HDF5 filename and re-run loading script'
                outFileList = None
    
    return outFileList


def roipac2multi_group_hdf5(fileType, fileList, hdf5File='unwrapIfgram.h5', extra_meta_dict=dict()):
    '''Load multiple ROI_PAC files into HDF5 file (Multi-group, one dataset and one attribute dict per group).
    Inputs:
        fileType : string, i.e. interferograms, coherence, snaphu_connect_component, etc.
        fileList : list of path, ROI_PAC .unw/.cor/.int/.byt file
        hdf5File : string, file name/path of the multi-group hdf5 PySAR file
        extra_meta_dict : dict, extra attribute dictionary 
    Outputs:
        hdf5File : output hdf5 file name
        fileList : list of string, files newly added

    '''
    ext = os.path.splitext(fileList[0])[1]
    print 'loading ROI_PAC '+ext+' files into '+fileType+' HDF5 file ...'
    print 'number of '+ext+' input: '+str(len(fileList))

    # Check width/length mode of input files
    fileList, mode_width, mode_length = check_file_size(fileList)
    if not fileList:
        return None, None

    # Check conflict with existing hdf5 file
    fileList2 = check_existed_hdf5_file(fileList, hdf5File)
    
    # Open(Create) HDF5 file with r+/w mode based on fileList2
    
    if fileList2 == fileList:
        # Create and open new hdf5 file with w mode
        print 'number of '+ext+' to add: '+str(len(fileList))
        print 'open '+hdf5File+' with w mode'
        h5file = h5py.File(hdf5File, 'w')
    elif fileList2:
        # Open existed hdf5 file with r+ mode
        print 'Continue by adding the following new epochs ...'
        print 'number of '+ext+' to add: '+str(len(fileList))
        print 'open '+hdf5File+' with r+ mode'
        h5file = h5py.File(hdf5File, 'r+')
        fileList = list(fileList2)
    else:
        print 'All input '+ext+' are included, no need to re-load.'
        fileList = None

    # Loop - Writing ROI_PAC files into hdf5 file
    if fileList:
        # Unwraped Interferograms
        if not fileType in h5file.keys():
            gg = h5file.create_group(fileType)     # new hdf5 file
        else:
            gg = h5file[fileType]                  # existing hdf5 file
        
        for file in fileList:
            print 'Adding ' + file
            data, rsc = readfile.read(file)
            
            # Dataset
            group = gg.create_group(os.path.basename(file))
            dset = group.create_dataset(os.path.basename(file), data=data, compression='gzip')
            
            # Attribute - *.unw.rsc
            for key,value in rsc.iteritems():
                group.attrs[key] = value
            # Attribute - *baseline.rsc
            d1, d2 = rsc['DATE12'].split('-')
            baseline_file = os.path.dirname(file)+'/'+d1+'_'+d2+'_baseline.rsc'
            baseline_rsc = readfile.read_roipac_rsc(baseline_file)
            for key,value in baseline_rsc.iteritems():
                group.attrs[key] = value
            # Attribute - PySAR
            try: group.attrs['PROJECT_NAME'] = extra_meta_dict['project_name']
            except: pass
        
        # End of Loop
        h5file.close()
        print 'finished writing to '+hdf5File

    return hdf5File, fileList


def roipac_nonzero_mask(unwFileList, maskFile='mask.h5'):
    '''Generate mask for non-zero amplitude pixel of ROI_PAC .unw file list.'''
    unwFileList, width, length = check_file_size(unwFileList)
    if unwFileList:
        # Initial mask value
        maskZero = np.ones([int(length), int(width)])

        # Update mask from input .unw file list
        fileNum = len(unwFileList)
        for i in range(fileNum):
            file = unwFileList[i]
            amp, unw, rsc = readfile.read_float32(file)
            
            maskZero *= amp
            ut.print_progress(i+1, fileNum, prefix='loading', suffix=os.path.basename(file))
        mask = np.ones([int(length), int(width)])
        mask[maskZero==0] = 0
        
        # write mask hdf5 file
        print 'writing >>> '+maskFile
        h5 = h5py.File(maskFile,'w')
        group = h5.create_group('mask')
        dset = group.create_dataset('mask', data=mask, compression='gzip')
        # Attribute - *.unw.rsc
        for key,value in rsc.iteritems():
            group.attrs[key] = value
        # Attribute - *baseline.rsc
        d1, d2 = rsc['DATE12'].split('-')
        baseline_file = os.path.dirname(file)+'/'+d1+'_'+d2+'_baseline.rsc'
        baseline_rsc = readfile.read_roipac_rsc(baseline_file)
        for key,value in baseline_rsc.iteritems():
            group.attrs[key] = value

    return maskFile, unwFileList


def roipac2single_dataset_hdf5(file_type, infile, outfile, extra_meta_dict=dict()):
    '''Convert ROI_PAC .dem / .hgt file to hdf5 file
    Based on load_dem.py written by Emre Havazli
    Inputs:
        file_type : string, group name of hdf5 file, i.e. dem, mask
        infile    : string, input ROI_PAC file name
        outfile   : string, output hdf5 file name
        extra_meta_dict : dict, extra attributes to output file
    Output:
        outfile   : string, output hdf5 file name
    '''
    # Read input file
    data, atr = readfile.read(infile)
    
    # Write output file - data
    print 'writing >>> '+outfile
    h5 = h5py.File(outfile, 'w')
    group = h5.create_group(file_type)
    dset = group.create_dataset(file_type, data=data, compression='gzip')
    
    # Write output file - attributes
    for key, value in atr.iteritems():
        group.attrs[key] = value
    try: group.attrs['PROJECT_NAME'] = extra_meta_dict['project_name']
    except: pass
    
    h5.close()
    return outfile


def copy_file(targetFile, destDir):
    '''Copy file and its .rsc/.par/.xml file to destination directory.'''
    #print '--------------------------------------------'
    destFile = destDir+'/'+os.path.basename(targetFile)
    if ut.update_file(destFile, targetFile):
        print 'copy '+targetFile+' to '+destDir
        shutil.copy2(targetFile, destDir)
        try: shutil.copy2(targetFile+'.rsc', destDir)
        except: pass
        try: shutil.copy2(targetFile+'.xml', destDir)
        except: pass
        try: shutil.copy2(targetFile+'.par', destDir)
        except: pass
    return destFile


def load_file(fileList, inps_dict=dict(), outfile=None, file_type=None):
    '''Load input file(s) into one HDF5 file 
    It supports ROI_PAC files only for now.
    Inputs:
        fileList  - list of string, path of files to load
        inps_dict - dict, including the following attributes
                    PROJECT_NAME   : KujuAlosAT422F650  (extra attribute dictionary to add to output file)
                    timeseries_dir : directory of time series analysis, e.g. KujuAlosAT422F650/PYSAR
        outfile   - string, output file name
        file_type - string, group name for output HDF5 file, interferograms, coherence, dem, etc.
    Output:
        outfile - string, output file name
    Example:
        unwrapIfgram.h5 = load_file('filt*.unw', inps_dict=vars(inps))
    '''
    # Get project_name from input template file
    if not 'project_name' in inps_dict.keys() and 'template_file' in inps_dict.keys():
        inps_dict['project_name'] = os.path.splitext(os.path.basename(inps_dict['template_file']))[0]

    # Input file(s) info
    fileList = ut.get_file_list(fileList, abspath=True)
    if not fileList:
        return None
    atr = readfile.read_attribute(fileList[0])
    k = atr['FILE_TYPE']
    print '--------------------------------------------'
    print 'Input file(s) is '+atr['PROCESSOR']+' '+k

    # Get output file type
    if not file_type:
        if k in ['.unw']:  file_type = 'interferograms'
        elif k in ['.cor']:  file_type = 'coherence'
        elif k in ['.int']:  file_type = 'wrapped'
        elif k in ['.byt']:  file_type = 'snap_connect_component'
        elif k in ['.msk']:  file_type = 'mask'
        elif k in ['.hgt','.dem','dem']:
            file_type = 'dem'
        elif k in ['.trans']:
            file_type = '.trans'
        else:
            file_type = k

    # Get output file name
    if not outfile:
        # output file basename
        if file_type == 'interferograms':  outfile = 'unwrapIfgram.h5'
        elif file_type == 'coherence':  outfile = 'coherence.h5'
        elif file_type == 'wrapped':  outfile = 'wrapIfgram.h5'
        elif file_type == 'snap_connect_component':  outfile = 'snapConnectComponent.h5'
        elif file_type == 'mask':  outfile = 'mask.h5'
        elif file_type == 'dem':
            if 'Y_FIRST' in atr.keys():
                outfile = 'demGeo.h5'
            else:
                outfile = 'demRadar.h5'
        elif file_type == '.trans':
            outfile = os.path.basename(fileList[0])
        else:
            warnings.warn('Un-recognized file type: '+file_type)
        
        # output directory
        if 'timeseries_dir' in inps_dict.keys():
            outdir = inps_dict['timeseries_dir']
        else:
            outdir = os.path.abspath(os.getcwd())
        outfile = outdir+'/'+outfile
    outfile = os.path.abspath(outfile)

    # Convert 
    if file_type in multi_group_hdf5_file:
        outfile = roipac2multi_group_hdf5(file_type, fileList, outfile, inps_dict)[0]

    elif file_type in single_dataset_hdf5_file:
        outfile = roipac2single_dataset_hdf5(file_type, fileList[0], outfile, inps_dict)

    elif file_type in ['.trans']:
        outfile = copy_file(fileList[0], os.path.dirname(outfile))
    else:
        raise ValueError('Un-supported file type: '+file_type)

    return outfile


def load_data_from_template(template_file, inps):
    '''Load dataset for PySAR time series using input template'''
    # Project Name
    if not inps.project_name:
        inps.project_name = os.path.splitext(os.path.basename(template_file))[0]
    
    ##------------------------------------ Read Input Path -------------------------------------##
    # Initial value
    inps.unw = None
    inps.cor = None
    inps.int = None
    inps.geomap = None
    inps.dem_radar = None
    inps.dem_geo = None

    # 1.1 Read template contents
    template_file = os.path.abspath(template_file)
    template_dict = readfile.read_template(template_file)
    for key, value in template_dict.iteritems():
        template_dict[key] = ut.check_variable_name(value)
    keyList = template_dict.keys()

    if 'pysar.unwrapFiles'    in keyList:   inps.unw       = template_dict['pysar.unwrapFiles']
    if 'pysar.corFiles'       in keyList:   inps.cor       = template_dict['pysar.corFiles']
    if 'pysar.wrapFiles'      in keyList:   inps.int       = template_dict['pysar.wrapFiles']
    if 'pysar.geomap'         in keyList:   inps.geomap    = template_dict['pysar.geomap']
    if 'pysar.dem.radarCoord' in keyList:   inps.dem_radar = template_dict['pysar.dem.radarCoord']
    if 'pysar.dem.geoCoord'   in keyList:   inps.dem_geo   = template_dict['pysar.dem.geoCoord']

    # 1.2 Auto Setting for Geodesy Lab - University of Miami 
    if pysar.miami_path and 'SCRATCHDIR' in os.environ:
        inps = auto_path_miami(inps, template_dict)

    # 1.3 get snap_connect.byt path if .unw is input
    if inps.unw:
        inps.snap_connect = inps.unw.split('.unw')[0]+'_snap_connect.byt'
    else:
        inps.snap_connect = None

    # PYSAR directory
    if not inps.timeseries_dir:
        inps.timeseries_dir = os.getcwd()
    if not os.path.isdir(inps.timeseries_dir):
        os.makedirs(inps.timeseries_dir)
    #print "PySAR working directory: "+inps.timeseries_dir
    
    # TEMPLATE file directory (to support relative path input)
    inps.template_dir = os.path.dirname(template_file)
    os.chdir(inps.template_dir)
    print 'Go to TEMPLATE directory: '+inps.template_dir
    print 'unwrapped interferograms to load: '+str(inps.unw)
    print 'wrapped   interferograms to load: '+str(inps.int)
    print 'spatial       coherences to load: '+str(inps.cor)
    print 'geomap*.trans      file to load: '+str(inps.geomap)
    print 'DEM file in radar coord to load: '+str(inps.dem_radar)
    print 'DEM file in geo   coord to load: '+str(inps.dem_geo)

    ##------------------------------------ Loading into HDF5 ---------------------------------------##
    # required - unwrapped interferograms
    inps.ifgram_file = load_file(inps.unw, vars(inps))

    # optional but recommended files - multi_group_hdf5_file
    inps.coherence_file = load_file(inps.cor, vars(inps))
    inps.wrap_ifgram_file = load_file(inps.int, vars(inps))
    if inps.snap_connect:
        inps.snap_connect_file = load_file(inps.snap_connect, vars(inps))

    # optional but recommend files - single_dataset file
    inps.geomap_file = load_file(inps.geomap, vars(inps))
    inps.dem_radar_file = load_file(inps.dem_radar, vars(inps))
    inps.dem_geo_file = load_file(inps.dem_geo, vars(inps))

    os.chdir(inps.timeseries_dir)
    print 'Go back to PYSAR directory: '+inps.timeseries_dir
    return inps


##########################  Usage  ###############################
EXAMPLE='''example:
  load_data.py  -f $SC/SanAndreasT356EnvD/PROCESS/DONE/IFG*/filt*.unw 
  load_data.py  -f $SC/SanAndreasT356EnvD/PROCESS/DONE/IFG*/filt*.unw  -o unwrapIfgram.h5
  load_data.py  -f $SC/SanAndreasT356EnvD/PROCESS/DONE/IFG*/filt*rlks.cor
  load_data.py  -f radar_4rlks.hgt  -o demRadar.h5
  load_data.py  -f srtm1.dem        -o demGeo.h5
  load_data.py  --template SanAndreasT356EnvD.tempalte
'''

TEMPLATE='''
pysar.unwrapFiles     = filt*.unw
pysar.corFiles        = filt*rlks.cor
pysar.wrapFiles       = filt*rlks.int
pysar.geomap          = geomap*.trans
pysar.dem.radarCoord  = radar*.hgt
pysar.dem.geoCoord    = srtm1.dem
'''


def cmdLineParse():
    parser = argparse.ArgumentParser(description='Load InSAR data into PySAR',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)
    parser.add_argument('--template', dest='template_file', help='template file, to get PROJECT_NAME')
    parser.add_argument('--project', dest='project_name', help='project name of dataset, used in INSARMAPS Web Viewer')

    singleFile = parser.add_argument_group('Load into single HDF5 file')
    singleFile.add_argument('-f','--file', nargs='*', help='file(s) to be loaded, processed by ROI_PAC, Gamma, DORIS or ISCE.')
    singleFile.add_argument('--file-type', dest='file_type', help='output file type, i.e.\n'+\
                            'interferograms, coherence, wrapped, dem, ...')
    singleFile.add_argument('-o','--output', dest='outfile', help='output file name')

    multiFile = parser.add_argument_group('Load whole dataset using template, i.e.',TEMPLATE)
    multiFile.add_argument('--dir', dest='timeseries_dir',\
                           help='directory for time series analysis, e.g. KujuAlosAT422F650/PYSAR')

    inps = parser.parse_args()
    # Print usage if no FILE and TEMPLATEF_FILE input
    if not inps.file and not inps.template_file:
        parser.print_usage()
        sys.exit(os.path.basename(sys.argv[0])+': error: empty FILE and TEMPLATE_FILE, at least one is needed.')
    return inps


#############################  Main Function  ################################
def main(argv):
    inps = cmdLineParse()

    if inps.file:
        # Load data into one hdf5 file
        inps.outfile = load_file(inps.file, vars(inps), inps.outfile)

    else:
        # Load the whole dataset for PySAR time series analysis, e.g. call from pysarApp.py
        inps = load_data_from_template(inps.template_file, inps)

        
    return inps.outfile

##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])


