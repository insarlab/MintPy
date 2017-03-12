#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# Yunjun, Jul 2015: Add check_num/check_file_size to .int/.cor file
# Yunjun, Jan 2017: Add auto_path_miami(), copy_roipac_file()
#                   Add load_roipac2multi_group_h5()
#                   Add r+ mode loading of multi_group hdf5 file


import os
import sys
import glob
import time
import argparse

import h5py
import numpy as np

import pysar
import pysar._readfile as readfile
import pysar._pysar_utilities as ut


############################ Sub Functions ###################################
##################################################################
def auto_path_miami(inps, template_dict={}):
    '''Auto File Path Setting for Geodesy Lab - University of Miami'''
    print 'Use auto path setting in University of Miami.'+\
          '(To turn it off, change miami_path value to False in pysar/__init__.py)'
    if not inps.timeseries_dir:
        inps.timeseries_dir = os.getenv('SCRATCHDIR')+'/'+inps.project_name+'/PYSAR'
    process_dir = os.getenv('SCRATCHDIR')+'/'+inps.project_name+'/PROCESS'
    print "PROCESS directory: "+process_dir

    if not inps.unw:   inps.unw = process_dir+'/DONE/IFGRAM*/filt_*.unw'
    if not inps.cor:   inps.cor = process_dir+'/DONE/IFGRAM*/filt_*rlks.cor'
    if not inps.int:   inps.int = process_dir+'/DONE/IFGRAM*/filt_*rlks.int'

    # Search PROCESS/GEO folder and use the first folder as master interferogram
    if not inps.geomap and not inps.dem_radar:
        try:
            master_igram_date12 = os.walk(process_dir+'/GEO').next()[1][0].split('geo_')[1]
            inps.geomap    = process_dir+'/GEO/*'+master_igram_date12+'*/geomap*.trans'
            inps.dem_radar = process_dir+'/DONE/*'+master_igram_date12+'*/radar*.hgt'
        except:
            print 'ERROR: do not find any folder in PROCESS/GEO as master interferogram'
    
    # Use DEMg/DEM option if dem_geo is not specified in pysar option
    if not inps.dem_geo and template_dict:
        if 'DEMg' in template_dict.keys():
            inps.dem_geo = template_dict['DEMg']
        elif 'DEM' in template_dict.keys():
            inps.dem_geo = template_dict['DEM']
        else:
            print 'No DEMg/DEM option found in template, continue without pysar.dem.geoCoord option.'

    return inps


########### Find Mode (most common) item in the list #############
def mode (thelist):
    counts = {}
    for item in thelist:
        counts [item] = counts.get (item, 0) + 1
    maxcount = 0
    maxitem  = None
    for k, v in counts.items ():
        if v > maxcount:
            maxitem  = k
            maxcount = v
    if maxcount == 1:                           print "All values only appear once"
    elif counts.values().count (maxcount) > 1:  print "List has multiple modes"
    else:                                       return maxitem


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
    if not mode_width and not mode_length:
        mode_width  = mode(widthList)
        mode_length = mode(lengthList)
    
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


def load_roipac2multi_group_h5(fileType, fileList, hdf5File='unwrapIfgram.h5', pysar_meta_dict=None):
    '''Load multiple ROI_PAC product into (Multi-group, one dataset and one attribute dict per group) HDF5 file.
    Inputs:
        fileType : string, i.e. interferograms, coherence, snaphu_connect_component, etc.
        fileList : list of path, ROI_PAC .unw/.cor/.int/.byt file
        hdf5File : string, file name/path of the multi-group hdf5 PySAR file
        pysar_meta_dict : dict, extra attribute dictionary 
    Outputs:
        hdf5File

    '''
    ext = os.path.splitext(fileList[0])[1]
    print '--------------------------------------------'
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
            if pysar_meta_dict:
                group.attrs['PROJECT_NAME'] = pysar_meta_dict['project_name']
        
        # End of Loop
        h5file.close()
        print 'finished writing to '+hdf5File

    return hdf5File, fileList


def roipac_nonzero_mask(unwFileList, maskFile='Mask.h5'):
    '''Generate mask for non-zero amplitude pixel of ROI_PAC .unw file list.'''
    unwFileList, width, length = check_file_size(unwFileList)
    if unwFileList:
        # Initial mask value
        if os.path.isfile(maskFile):
            maskZero, atr = readfile.read(maskFile)
            print 'update existing mask file: '+maskFile
        else:
            maskZero = np.ones([int(length), int(width)])
            atr = None
            print 'create initial mask matrix'

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
        # Attribute - existed file
        if atr:
            for key, value in atr.iteritems():
                group.attrs[key] = value

    return maskFile, unwFileList


def copy_roipac_file(targetFile, destDir):
    '''Copy ROI_PAC file and its .rsc file to destination directory.'''
    print '--------------------------------------------'
    if os.path.isfile(destDir+'/'+os.path.basename(targetFile)):
        print os.path.basename(targetFile)+'\t already exists, no need to re-load.'
    else:
        cpCmd="cp "+targetFile+" "+destDir;       print cpCmd;   os.system(cpCmd)
        cpCmd="cp "+targetFile+".rsc "+destDir;   print cpCmd;   os.system(cpCmd)


##########################  Usage  ###############################
EXAMPLE='''example:
  load_data_roipac.py  $TE/SanAndreasT356EnvD.template
  load_data_roipac.py  $TE/SanAndreasT356EnvD.template  --dir $SC/SanAndreasT356EnvD/PYSAR
'''

TEMPLATE='''template:
  pysar.unwrapFiles    = $SC/SanAndreasT356EnvD/PROCESS/DONE/IFG*/filt*.unw
  pysar.corFiles       = $SC/SanAndreasT356EnvD/PROCESS/DONE/IFG*/filt*rlks.cor
  pysar.wrapFiles      = $SC/SanAndreasT356EnvD/PROCESS/DONE/IFG*/filt*rlks.int                       #optional
  pysar.geomap         = $SC/SanAndreasT356EnvD/PROCESS/GEO/*050102-070809*/geomap*.trans
  pysar.dem.radarCoord = $SC/SanAndreasT356EnvD/PROCESS/DONE/*050102-070809*/radar*.hgt
  pysar.dem.geoCoord   = $SC/SanAndreasT356EnvD/DEM/srtm1_30m.dem                                     #optional
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Load ROI_PAC data.\n'\
                                     'Load ROI_PAC product (from process_dir to timeseries_dir) for PySAR analysis.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=TEMPLATE+'\n'+EXAMPLE)
    parser.add_argument('template_file', help='template file with path of ROI_PAC products.')
    parser.add_argument('--dir', dest='timeseries_dir', help='output directory for PySAR time series analysis.'\
                                                        'Use current directory if not assigned.')
    parser.add_argument('--nomiami', dest='auto_path_miami', action='store_false',\
                        help='Disable updating file path based on University of Miami processing structure.')

    infile_group = parser.add_argument_group('Manually input file path')
    infile_group.add_argument('--unw', nargs='*', help='ROI_PAC unwrapped interferogram files (.unw)')
    infile_group.add_argument('--cor', nargs='*', help='ROI_PAC spatial   coherence     files (.cor)')
    infile_group.add_argument('--int', nargs='*', help='ROI_PAC wrapped   interferogram files (.int)')
    infile_group.add_argument('--geomap', help='ROI_PAC geomap_*.trans file for geocoding (.trans)')
    infile_group.add_argument('--dem-radar', dest='dem_radar', help='DEM file in radar coordinate (.hgt)')
    infile_group.add_argument('--dem-geo', dest='dem_geo', help='DEM file in geo coordinate (.dem)')

    inps = parser.parse_args()
    return inps


#############################  Main Function  ################################
def main(argv):
    inps = cmdLineParse()
    #print '\n*************** Loading ROI_PAC Data into PySAR ****************'
    inps.project_name = os.path.splitext(os.path.basename(inps.template_file))[0]
    print 'project: '+inps.project_name
    
    ##### 1. Read file path
    # Priority: command line input > template > auto setting
    # Read template contents into inps Namespace
    inps.template_file = os.path.abspath(inps.template_file)
    template_dict = readfile.read_template(inps.template_file)
    for key, value in template_dict.iteritems():
        template_dict[key] = ut.check_variable_name(value)
    keyList = template_dict.keys()
    
    if not inps.unw and 'pysar.unwrapFiles'     in keyList:   inps.unw = template_dict['pysar.unwrapFiles']
    if not inps.cor and 'pysar.corFiles'        in keyList:   inps.cor = template_dict['pysar.corFiles']
    if not inps.int and 'pysar.wrapFiles'       in keyList:   inps.int = template_dict['pysar.wrapFiles']
    if not inps.geomap    and 'pysar.geomap'    in keyList:   inps.geomap    = template_dict['pysar.geomap']
    if not inps.dem_radar and 'pysar.dem.radarCoord' in keyList:   inps.dem_radar = template_dict['pysar.dem.radarCoord']
    if not inps.dem_geo   and 'pysar.dem.geoCoord'   in keyList:   inps.dem_geo   = template_dict['pysar.dem.geoCoord']

    # Auto Setting for Geodesy Lab - University of Miami 
    if pysar.miami_path and 'SCRATCHDIR' in os.environ:
        inps = auto_path_miami(inps, template_dict)

    # TIMESERIES directory for PySAR
    if not inps.timeseries_dir:
        inps.timeseries_dir = os.getcwd()
    if not os.path.isdir(inps.timeseries_dir):
        os.mkdir(inps.timeseries_dir)
    print "PySAR working directory: "+inps.timeseries_dir
    
    # TEMPLATE file directory (to support relative path input)
    inps.template_dir = os.path.dirname(inps.template_file)
    os.chdir(inps.template_dir)
    print 'Go to TEMPLATE directory: '+inps.template_dir

    # Get all file list
    inps.snap_connect = []
    if inps.unw:
        print 'unwrapped interferograms: '+str(inps.unw)
        inps.snap_connect = inps.unw.split('.unw')[0]+'_snap_connect.byt'
        inps.snap_connect = sorted(glob.glob(inps.snap_connect))
        inps.unw = sorted(glob.glob(inps.unw))
    if inps.int:
        print 'wrapped   interferograms: '+str(inps.int)
        inps.int = sorted(glob.glob(inps.int))
    if inps.cor:
        print 'coherence files: '+str(inps.cor)
        inps.cor = sorted(glob.glob(inps.cor))
    
    try:    inps.geomap = glob.glob(inps.geomap)[0]
    except: inps.geomap = None
    try:    inps.dem_radar = glob.glob(inps.dem_radar)[-1]
    except: inps.dem_radar = None
    try:    inps.dem_geo = glob.glob(inps.dem_geo)[0]
    except: inps.dem_geo = None
    print 'geomap file: '+str(inps.geomap)
    print 'DEM file in radar coord: '+str(inps.dem_radar)
    print 'DEM file in geo   coord: '+str(inps.dem_geo)

    ##### 2. Load data into hdf5 file
    inps.ifgram_file     = inps.timeseries_dir+'/unwrapIfgram.h5'
    inps.coherence_file  = inps.timeseries_dir+'/coherence.h5'
    inps.wrapIfgram_file = inps.timeseries_dir+'/wrapIfgram.h5'
    inps.snap_connect_file = inps.timeseries_dir+'/snaphuConnectComponent.h5'
    inps.mask_file = inps.timeseries_dir+'/Mask.h5'
    inps.spatial_coherence_file = inps.timeseries_dir+'/average_spatial_coherence.h5'
    
    # 2.1 multi_group_hdf5_file
    # Unwrapped Interferograms
    if inps.unw:
        unwList = load_roipac2multi_group_h5('interferograms', inps.unw, inps.ifgram_file, vars(inps))[1]
        # Update mask only when update unwrapIfgram.h5
        if unwList:
            print 'Generate mask from amplitude of interferograms'
            roipac_nonzero_mask(inps.unw, inps.mask_file)
    elif os.path.isfile(inps.ifgram_file):
        print os.path.basename(inps.ifgram_file)+' already exists, no need to re-load.'
    else:
        sys.exit('ERROR: Cannot load/find unwrapped interferograms!')

    # Optional
    if inps.snap_connect:
        load_roipac2multi_group_h5('snaphu_connect_component', inps.snap_connect, inps.snap_connect_file, vars(inps))

    # Coherence
    if inps.cor:
        cohFile,corList = load_roipac2multi_group_h5('coherence', inps.cor, inps.coherence_file, vars(inps))
        if corList:
            meanCohCmd = 'temporal_average.py '+cohFile+' '+inps.spatial_coherence_file
            print meanCohCmd
            os.system(meanCohCmd)
    elif os.path.isfile(inps.coherence_file):
        print os.path.basename(inps.coherence_file)+' already exists, no need to re-load.'
    else:
        print 'WARNING: Cannot load/find coherence.'

    # Wrapped Interferograms
    if inps.int:
        load_roipac2multi_group_h5('wrapped', inps.int, inps.wrapIfgram_file, vars(inps))
    elif os.path.isfile(inps.wrapIfgram_file):
        print os.path.basename(inps.wrapIfgram_file)+' already exists, no need to re-load.'
    else:
        print "WARNING: Cannot load/find wrapped interferograms. It's okay, continue without it ..."

    # 2.2 single dataset file
    if inps.geomap:
        copy_roipac_file(inps.geomap, inps.timeseries_dir)

    if inps.dem_radar:
        copy_roipac_file(inps.dem_radar, inps.timeseries_dir)

    if inps.dem_geo:
        copy_roipac_file(inps.dem_geo, inps.timeseries_dir)


##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
