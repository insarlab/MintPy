#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################


import os
import sys
import glob
import argparse
import warnings
import re

import h5py
import numpy as np
import shutil

import pysar
import pysar._datetime as ptime
import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._pysar_utilities as ut
from pysar._readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file


sensorList = ['ers','env','sen','s1',\
              'jers','alos','palsar','alos2','palsar2',\
              'tsx','tdx','terra','csk','cosmo',\
              'rsat','radarsat','rsat2','radarsat2','kmps5','gaofen3']

############################ Sub Functions ###################################
def project_name2sensor(projectName):
    if not projectName:
        return None
    if    re.search('Ers'    , projectName):  sensor = 'Ers'
    elif  re.search('Env'    , projectName):  sensor = 'Env'
    elif  re.search('Jers'   , projectName):  sensor = 'Jers'
    elif  re.search('Alos'   , projectName):  sensor = 'Alos'
    elif  re.search('Alos2'  , projectName):  sensor = 'Alos2' 
    elif  re.search('Tsx'    , projectName):  sensor = 'Tsx'
    elif  re.search('Tdm'    , projectName):  sensor = 'Tsx'
    elif  re.search('Csk'    , projectName):  sensor = 'Csk'
    elif  re.search('Rsat'   , projectName):  sensor = 'Rsat'
    elif  re.search('Rsat2'  , projectName):  sensor = 'Rsat2'
    elif  re.search('Sen'    , projectName):  sensor = 'Sen'
    elif  re.search('Kmps5'  , projectName):  sensor = 'Kmps5'
    elif  re.search('Gaofen3', projectName):  sensor = 'Gaofen3'
    else: print 'satellite not found';  sensor = None
    return sensor


##################################################################
def auto_path_miami(inps, template={}):
    '''Auto File Path Setting for Geodesy Lab - University of Miami'''
    print 'Use auto path setting in University of Miami.'+\
          '(To turn it off, change miami_path value to False in pysar/__init__.py)'
    # PYSAR working directory
    if not inps.timeseries_dir:
        inps.timeseries_dir = os.getenv('SCRATCHDIR')+'/'+inps.project_name+'/PYSAR'

    ##### .unw/.cor/.int files
    process_dir = os.getenv('SCRATCHDIR')+'/'+inps.project_name+'/PROCESS'
    print "PROCESS directory: "+process_dir
    if inps.insarProcessor == 'roipac':
        if not inps.unw or inps.unw == 'auto':   inps.unw = process_dir+'/DONE/IFGRAM*/filt_*.unw'
        if not inps.cor or inps.cor == 'auto':   inps.cor = process_dir+'/DONE/IFGRAM*/filt_*rlks.cor'
        #if not inps.int or inps.int == 'auto':   inps.int = process_dir+'/DONE/IFGRAM*/filt_*rlks.int'
    elif inps.insarProcessor == 'gamma':
        if not inps.unw or inps.unw == 'auto':   inps.unw = process_dir+'/DONE/IFGRAM*/diff_*rlks.unw'
        if not inps.cor or inps.cor == 'auto':   inps.cor = process_dir+'/DONE/IFGRAM*/filt_*rlks.cor'
        #if not inps.int or inps.int == 'auto':   inps.int = process_dir+'/DONE/IFGRAM*/diff_*rlks.int'

    ##### master interferogram for lookup table and DEM in radar coord
    if all(fname and fname != 'auto' for fname in [inps.lut, inps.dem_radar, inps.dem_geo]):
        return inps

    try:     m_date12 = np.loadtxt(process_dir+'/master_ifgram.txt', dtype=str).tolist()
    except:
        try: m_date12 = os.walk(process_dir+'/GEO').next()[1][0].split('geo_')[1]
        except: pass

    if not inps.lut or inps.lut == 'auto':
        try:
            if inps.insarProcessor == 'roipac':
                inps.lut = process_dir+'/GEO/*'+m_date12+'*/geomap*.trans'
            elif inps.insarProcessor == 'gamma':
                inps.lut = process_dir+'/SIM/sim_'+m_date12+'/sim_*.UTM_TO_RDC'
        except:
            warnings.warn('No master interferogram found! Can not locate mapping transformation file for geocoding!')

    if not inps.dem_radar or inps.dem_radar == 'auto':
        try:
            if inps.insarProcessor == 'roipac':
                inps.dem_radar = process_dir+'/DONE/*'+m_date12+'*/radar*.hgt'
            elif inps.insarProcessor == 'gamma':
                inps.dem_radar = process_dir+'/SIM/sim_'+m_date12+'/sim_*.hgt_sim'
        except:
            warnings.warn('No master interferogram found! Can not locate DEM in radar coord!')

    # Use DEMg/DEM option if dem_geo is not specified in pysar option
    dem_dir = os.getenv('SCRATCHDIR')+'/'+inps.project_name+'/DEM'
    if not inps.dem_geo or inps.dem_geo == 'auto':
        inps.dem_geo = []
        if os.path.isdir(dem_dir):
            inps.dem_geo = [dem_dir+'/*.dem']
        elif inps.insarProcessor == 'gamma':
            inps.dem_geo = [process_dir+'/SIM/sim_'+m_date12+'/sim_*.utm.dem']

        if   'DEMg' in template.keys():  inps.dem_geo.append(template['DEMg'])
        elif 'DEM'  in template.keys():  inps.dem_geo.append(template['DEM'])
        try:    inps.dem_geo = ut.get_file_list(inps.dem_geo)[0]
        except: inps.dem_geo = None

        if not inps.dem_geo:
            warnings.warn('Can not locate DEM in geo coord!')

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
                print '%s    width: %s    length: %s' % (os.path.basename(fileList[i]), widthList[i], lengthList[i])
                fileListOut.remove(fileList[i])
        print '\nNumber of '+ext+' left: '+str(len(fileListOut))
        print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    return fileListOut, mode_width, mode_length


def check_existed_hdf5_file(inFiles, hdf5File):
    '''Check file list with existed hdf5 file
    Return list of files that are not included in the existed readable hdf5 file.
    If all included, return None.
    '''
    outFiles = list(inFiles)
    if not outFiles:
        return outFiles

    # if previous hdf5 file existed
    if os.path.isfile(hdf5File):
        print os.path.basename(hdf5File)+'  already exists.'
        try:
            atr = readfile.read_attribute(hdf5File)
        except:
            print 'File exists but not readable, delete it.'
            rmCmd = 'rm '+hdf5File
            print rmCmd
            os.system(rmCmd)
            return outFiles

        k = atr['FILE_TYPE']
        h5 = h5py.File(hdf5File, 'r')
        h5DsetNames = sorted(h5[k].keys())
        h5.close()

        # Remove file/epoch that already existed
        for h5DsetName in h5DsetNames:
            for inFile in inFiles:
                if h5DsetName in inFile:
                    outFiles.remove(inFile)

        # Check mode length/width with existed hdf5 file
        if outFiles:
            ext = os.path.splitext(outFiles[0])[1]
            outFiles, mode_width, mode_length = check_file_size(outFiles)
            if mode_width != atr['WIDTH'] or mode_length != atr['FILE_LENGTH']:
                print 'WARNING: input files have different size than existed hdf5 file:'
                print 'Input file size: '+mode_length+', '+mode_width
                print 'HDF5  file size: '+atr['FILE_LENGTH']+', '+atr['WIDTH']
                print 'Continue WITHOUT loading '+ext+' file'
                print 'To enforse loading, change/remove existed HDF5 filename and re-run loading script'
                outFiles = None
    return outFiles


def load_multi_group_hdf5(fileType, fileList, outfile='unwrapIfgram.h5', exDict=dict()):
    '''Load multiple ROI_PAC files into HDF5 file (Multi-group, one dataset and one attribute dict per group).
    Inputs:
        fileType : string, i.e. interferograms, coherence, snaphu_connect_component, etc.
        fileList : list of path, ROI_PAC .unw/.cor/.int/.byt file
        outfile : string, file name/path of the multi-group hdf5 PySAR file
        exDict : dict, extra attribute dictionary 
    Outputs:
        outfile : output hdf5 file name
        fileList : list of string, files newly added
    '''
    ext = os.path.splitext(fileList[0])[1]
    print 'loading '+ext+' files into '+fileType+' HDF5 file ...'
    print 'number of '+ext+' input: '+str(len(fileList))

    # Check width/length mode of input files
    fileList, mode_width, mode_length = check_file_size(fileList)
    if not fileList:
        return None, None

    # Check conflict with existing hdf5 file
    fileList2 = check_existed_hdf5_file(fileList, outfile)

    # Open(Create) HDF5 file with r+/w mode based on fileList2
    if fileList2 == fileList:
        # Create and open new hdf5 file with w mode
        print 'number of '+ext+' to add: '+str(len(fileList))
        print 'open '+outfile+' with w mode'
        h5file = h5py.File(outfile, 'w')
    elif fileList2:
        # Open existed hdf5 file with r+ mode
        print 'Continue by adding the following new epochs ...'
        print 'number of '+ext+' to add: '+str(len(fileList2))
        print 'open '+outfile+' with r+ mode'
        h5file = h5py.File(outfile, 'r+')
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
            # Read data and attributes
            print 'Adding ' + file
            data, atr = readfile.read(file)

            # PySAR attributes
            atr['drop_ifgram'] = 'no'
            try:     atr['PROJECT_NAME'] = exDict['project_name']
            except:  atr['PROJECT_NAME'] = 'PYSAR'
            key = 'INSAR_PROCESSOR'
            if key not in atr.keys():
                try:  atr[key] = exDict['insarProcessor']
                except:  pass
            key = 'PLATFORM'
            if ((key not in atr.keys() or not any(re.search(i, atr[key].lower()) for i in sensorList))\
                and exDict['PLATFORM']):
                atr[key] = exDict['PLATFORM']

            # Write dataset
            group = gg.create_group(os.path.basename(file))
            dset = group.create_dataset(os.path.basename(file), data=data, compression='gzip')

            # Write attributes
            for key, value in atr.iteritems():
                group.attrs[key] = str(value)

        # End of Loop
        h5file.close()
        print 'finished writing to '+outfile

    return outfile, fileList


def load_geometry_hdf5(fileType, fileList, outfile=None, exDict=dict()):
    '''Load multiple geometry files into hdf5 file: geometryGeo.h5 or geometryRadar.h5.
    File structure:
        /geometry.attrs
        /geometry/latitude          #for geometryRadar.h5 only, from ISCE/Doris lookup table
        /geometry/longitude         #for geometryRadar.h5 only, from ISCE/Doris lookup table
        /geometry/rangeCoord        #for geometryGeo.h5 only, from ROI_PAC/Gamma lookup table
        /geometry/azimuthCoord      #for geometryGeo.h5 only, from ROI_PAC/Gamma lookup table
        /geometry/height
        /geometry/incidenceAngle
        /geometry/headingAngle
        /geometry/slantRangeDistance
        /geometry/shadowMask
        /geometry/waterMask
    '''
    ext = os.path.splitext(fileList[0])[1]
    atr = readfile.read_attribute(fileList[0])
    if not outfile:
        if 'Y_FIRST' in atr.keys():
            outfile = 'geometryGeo.h5'
        else:
            outfile = 'geometryRadar.h5'
        # output directory
        if 'timeseries_dir' in exDict.keys() and exDict['timeseries_dir']:
            outdir = exDict['timeseries_dir']
        else:
            outdir = os.path.abspath(os.getcwd())
        outfile = os.path.join(outdir, outfile)
    outfile = os.path.abspath(outfile)

    #####Check overlap with existing hdf5 file
    h5dnameList = []
    if os.path.isfile(outfile):
        print os.path.basename(outfile)+'  already exists.'
        try:
            atr = readfile.read_attribute(outfile)
        except:
            print 'File exists but not readable, delete it.'
            rmCmd = 'rm '+outfile
            print rmCmd
            os.system(rmCmd)
        h5 = h5py.File(outfile, 'r')
        h5dnameList = sorted(h5['geometry'].keys())
        h5.close()

    dnameList = []
    fileList2 = list(fileList)
    for fname in fileList2:
        fbase = os.path.basename(fname).lower()
        if ((fbase.startswith('lat') and 'latitude' in h5dnameList) or\
            (fbase.startswith('lon') and 'longitude' in h5dnameList) or\
            (fbase.startswith('los') and 'incidenceAngle' in h5dnameList) or\
            (fbase.startswith('shadowmask') and 'shadowMask' in h5dnameList) or\
            (fbase.startswith('watermask') and 'waterMask' in h5dnameList) or\
            (fbase.startswith('incidenceang') and 'incidenceAngle' in h5dnameList) or\
            (fbase.endswith(('.trans','.utm_to_rdc')) and 'rangeCoord' in h5dnameList) or\
            ((fbase.startswith(('hgt','dem')) or fbase.endswith(('.hgt','.dem','wgs84'))) and 'height' in h5dnameList) or\
            (os.path.abspath(fname) == outfile) or\
            #(fbase.startswith('geometry') and any(i in h5dnameList for i in ['rangeCoord','longitude'])) or \
            (fbase.startswith('rangedist') and 'slantRangeDistance' in h5dnameList)):
            fileList.remove(fname)

    # Loop - Writing files into hdf5 file
    if fileList:
        print 'number of '+ext+' to add: '+str(len(fileList))
        ##Open HDF5 file
        if os.path.isfile(outfile):
            print 'open '+outfile+' with r+ mode'
            h5 = h5py.File(outfile, 'r+')
        else:
            print 'open '+outfile+' with w mode'
            h5 = h5py.File(outfile, 'w')            

        ##top level group
        if fileType not in h5.keys():
            group = h5.create_group(fileType)
            print 'create group: '+fileType
        else:
            group = h5[fileType]
            print 'open group: '+fileType
        ##datasets
        for fname in fileList:
            fbase = os.path.basename(fname).lower()
            print 'Add '+fname
            if fbase.startswith('lat'):
                data, atr = readfile.read(fname)
                dset = group.create_dataset('latitude', data=data, compression='gzip')

            elif fbase.startswith('lon'):
                data, atr = readfile.read(fname)
                dset = group.create_dataset('longitude', data=data, compression='gzip')

            elif fbase.startswith('los'):
                d0, d1, atr = readfile.read(fname)
                dset = group.create_dataset('incidenceAngle', data=d0, compression='gzip')
                dset = group.create_dataset('headingAngle', data=d1, compression='gzip')

            elif 'shadowmask' in fbase:
                data, atr = readfile.read(fname)
                dset = group.create_dataset('shadowMask', data=data, compression='gzip')

            elif 'watermask' in fbase:
                data, atr = readfile.read(fname)
                dset = group.create_dataset('waterMask', data=data, compression='gzip')

            elif 'incidenceang' in fbase:
                data, atr = readfile.read(fname)
                dset = group.create_dataset('incidenceAngle', data=data, compression='gzip')

            elif 'rangedist' in fbase:
                data, atr = readfile.read(fname)
                dset = group.create_dataset('slantRangeDistance', data=data, compression='gzip')

            elif fbase.endswith(('.trans','.utm_to_rdc')):
                d0, d1, atr = readfile.read(fname)
                dset = group.create_dataset('rangeCoord', data=d0, compression='gzip')
                dset = group.create_dataset('azimuthCoord', data=d1, compression='gzip')

            elif fbase.startswith(('hgt','dem')) or fbase.endswith(('.hgt','.dem','wgs84')):
                data, atr = readfile.read(fname)
                dset = group.create_dataset('height', data=data, compression='gzip')

            else:
                print 'Un-recognized file type: '+fbase

            # PySAR attributes
            try:     atr['PROJECT_NAME'] = exDict['project_name']
            except:  atr['PROJECT_NAME'] = 'PYSAR'
            key = 'INSAR_PROCESSOR'
            if key not in atr.keys():
                try:  atr[key] = exDict['insarProcessor']
                except:  pass
            # Write attributes
            for key,value in atr.iteritems():
                if key not in group.attrs.keys():
                    group.attrs[key] = str(value)
        h5.close()
    else:
        print 'All input '+ext+' are included, no need to re-load.'
        fileList = None
    return outfile


def load_single_dataset_hdf5(file_type, infile, outfile=None, exDict=dict()):
    '''Convert ROI_PAC .dem / .hgt file to hdf5 file
    Based on load_dem.py written by Emre Havazli
    Inputs:
        file_type : string, group name of hdf5 file, i.e. dem, mask
        infile    : string, input ROI_PAC file name
        outfile   : string, output hdf5 file name
        exDict : dict, extra attributes to output file
    Output:
        outfile   : string, output hdf5 file name
    '''
    atr = readfile.read_attribute(infile)

    if ut.update_file(outfile, infile):
        if (os.path.dirname(infile) == os.path.dirname(outfile) and \
            os.path.splitext(infile)[1] == os.path.splitext(outfile)[1]):
            print infile+' already in working directory with recommended format, no need to re-load.'
            outfile = infile

        else:
            # Read input file
            print 'loading file: '+infile
            data = readfile.read(infile)[0]

            # Write output file - data
            print 'writing >>> '+outfile
            h5 = h5py.File(outfile, 'w')
            group = h5.create_group(file_type)
            dset = group.create_dataset(file_type, data=data, compression='gzip')

            # Write output file - attributes
            for key, value in atr.iteritems():
                group.attrs[key] = value
            try: group.attrs['PROJECT_NAME'] = exDict['project_name']
            except: pass
            key = 'INSAR_PROCESSOR'
            if key not in atr.keys():
                try:  atr[key] = exDict['insarProcessor']
                except:  pass
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
        fileList  - string / list of string, path of files to load
        inps_dict - dict, including the following attributes
                    PROJECT_NAME   : KujuAlosAT422F650  (extra attribute dictionary to add to output file)
                    sensor         : (optional)
                    timeseries_dir : directory of time series analysis, e.g. KujuAlosAT422F650/PYSAR
                    insarProcessor: InSAR processor, roipac, isce, gamma, doris
        outfile   - string, output file name
        file_type - string, group name for output HDF5 file, interferograms, coherence, dem, etc.
    Output:
        outfile - string, output file name
    Example:
        unwrapIfgram.h5 = load_file('filt*.unw', inps_dict=vars(inps))
    '''
    # Get project_name from input template file
    if not 'project_name' in inps_dict.keys() and 'template_file' in inps_dict.keys():
        template_filename_list = [os.path.basename(i) for i in inps_dict['template_file']]
        try:  template_filename_list.remove('pysarApp_template.txt')
        except:  pass
        if template_filename_list:
            inps_dict['project_name'] = os.path.splitext(template_filename_list[0])[0]

    #Sensor
    inps_dict['PLATFORM'] = project_name2sensor(inps_dict['project_name'])        

    # Input file(s) info
    fileList = ut.get_file_list(fileList, abspath=True)
    if not fileList:
        return None

    ##### Prepare attributes file
    processor = inps_dict['insarProcessor']
    print '--------------------------------------------'
    print 'preparing attributes files using prep_%s.py ...' % processor
    # prepare multiple files input for cmd calling
    files_input = ''
    for x in fileList:
        files_input += x+' '
    # call prepare_*.py
    if   processor == 'gamma' :  prepCmd = 'prep_gamma.py ' +files_input;   os.system(prepCmd)
    elif processor == 'roipac':  prepCmd = 'prep_roipac.py '+files_input;   os.system(prepCmd)
    elif processor == 'isce'  :  prepCmd = 'prep_isce.py '+files_input;     #os.system(prepCmd)
    else:
        print 'Un-supported InSAR processor: '+processor
        print 'Skip preparing attributes files'

    print '----------------------------'
    print 'loading files ...'
    atr = readfile.read_attribute(fileList[0])
    k = atr['FILE_TYPE']
    print 'Input file(s) is '+atr['PROCESSOR']+' '+k

    # Get output file type
    if not file_type:
        if k in ['.unw']:  file_type = 'interferograms'
        elif k in ['.cor']:  file_type = 'coherence'
        elif k in ['.int']:  file_type = 'wrapped'
        elif k in ['.byt']:  file_type = 'snaphu_connect_component'
        elif k in ['.msk']:  file_type = 'mask'
        elif k in ['.hgt','.dem','dem','.hgt_sim']:
            file_type = 'dem'
        elif k in ['.trans','.utm_to_rdc','geometry']:
            file_type = 'geometry'
        else:
            file_type = k

    # Get output file name
    if not outfile:
        # output file basename
        if file_type == 'interferograms':  outfile = 'unwrapIfgram.h5'
        elif file_type == 'coherence':  outfile = 'coherence.h5'
        elif file_type == 'wrapped':  outfile = 'wrapIfgram.h5'
        elif file_type == 'snaphu_connect_component':  outfile = 'snaphuConnectComponent.h5'
        elif file_type == 'mask':  outfile = 'mask.h5'
        elif file_type == 'dem':
            if 'Y_FIRST' in atr.keys():
                outfile = 'demGeo.h5'
            else:
                outfile = 'demRadar.h5'

        # output directory
        if 'timeseries_dir' in inps_dict.keys() and inps_dict['timeseries_dir']:
            outdir = inps_dict['timeseries_dir']
        else:
            outdir = os.path.abspath(os.getcwd())
        if outfile:
            outfile = outdir+'/'+outfile
    if outfile:
        outfile = os.path.abspath(outfile)

    # Convert
    if file_type in multi_group_hdf5_file:
        outfile = load_multi_group_hdf5(file_type, fileList, outfile=outfile, exDict=inps_dict)[0]

    elif file_type in single_dataset_hdf5_file:
        outfile = load_single_dataset_hdf5(file_type, fileList[-1], outfile=outfile, exDict=inps_dict)

    elif file_type in ['geometry','.trans','.utm_to_rdc','.UTM_TO_RDC']:
        outfile = load_geometry_hdf5(file_type, fileList, outfile=outfile, exDict=inps_dict)
    else:
        warnings.warn('Un-supported file type: '+file_type)

    return outfile


def load_data_from_template(inps):
    '''Load dataset for PySAR time series using input template'''
    ##------------------------------------ Read Input Path -------------------------------------##
    # Initial value
    inps.unw = None
    inps.cor = None
    #inps.int = None
    inps.lut = None
    inps.dem_radar = None
    inps.dem_geo = None

    # 1.1 Read template contents (support multiple template files input)
    inps.template_file = [os.path.abspath(i) for i in inps.template_file]
    # Move default template file pysarApp_template.txt to the end of list, so that it has highest priority
    default_template_file = [i for i in inps.template_file if 'pysarApp' in i]
    if default_template_file:
        inps.template_file.remove(default_template_file[0])
        inps.template_file.append(default_template_file[0])
    template = dict()
    # Read file by file
    for File in inps.template_file:
        temp_dict = readfile.read_template(File)
        for key, value in temp_dict.iteritems():
            temp_dict[key] = ut.check_variable_name(value)
        template.update(temp_dict)
    keyList = template.keys()

    # Project Name
    if not inps.project_name:
        inps.template_filename = [os.path.basename(i) for i in inps.template_file]
        try:  inps.template_filename.remove('pysarApp_template.txt')
        except:  pass
        if inps.template_filename:
            inps.project_name = os.path.splitext(inps.template_filename[0])[0]

    for key in ['processor','processing_software','unavco.processing_software','pysar.insarProcessor']:
        if key in keyList:
            value = template[key]
            if value == 'auto':
                inps.insarProcessor = 'roipac'
            else:
                inps.insarProcessor = value

    print '--------------------------------------------'
    print 'InSAR processing software: '+inps.insarProcessor
    if 'pysar.unwrapFiles'        in keyList:   inps.unw       = template['pysar.unwrapFiles']
    if 'pysar.corFiles'           in keyList:   inps.cor       = template['pysar.corFiles']
    if 'pysar.lookupFile'         in keyList:   inps.lut       = template['pysar.lookupFile']
    if 'pysar.demFile.radarCoord' in keyList:   inps.dem_radar = template['pysar.demFile.radarCoord']
    if 'pysar.demFile.geoCoord'   in keyList:   inps.dem_geo   = template['pysar.demFile.geoCoord']

    # Check existed single dataset files
    inps_tmp = argparse.Namespace()
    inps_tmp = ut.check_loaded_dataset(inps.timeseries_dir, inps_tmp, print_msg=False)
    if (not inps.lut       or inps.lut       == 'auto') and inps_tmp.lookup_file   :  inps.lut       = inps_tmp.lookup_file
    if (not inps.dem_radar or inps.dem_radar == 'auto') and inps_tmp.dem_radar_file:  inps.dem_radar = inps_tmp.dem_radar_file
    if (not inps.dem_geo   or inps.dem_geo   == 'auto') and inps_tmp.dem_geo_file  :  inps.dem_geo   = inps_tmp.dem_geo_file

    # 1.2 Auto Setting for Geodesy Lab - University of Miami 
    if pysar.miami_path and 'SCRATCHDIR' in os.environ and inps.project_name:
        inps = auto_path_miami(inps, template)

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
    inps.template_dir = os.path.dirname(inps.template_file[-1])
    os.chdir(inps.template_dir)
    print 'Go to TEMPLATE directory: '+inps.template_dir
    print 'unwrapped interferograms to load: '+str(inps.unw)
    #print 'wrapped   interferograms to load: '+str(inps.int)
    print 'spatial coherence  files to load: '+str(inps.cor)
    print 'lookup table        file to load: '+str(inps.lut)
    print 'DEM file in radar  coord to load: '+str(inps.dem_radar)
    print 'DEM file in geo    coord to load: '+str(inps.dem_geo)

    ##------------------------------------ Loading into HDF5 ---------------------------------------##
    # required - unwrapped interferograms
    inps.ifgram_file      = load_file(inps.unw, vars(inps), file_type='interferograms')
    inps.coherence_file   = load_file(inps.cor, vars(inps), file_type='coherence')
    #inps.wrap_ifgram_file = load_file(inps.int, vars(inps), file_type='wrapped')
    if inps.snap_connect:
        inps.snap_connect_file = load_file(inps.snap_connect, vars(inps))

    # optional but recommend files - single_dataset file
    inps.lookup_file    = load_file(inps.lut,       vars(inps), file_type='geometry')
    inps.dem_radar_file = load_file(inps.dem_radar, vars(inps), file_type='dem')
    inps.dem_geo_file   = load_file(inps.dem_geo,   vars(inps), file_type='dem')

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
  load_data.py  --template pysarApp_template.txt SanAndreasT356EnvD.tempalte

  load_data.py -f demRadar.h5 ./merged/geom_master/*.rdr rangeDistance.h5       --file-type geometry
  load_data.py -f sim*.UTM_TO_RDC demGeo.h5 geomap*.trans geo_incidenceAngle.h5 --file-type geometry
'''

TEMPLATE='''
## 1. Load Data (--load to exit after this step)
## recommended input files for data in radar coordinates:
##     pysar.insarProcessor     = InSAR processor
##     pysar.unwrapFiles        = 'path of all unwrapped interferograms'
##     pysar.corFiles           = 'path of all coherence files'
##     pysar.lookupFile         = 'path of lookup table / mapping transformation file'
##     pysar.demFile.geoCoord   = 'path of DEM in geo   coordinates'
##     pysar.demFile.radarCoord = 'path of DEM in radar coordinates'
## recommended input files for data in geo coordinates:
##     pysar.insarProcessor
##     pysar.unwrapFiles 
##     pysar.corFiles    
##     pysar.dem.geoCoord
## auto - automatic path pattern for Univ of Miami file structure, which are:
##     pysar.insarProcessor     = roipac
##     pysar.unwrapFiles        = $SCRATCHDIR/$PROJECT_NAME/DONE/IFGRAM*/filt_*.unw
##     pysar.corFiles           = $SCRATCHDIR/$PROJECT_NAME/DONE/IFGRAM*/filt_*rlks.cor
##     pysar.lookupFile         = $SCRATCHDIR/$PROJECT_NAME/GEO/*master_date12*/geomap*.trans
##     pysar.demFile.geoCoord   = $SCRATCHDIR/$PROJECT_NAME/DEM/*.dem
##     pysar.demFile.radarCoord = $SCRATCHDIR/$PROJECT_NAME/DONE/*master_date12*/radar*.hgt
pysar.insarProcessor     = auto  #[roipac,        gamma,           isce], auto for roipac
pysar.unwrapFiles        = auto  #[filt*rlks.unw, diff_*rlks.unw,  filt*.unw]
pysar.corFiles           = auto  #[filt*rlks.cor, filt_*rlks.cor,  filt*.cor]
pysar.lookupFile         = auto  #[geomap*.trans, sim*.UTM_TO_RDC, l*.rdr]
pysar.demFile.radarCoord = auto  #[radar*.hgt,    sim*.hgt_sim,    hgt.rdr]
pysar.demFile.geoCoord   = auto  #[*.dem,         sim*.utm.dem,    None]
'''


def cmdLineParse():
    parser = argparse.ArgumentParser(description='Load InSAR data into PySAR',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('--template','-t', dest='template_file', nargs='*', help='template file, to get PROJECT_NAME')
    parser.add_argument('--project', dest='project_name', help='project name of dataset, used in INSARMAPS Web Viewer')
    parser.add_argument('--processor','-p', dest='insarProcessor',\
                        default='roipac', choices={'roipac','gamma','isce','doris','gmtsar'},\
                        help='InSAR processor/software of the file')

    singleFile = parser.add_argument_group('Load into single HDF5 file')
    singleFile.add_argument('-f','--file', nargs='*', help='file(s) to be loaded, processed by ROI_PAC, Gamma, DORIS or ISCE.')
    singleFile.add_argument('--file-type', dest='file_type', help='output file type, i.e.\n'+\
                            'interferograms, coherence, wrapped, dem, ...')
    singleFile.add_argument('-o','--output', dest='outfile', help='output file name')

    multiFile = parser.add_argument_group('Load whole dataset using template, i.e.',TEMPLATE)
    multiFile.add_argument('--dir', dest='timeseries_dir',\
                           help='directory for time series analysis, e.g. KujuAlosAT422F650/PYSAR\n'+\
                                'use current directory by default if pysar.miami_path is False')

    inps = parser.parse_args()
    # Print usage if no FILE and TEMPLATEF_FILE input
    if not inps.file and not inps.template_file:
        parser.print_usage()
        sys.exit(os.path.basename(sys.argv[0])+': error: empty FILE and TEMPLATE_FILE, at least one is needed.')
    return inps


#############################  Main Function  ################################
def main(argv):
    inps = cmdLineParse()

    # Load data into one hdf5 file
    if inps.file:
        inps.outfile = load_file(inps.file, vars(inps), outfile=inps.outfile, file_type=inps.file_type)

    # Load the whole dataset for PySAR time series analysis, e.g. call from pysarApp.py
    else:
        inps = load_data_from_template(inps)

    return inps.outfile

##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])


