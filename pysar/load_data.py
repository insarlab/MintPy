#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# Yunjun, Jul 2015: Add check_num/check_file_size to .int/.cor file
#

import os
import sys
import glob
import time
import argparse

import h5py
import numpy as np

import pysar._readfile as readfile
from pysar._pysar_utilities import check_variable_name



############################ Sub Functions ###################################
##################################################################
def auto_path_miami(inps):
    '''Auto File Path Setting for Geodesy Lab - University of Miami'''
    if not inps.tssar_dir:
        inps.tssar_dir = os.getenv('SCRATCHDIR')+'/'+inps.project_name+'/TSSAR'
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
        rsc = readfile.read_roipac_rsc(file+'.rsc')
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


def check_existed_hdf5_file(roipacFileList, pysarFile):
    '''Check file list with existed hdf5 file'''
    # If input file list is empty
    if not roipacFileList:
        return roipacFileList
    
    # if previous hdf5 file existed
    if os.path.isfile(pysarFile):
        print os.path.basename(pysarFile)+'  already exists.'
        h5 = h5py.File(pysarFile, 'r')
        epochList = sorted(h5[k].keys())
        atr = h5[k][epochList[0]].attrs
        h5.close()
        
        # Remove file/epoch that already existed
        roipacFileList = list(set(roipacFileList) - set(epochList))

        # Check mode length/width with existed hdf5 file
        ext = os.path.splitext(roipacFileList[0])[1]
        roipacFileList, mode_width, mode_length = check_file_size(roipacFileList)
        if mode_width != atr['WIDTH'] or mode_length != atr['FILE_LENGTH']:
            print 'WARNING: input ROI_PAC files have different size than existed hdf5 file:'
            print 'ROI_PAC file size: '+mode_length+', '+mode_width
            print 'HDF5    file size: '+atr['FILE_LENGTH']+', '+atr['WIDTH']
            print 'Continue WITHOUT loading '+ext+' file'
            print 'To enforse loading, change/remove existed HDF5 filename and re-run loading script'
            roipacFileList = None
    
    return roipacFileList


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
    print 'loading ROI_PAC'+ext+' files into '+fileType+' HDF5 file ...'
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
        
        for igram in fileList:
            print 'Adding ' + igram
            amp,unw,unwrsc = readfile.read_float32(igram)
            
            # Dataset
            group = gg.create_group(os.path.basename(igram))
            dset = group.create_dataset(os.path.basename(igram), data=unw, compression='gzip')
            
            # Attribute - *.unw.rsc
            for key,value in unwrsc.iteritems():
                group.attrs[key] = value
            # Attribute - *baseline.rsc
            d1, d2 = unwrsc['DATE12'].split('-')
            baseline_file = os.path.dirname(igram)+'/'+d1+'_'+d2+'_baseline.rsc'
            baseline_rsc = readfile.read_roipac_rsc(baseline_file)
            for key,value in baseline_rsc.iteritems():
                group.attrs[key] = value
            # Attribute - PySAR
            if pysar_meta_dict:
                group.attrs['PROJECT_NAME'] = pysar_meta_dict['project_name']
        
        # End of Loop
        h5file.close()
        print 'finished writing to '+hdf5File

    return fileList, hdf5File


##########################  Usage  ###############################
def usage():
    print '''
************************************************************************
   loading the processed data for PySAR:
       interferograms (unwrapped and wrapped)
       coherence files
       geomap.trans file
       DEM (radar and geo coordinate)
   
   Usage: load_data.py TEMPLATEFILE  [inDir outDir]

   Example:
       load_data.py $TE/SanAndreasT356EnvD.template
       load_data.py $TE/SanAndreasT356EnvD.template $SC/PROCESS/SanAndreasT356EnvD $SC/TSSAR/SanAndreasT356EnvD

************************************************************************
    '''
    return

EXAMPLE='''example:
  load_data_roipac.py  $TE/SanAndreasT356EnvD.template
  load_data_roipac.py  $TE/SanAndreasT356EnvD.template  $SC/SanAndreasT356EnvD/PROCESS  $SC/SanAndreasT356EnvD/TSSAR
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
                                     'Load ROI_PAC product (from process_dir to tssar_dir) for PySAR analysis.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE+'\n'+TEMPLATE)
    parser.add_argument('template_file', help='template file with path of ROI_PAC products.')
    parser.add_argument('--dir', dest='tssar_dir', help='output directory for PySAR time series analysis.'\
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
    print '\n*************** Loading ROI_PAC Data into PySAR ****************'
    inps.project_name = os.path.basename(inps.template_file).partition('.')[0]
    print 'project: '+inps.project_name
    
    ##### 1. Read file path
    # Priority: command line input > template > auto setting
    # Read template contents into inps Namespace
    template_dict = readfile.read_template(inps.template_file)
    for key, value in template_dict.iteritems():
        template_dict[key] = check_variable_name(value)
    keyList = template_dict.keys()
    
    if not inps.unw and 'pysar.unwrapFiles'     in keyList:   inps.unw = template_dict['pysar.unwrapFiles']
    if not inps.cor and 'pysar.corFiles'        in keyList:   inps.cor = template_dict['pysar.corFiles']
    if not inps.int and 'pysar.wrapFiles'       in keyList:   inps.int = template_dict['pysar.wrapFiles']
    if not inps.geomap    and 'pysar.geomap'    in keyList:   inps.geomap    = template_dict['pysar.geomap']
    if not inps.dem_radar and 'pysar.dem.radarCoord' in keyList:   inps.dem_radar = template_dict['pysar.dem.radarCoord']
    if not inps.dem_geo   and 'pysar.dem.geoCoord'   in keyList:   inps.dem_geo   = template_dict['pysar.dem.geoCoord']

    # Auto Setting for Geodesy Lab - University of Miami 
    if inps.auto_path_miami:
        inps = auto_path_miami(inps)

    # Working directory for PySAR
    if not inps.tssar_dir:
        inps.tssar_dir = os.getcwd()
    if not os.path.isdir(inps.tssar_dir):
        os.mkdir(inps.tssar_dir)
    print "work    directory: "+inps.tssar_dir

    # display message
    print 'unwrapped interferograms: '+inps.unw
    print 'wrapped   interferograms: '+inps.int
    print 'coherence files         : '+inps.cor
    print 'geomap    file          : '+inps.geomap
    print 'DEM file in radar coord : '+inps.dem_radar
    print 'DEM file in geo   coord : '+inps.dem_geo
    
    # Get all file list
    if inps.unw:  inps.unw = sorted(glob.glob(inps.unw))
    if inps.cor:  inps.cor = sorted(glob.glob(inps.cor))
    if inps.int:  inps.int = sorted(glob.glob(inps.int))
    if inps.geomap:     inps.geomap    = glob.glob(inps.geomap)[0]
    if inps.dem_radar:  inps.dem_radar = glob.glob(inps.dem_radar)[0]
    if inps.dem_geo:    inps.dem_geo   = glob.glob(inps.dem_geo)[0]
    import pdb; pdb.set_trace()
    inps.snap_connect = []
        
    ##### 2. Load data into hdf5 file
    # 2.1 Unwrapped Interferograms - unwrapIfgram.h5
    if inps.unw:
        load_roipac2multi_group_h5('interferograms', inps.unw, inps.tssar_dir+'/unwrapIfgram.h5', vars(inps))
        
        
        import pdb; pdb.set_trace()
    
    
    
    
    
    
    if inps.unw:
        print '--------------------------------------------'
        print 'loading unwrapped interferograms ...'
        k = 'interferograms'
        outfile = inps.tssar_dir+'/unwrapIfgram.h5'
        maskfile = inps.tssar_dir+'/Mask.h5'
        print 'number of '+k+' input: '+str(len(inps.unw))

        # Check width/length mode of input files
        inps.unw, mode_width, mode_length = check_file_size(k, inps.unw)
        # Check conflict with existing hdf5 file and Open hdf5 file
        inps.unw, h5file = check_with_existed_hdf5_file(inps.unw, outfile)

        # Loop - Writing ROI_PAC files into hdf5 file
        if inps.unw:
            # Mask
            try:
                MaskZero = readfile.read(maskfile)[0]
                print 'updating existing '+maskfile+' ...'
            except:
                MaskZero = np.ones([int(mode_length), int(mode_width)])
                print 'creating new '+maskfile+' ...'
            
            # Unwraped Interferograms
            if not k in h5file.keys():
                gg = h5file.create_group(k)
            else:
                gg = h5file[k]
            for igram in inps.unw:
                print 'Adding ' + igram
                amp,unw,unwrsc = readfile.read_float32(igram)
                
                # Dataset
                group = gg.create_group(os.path.basename(igram))
                dset = group.create_dataset(os.path.basename(igram), data=unw, compression='gzip')
                
                # Attribute - *.unw.rsc
                for key,value in unwrsc.iteritems():
                    group.attrs[key] = value
                # Attribute - *baseline.rsc
                d1, d2 = unwrsc['DATE12'].split('-')
                baseline_file = os.path.dirname(igram)+'/'+d1+'_'+d2+'_baseline.rsc'
                baseline_rsc = readfile.read_roipac_rsc(baseline_file)
                for key,value in baseline_rsc.iteritems():
                    group.attrs[key] = value
                # Attribute - PySAR
                group.attrs['PROJECT_NAME'] = project_name
                group.attrs['UNIT'] = 'radian'
                
                # Mask
                MaskZero *= amp
            # End of Loop
            f.close()
            print 'finished writing to '+outfile

            # writing Mask file
            Mask = np.ones([int(mode_length), int(mode_width)])
            Mask[MaskZero==0] = 0
            h5mask = h5py.File(maskfile,'w')
            group = h5mask.create_group('mask')
            dset = group.create_dataset('mask', data=Mask, compression='gzip')
            for key,value in unwrsc.iteritems():
                group.attrs[key] = value
            h5mask.close()

    ########################################################################
    ############################# Coherence ################################
    try:
        if os.path.isfile(inps.tssar_dir+'/Coherence.h5'):
            print '\nCoherence.h5'+ '  already exists.'
            sys.exit(1)
        corList = glob.glob(corPath)
        corList = sorted(corList)
        k = 'coherence'
        check_number(k,optionName[k],corList)   # number check 
        print 'loading coherence files ...'
        corList,mode_width,mode_length = check_file_size(k,corList)     # size check
        corList = sorted(corList)
    
        h5file = inps.tssar_dir+'/Coherence.h5'
        fcor = h5py.File(h5file,'w')
        gg = fcor.create_group('coherence')
        meanCoherence=np.zeros([int(mode_length),int(mode_width)])
        for cor in corList:
            if not os.path.basename(cor) in fcor:
                print 'Adding ' + cor
                group = gg.create_group(os.path.basename(cor))
                amp,unw,unwrsc = readfile.read_float32(cor)
    
                meanCoherence += unw
                dset = group.create_dataset(os.path.basename(cor), data=unw, compression='gzip')
                for key,value in unwrsc.iteritems():    group.attrs[key] = value
    
                d1,d2=unwrsc['DATE12'].split('-')
                baseline_file=os.path.dirname(cor)+'/'+d1+'_'+d2+'_baseline.rsc'
                baseline=readfile.read_roipac_rsc(baseline_file)
                for key,value in baseline.iteritems():   group.attrs[key] = value
                group.attrs['PROJECT_NAME'] = project_name
                group.attrs['UNIT']         = '1'
            else:
                print os.path.basename(h5file) + " already contains " + os.path.basename(cor)
        #fcor.close()
    
        ########### mean coherence file ###############
        meanCoherence=meanCoherence/(len(corList))
        print 'writing meanCoherence group to the coherence h5 file'
        gc = fcor.create_group('meanCoherence')
        dset = gc.create_dataset('meanCoherence', data=meanCoherence, compression='gzip')
    
        print 'writing average_spatial_coherence.h5\n'
        h5file_CorMean = inps.tssar_dir+'/average_spatial_coherence.h5'
        fcor_mean = h5py.File(h5file_CorMean,'w')
        group=fcor_mean.create_group('mask')
        dset = group.create_dataset(os.path.basename('mask'), data=meanCoherence, compression='gzip')
        for key,value in unwrsc.iteritems():    group.attrs[key] = value
        fcor_mean.close()
    
        fcor.close()
  
    except:
        print 'No correlation file is loaded.\n'


    ##############################################################################
    ########################## Wrapped Interferograms ############################

    try:
        if os.path.isfile(inps.tssar_dir+'/Wrapped.h5'):
            print '\nWrapped.h5'+ '  already exists.'
            sys.exit(1)
        wrapList = glob.glob(wrapPath)
        wrapList = sorted(wrapList)
        k = 'wrapped'
        check_number(k,optionName[k],wrapList)   # number check 
        print 'loading wrapped phase ...'
        wrapList,mode_width,mode_length = check_file_size(k,wrapList)     # size check
        wrapList = sorted(wrapList)
    
        h5file = inps.tssar_dir+'/Wrapped.h5'
        fw = h5py.File(h5file,'w')
        gg = fw.create_group('wrapped')
        for wrap in wrapList:
            if not os.path.basename(wrap) in fw:
                print 'Adding ' + wrap
                group = gg.create_group(os.path.basename(wrap))
                amp,unw,unwrsc = readfile.read_complex_float32(wrap)
    
                dset = group.create_dataset(os.path.basename(wrap), data=unw, compression='gzip')
                for key,value in unwrsc.iteritems():   group.attrs[key] = value
    
                d1,d2=unwrsc['DATE12'].split('-')
                baseline_file=os.path.dirname(wrap)+'/'+d1+'_'+d2+'_baseline.rsc'
                baseline=readfile.read_roipac_rsc(baseline_file)
                for key,value in baseline.iteritems():    group.attrs[key] = value
                group.attrs['PROJECT_NAME'] = project_name
                group.attrs['UNIT']         = 'radian'
            else:
                print os.path.basename(h5file) + " already contains " + os.path.basename(wrap)
        fw.close()
        print 'Writed '+str(len(wrapList))+' wrapped interferograms to '+h5file+'\n'
  
    except:
        print 'No wrapped interferogram is loaded.\n'


    ##############################################################################
    ################################# geomap.trans ###############################

    try:
        geomapPath = inps.tssar_dir+'/geomap*.trans'
        geomapList = glob.glob(geomapPath)
        if len(geomapList)>0:
            print '\ngeomap*.trans'+ '  already exists.'
            sys.exit(1)
    
        geomapPath=template_dict['pysar.geomap']
        geomapPath=check_variable_name(geomapPath)
        geomapList=glob.glob(geomapPath)
    
        cpCmd="cp " + geomapList[0] + " " + inps.tssar_dir
        print cpCmd
        os.system(cpCmd)
        cpCmd="cp " + geomapList[0] + ".rsc " + inps.tssar_dir
        print cpCmd+'\n'
        os.system(cpCmd)
    except:
        #print "*********************************"
        print "no geomap file is loaded.\n"
        #print "*********************************\n"


    ##############################################################################
    ##################################  DEM  #####################################

    try:
        demRdrPath = inps.tssar_dir+'/radar*.hgt'
        demRdrList = glob.glob(demRdrPath)
        if len(demRdrList)>0:
            print '\nradar*.hgt'+ '  already exists.'
            sys.exit(1)
    
        demRdrPath=template_dict['pysar.dem.radarCoord']
        demRdrPath=check_variable_name(demRdrPath)
        demRdrList=glob.glob(demRdrPath)
    
        cpCmd="cp " + demRdrList[0] + " " + inps.tssar_dir
        print cpCmd
        os.system(cpCmd)
        cpCmd="cp " + demRdrList[0] + ".rsc " + inps.tssar_dir
        print cpCmd+'\n'
        os.system(cpCmd)
    except:
        #print "*********************************"
        print "no DEM (radar coordinate) file is loaded.\n"
        #print "*********************************"
  
    try:
        demGeoPath = inps.tssar_dir+'/*.dem'
        demGeoList = glob.glob(demGeoPath)
        if len(demGeoList)>0:
            print '\n*.dem'+ '  already exists.'
            sys.exit(1)
    
        demGeoPath=template_dict['pysar.dem.geoCoord']
        demGeoPath=check_variable_name(demGeoPath)
        demGeoList=glob.glob(demGeoPath)
    
        cpCmd="cp " + demGeoList[0] + " " + inps.tssar_dir
        print cpCmd
        os.system(cpCmd)
        cpCmd="cp " + demGeoList[0] + ".rsc " + inps.tssar_dir
        print cpCmd+'\n'
        os.system(cpCmd)
    except:
        #print "*********************************"
        print "no DEM (geo coordinate) file is loaded.\n"
        #print "*********************************\n"


##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
