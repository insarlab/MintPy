#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2018, Zhang Yunjun                          #
# Author:  Zhang Yunjun, 2018 Mar                          #
############################################################


import os
import numpy as np

## Auto setting for file structure of Univ. of Miami, as shown below. 
# It required 3 conditions: 1) autoPath = True
#                           2) $SCRATCHDIR is defined in environmental variable
#                           3) input custom template with basename same as projectName
# Change it to False if you are not using it.
autoPath = True


###### Default path of data files from different InSAR processors to be loaded into PySAR
isceAutoPath = '''##----------Default file path of ISCE/SentinelStack products
pysar.load.processor      = isce
pysar.load.unwFile        = $PROJECT_DIR/merged/interferograms/*/filt*.unw
pysar.load.corFile        = $PROJECT_DIR/merged/interferograms/*/filt*.cor
pysar.load.connCompFile   = $PROJECT_DIR/merged/interferograms/*/filt*.unw.conn*
pysar.load.intFile        = $PROJECT_DIR/merged/interferograms/*/filt*.int

pysar.load.demFile        = $PROJECT_DIR/merged/geom_master/hgt.rdr
pysar.load.lookupYFile    = $PROJECT_DIR/merged/geom_master/lat.rdr
pysar.load.lookupXFile    = $PROJECT_DIR/merged/geom_master/lon.rdr
pysar.load.incAngleFile   = $PROJECT_DIR/merged/geom_master/los.rdr
pysar.load.headAngleFile  = $PROJECT_DIR/merged/geom_master/los.rdr
pysar.load.shadowMaskFile = $PROJECT_DIR/merged/geom_master/shadowMask.rdr
pysar.load.bperpFile      = $PROJECT_DIR/merged/baselines/*/bperp
'''

roipacAutoPath = '''##----------Default file path of ROI_PAC products
pysar.load.processor      = roipac
pysar.load.unwFile        = $PROJECT_DIR/PROCESS/DONE/IFG*/filt*rlks_c*.unw
pysar.load.corFile        = $PROJECT_DIR/PROCESS/DONE/IFG*/filt*rlks.cor
pysar.load.connCompFile   = $PROJECT_DIR/PROCESS/DONE/IFG*/filt*rlks_snap*.byt
pysar.load.intFile        = $PROJECT_DIR/PROCESS/DONE/IFG*/filt*rlks.int

pysar.load.demFile        = $PROJECT_DIR/PROCESS/DONE/*${m_date12}*/radar*.hgt
pysar.load.lookupYFile    = $PROJECT_DIR/PROCESS/GEO/geo_${m_date12}/geomap*.trans
pysar.load.lookupXFile    = $PROJECT_DIR/PROCESS/GEO/geo_${m_date12}/geomap*.trans
pysar.load.incAngleFile   = None
pysar.load.headAngleFile  = None
pysar.load.shadowMaskFile = None
pysar.load.bperpFile      = None
'''

gammaAutoPath = '''##----------Default file path of GAMMA products
pysar.load.processor      = gamma
pysar.load.unwFile        = $PROJECT_DIR/PROCESS/DONE/IFG*/diff*rlks.unw
pysar.load.corFile        = $PROJECT_DIR/PROCESS/DONE/IFG*/*filt*rlks.cor
pysar.load.connCompFile   = None
pysar.load.intFile        = $PROJECT_DIR/PROCESS/DONE/IFG*/diff*rlks.int

pysar.load.demFile        = $PROJECT_DIR/PROCESS/SIM/sim_${m_date12}/sim*rlks.rdc.dem
pysar.load.lookupYFile    = $PROJECT_DIR/PROCESS/SIM/sim_${m_date12}/sim*rlks.UTM_TO_RDC
pysar.load.lookupXFile    = $PROJECT_DIR/PROCESS/SIM/sim_${m_date12}/sim*rlks.UTM_TO_RDC
pysar.load.incAngleFile   = None
pysar.load.headAngleFile  = None
pysar.load.shadowMaskFile = None
pysar.load.bperpFile      = $PROJECT_DIR/merged/baselines/*/*.base_perp
'''

prefix = 'pysar.load.'
#config = configparser.ConfigParser()
#config.optionxform = str


##----------------- Functions from pysar.utils.readfile to be independnt module ---------##
def check_variable_name(path, print_msg=True):
    s = path.split("/")[0]
    if len(s)>0 and s[0]=="$":
        try:
            p0 = os.getenv(s[1:])
            path = path.replace(path.split("/")[0], p0)
        except:
            if print_msg:
                print('WARNING: Un-recognized environmental variable: '+s)
    return path


def read_str2dict(inString, delimiter='=', print_msg=False):
    '''Read multiple lines of string into dict
    Based on pysar.utils.readfile.read_template()
    '''
    strDict = {}
    lines = inString.split('\n')
    for line in lines:
        c = [i.strip() for i in line.strip().split(delimiter, 1)]
        if len(c) < 2 or line.startswith(('%','#')):
            next
        else:
            key  = c[0]
            value = str.replace(c[1],'\n','').split("#")[0].strip()
            value = check_variable_name(value, print_msg=print_msg)
            if value != '':
                strDict[key] = value
    return strDict
                

##----------------------------------------------------------------------------------------##
def get_auto_path4sentinel_stack(projectName, template=dict()):
    ## default file pattern
    defDict = read_str2dict(isceAutoPath, print_msg=False)
    for key, value in defDict.item():
        defDict[key] = os.path.basename(value)

    projectDir = os.path.join(os.getenv('SCRATCHDIR'),projectName)
    ifgramDir = os.path.join(projectDir,'merged','interferograms','*')
    geomDir = os.path.join(projectDir,'merged','geom_master')

    ## ifgramStack
    for suffix in ['unwFile','corFile','connCompFile','intFile']:
        key = prefix+suffix
        if template[key] == 'auto':
            template[key] = os.path.join(ifgramDir, defDict[key])

    ## geometry
    for suffix in ['demFile','lookupYFile','lookupXFile','incAngleFile','headAngleFile','shadowMaskFile']:
        key = prefix+suffix
        if template[key] == 'auto':
            template[key] = os.path.join(geomDir, defDict[key])

    return template


def get_auto_path4roipac(projectName, template=dict()):
    ## default file pattern
    defDict = read_str2dict(roipacAutoPath, print_msg=False)
    for key, value in defDict.item():
        defDict[key] = os.path.basename(value)

    projectDir = os.path.join(os.getenv('SCRATCHDIR'),projectName)
    ifgramDir = os.path.join(projectDir,'PROCESS','DONE','IFG*')
    geomDir = os.path.join(projectDir,'PROCESS','GEO')

    ## ifgramStack
    for suffix in ['unwFile','corFile','connCompFile','intFile']:
        key = prefix+suffix
        if template[key] == 'auto':
            template[key] = os.path.join(ifgramDir, defDict[key])

    ## geometry
    m_date12 = None
    masterIfgramTxtFile = os.path.join(projectDir,'PROCESS','master_ifgram.txt')
    if os.path.isfile(masterIfgramTxtFile):
        m_date12 = str(np.loadtxt(masterIfgramTxtFile, dtype=bytes).astype(str))
    else:
        try:
            m_date12 = os.walk(geomDir).next()[1][0].split('geo_')[1]
        except:
            print("No master interferogram found! Check the {} folder".format(os.path.join(geomDir,'geo_')))
            m_date12 = None

    for suffix in ['demFile','incAngleFile','headAngleFile','shadowMaskFile']:
        key = prefix+suffix
        if template[key] == 'auto':
            if m_date12 and defDict[key] != 'None':
                template[key] = os.path.join(projectDir,'PROCESS','DONE','*'+m_date12+'*',defDict[key])
            else:
                template[key] = None

    for suffix in ['lookupYFile','lookupXFile']:
        key = prefix+suffix
        if template[key] == 'auto':
            if m_date12:
                template[key] = os.path.join(geomDir,'*'+m_date12+'*',defDict[key])
            else:
                template[key] = None

    return template


def get_auto_path4gamma(projectName, template=dict()):
    ## default file pattern
    defDict = read_str2dict(gammaAutoPath, print_msg=False)
    for key, value in defDict.item():
        defDict[key] = os.path.basename(value)

    projectDir = os.path.join(os.getenv('SCRATCHDIR'),projectName)
    ifgramDir = os.path.join(projectDir,'PROCESS','DONE','IFG*')
    geomDir = os.path.join(projectDir,'PROCESS','SIM')

    ## ifgramStack
    for suffix in ['unwFile','corFile','connCompFile','intFile']:
        key = prefix+suffix
        if template[key] == 'auto':
            if defDict[key] != 'None':
                template[key] = os.path.join(ifgramDir, defDict[key])
            else:
                template[key] = None

    ## geometry
    m_date12 = None
    masterIfgramTxtFile = os.path.join(projectDir,'PROCESS','master_ifgram.txt')
    if os.path.isfile(masterIfgramTxtFile):
        m_date12 = str(np.loadtxt(masterIfgramTxtFile, dtype=bytes).astype(str))
    else:
        try:
            m_date12 = os.walk(geomDir).next()[1][0].split('sim_')[1]
        except:
            print("No master interferogram found! Check the {} folder".format(os.path.join(geomDir,'sim_')))
            m_date12 = None

    for suffix in ['demFile','lookupYFile','lookupXFile','incAngleFile','headAngleFile','shadowMaskFile']:
        key = prefix+suffix
        if template[key] == 'auto':
            if m_date12 and defDict[key] != 'None':
                template[key] = os.path.join(geomDir,'*'+m_date12+'*',defDict[key])
            else:
                template[key] = None

    return template

