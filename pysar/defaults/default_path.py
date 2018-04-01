#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2018, Zhang Yunjun                          #
# Author:  Zhang Yunjun, 2018 Mar                          #
############################################################


import os
import numpy as np
from pysar.utils import readfile


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
pysar.load.bperpFile      = $PROJECT_DIR/baselines/bperp_*.rdr
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
'''

gammaAutoPath = '''##----------Default file path of GAMMA products
pysar.load.processor      = gamma
pysar.load.unwFile        = $PROJECT_DIR/PROCESS/DONE/IFG*/diff*rlks.unw
pysar.load.corFile        = $PROJECT_DIR/PROCESS/DONE/IFG*/*filt*rlks.cor
pysar.load.connCompFile   = None
pysar.load.intFile        = $PROJECT_DIR/PROCESS/DONE/IFG*/diff*rlks.int

pysar.load.demFile        = $PROJECT_DIR/PROCESS/SIM/sim_${m_date12}/sim*rlks.rdc.hgt
pysar.load.lookupYFile    = $PROJECT_DIR/PROCESS/SIM/sim_${m_date12}/sim*rlks.UTM_TO_RDC
pysar.load.lookupXFile    = $PROJECT_DIR/PROCESS/SIM/sim_${m_date12}/sim*rlks.UTM_TO_RDC
pysar.load.incAngleFile   = None
pysar.load.headAngleFile  = None
pysar.load.shadowMaskFile = None
'''

prefix = 'pysar.load.'
#config = configparser.ConfigParser()
#config.optionxform = str

##----------------------------------------------------------------------------------------##
def default_path4sentinel_stack(projectName, template=dict()):
    ## default file pattern
    defDict = read_file.read_template(isceAutoPath, printMsg=False)
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


def default_path4roipac(projectName, template=dict()):
    ## default file pattern
    defDict = read_file.read_template(roipacAutoPath, printMsg=False)
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


def default_path4gamma(projectName, template=dict()):
    ## default file pattern
    defDict = read_file.read_template(gammaAutoPath, printMsg=False)
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

