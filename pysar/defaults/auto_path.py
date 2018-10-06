#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2018, Zhang Yunjun                          #
# Author:  Zhang Yunjun, 2018 Mar                          #
############################################################


import os
import re
import glob
import numpy as np

# Auto setting for file structure of Univ. of Miami, as shown below.
# It required 3 conditions: 1) autoPath = True
#                           2) $SCRATCHDIR is defined in environmental variable
#                           3) input custom template with basename same as project_name
# Change it to False if you are not using it.
autoPath = True


# Default path of data files from different InSAR processors to be loaded into PySAR
isceAutoPath = '''##----------Default file path of ISCE-topsStack products
pysar.load.processor      = isce
pysar.load.metaFile       = ${PROJECT_DIR}/master/IW*.xml
pysar.load.baselineDir    = ${PROJECT_DIR}/baselines

pysar.load.unwFile        = ${PROJECT_DIR}/merged/interferograms/*/filt*.unw
pysar.load.corFile        = ${PROJECT_DIR}/merged/interferograms/*/filt*.cor
pysar.load.connCompFile   = ${PROJECT_DIR}/merged/interferograms/*/filt*.unw.conncomp
pysar.load.intFile        = None

pysar.load.demFile        = ${PROJECT_DIR}/merged/geom_master/hgt.rdr
pysar.load.lookupYFile    = ${PROJECT_DIR}/merged/geom_master/lat.rdr
pysar.load.lookupXFile    = ${PROJECT_DIR}/merged/geom_master/lon.rdr
pysar.load.incAngleFile   = ${PROJECT_DIR}/merged/geom_master/los.rdr
pysar.load.azAngleFile    = ${PROJECT_DIR}/merged/geom_master/los.rdr
pysar.load.shadowMaskFile = ${PROJECT_DIR}/merged/geom_master/shadowMask.rdr
pysar.load.bperpFile      = ${PROJECT_DIR}/merged/baseline_grid/*/bperp.rdr
'''

roipacAutoPath = '''##----------Default file path of ROI_PAC products
pysar.load.processor      = roipac
pysar.load.unwFile        = ${PROJECT_DIR}/PROCESS/DONE/IFG*/filt*.unw
pysar.load.corFile        = ${PROJECT_DIR}/PROCESS/DONE/IFG*/filt*.cor
pysar.load.connCompFile   = ${PROJECT_DIR}/PROCESS/DONE/IFG*/filt*snap_connect.byt
pysar.load.intFile        = ${PROJECT_DIR}/PROCESS/DONE/IFG*/filt*rlks.int

pysar.load.demFile        = ${PROJECT_DIR}/PROCESS/DONE/*${m_date12}*/radar_*rlks.hgt
pysar.load.lookupYFile    = ${PROJECT_DIR}/PROCESS/GEO/geo_${m_date12}/geomap_*rlks.trans
pysar.load.lookupXFile    = ${PROJECT_DIR}/PROCESS/GEO/geo_${m_date12}/geomap_*rlks.trans
pysar.load.incAngleFile   = None
pysar.load.azAngleFile    = None
pysar.load.shadowMaskFile = None
pysar.load.bperpFile      = None
'''

gammaAutoPath = '''##----------Default file path of GAMMA products
pysar.load.processor      = gamma
pysar.load.unwFile        = ${PROJECT_DIR}/PROCESS/DONE/IFG*/diff*rlks.unw
pysar.load.corFile        = ${PROJECT_DIR}/PROCESS/DONE/IFG*/*filt*rlks.cor
pysar.load.connCompFile   = None
pysar.load.intFile        = ${PROJECT_DIR}/PROCESS/DONE/IFG*/diff*rlks.int

pysar.load.demFile        = ${PROJECT_DIR}/PROCESS/SIM/sim_${m_date12}/sim*rlks.rdc.dem
pysar.load.lookupYFile    = ${PROJECT_DIR}/PROCESS/SIM/sim_${m_date12}/sim*rlks.UTM_TO_RDC
pysar.load.lookupXFile    = ${PROJECT_DIR}/PROCESS/SIM/sim_${m_date12}/sim*rlks.UTM_TO_RDC
pysar.load.incAngleFile   = None
pysar.load.azAngleFile    = None
pysar.load.shadowMaskFile = None
pysar.load.bperpFile      = ${PROJECT_DIR}/merged/baselines/*/*.base_perp
'''

prefix = 'pysar.load.'
#config = configparser.ConfigParser()
#config.optionxform = str


##----------------- Functions from pysar.utils.readfile to be independnt module ---------##
def read_str2dict(inString, delimiter='=', print_msg=False):
    '''Read multiple lines of string into dict
    Based on pysar.utils.readfile.read_template()
    '''
    strDict = {}
    lines = inString.split('\n')
    for line in lines:
        c = [i.strip() for i in line.strip().split(delimiter, 1)]
        if len(c) < 2 or line.startswith(('%', '#')):
            next
        else:
            key = c[0]
            value = str.replace(c[1], '\n', '').split("#")[0].strip()
            if value != '':
                strDict[key] = value

    for key, value in strDict.items():
        if value.lower() == 'none':
            strDict[key] = None
    return strDict


##----------------------------------------------------------------------------------------##
def get_auto_path4isce(project_name, template=dict()):
    project_dir = os.path.join(os.getenv('SCRATCHDIR'), project_name)
    auto_dict = read_str2dict(isceAutoPath, print_msg=False)
    for key, value in auto_dict.items():
        if value and template[key] == 'auto':
            value = value.replace
            template[key] = value.replace('${PROJECT_DIR}', project_dir)
    return template


def get_auto_path4roipac(project_name, template=dict()):
    # default file pattern
    auto_dict = read_str2dict(roipacAutoPath, print_msg=False)
    for key, value in auto_dict.items():
        auto_dict[key] = os.path.basename(value)

    project_dir = os.path.join(os.getenv('SCRATCHDIR'), project_name)
    ifgram_dir = os.path.join(project_dir, 'PROCESS', 'DONE', 'IFG*')
    geom_dir = os.path.join(project_dir, 'PROCESS', 'GEO')

    # ifgramStack
    for suffix in ['unwFile', 'corFile', 'connCompFile', 'intFile']:
        key = prefix+suffix
        if template[key] == 'auto':
            template[key] = os.path.join(ifgram_dir, auto_dict[key])

    # m_date12
    m_date12 = None     #dates of master interferogram in YYMMDD-YYMMDD format
    masterIfgramTxtFile = os.path.join(project_dir, 'PROCESS', 'master_ifgram.txt')
    if os.path.isfile(masterIfgramTxtFile):
        m_date12 = str(np.loadtxt(masterIfgramTxtFile, dtype=bytes).astype(str))
    else:
        try:
            lookup_file = glob.glob(os.path.join(geom_dir, 'geo_*/geomap*.trans'))[0]
            m_date12 = re.findall('\d{6}-\d{6}', lookup_file)[0]
        except:
            print(("No master interferogram found!"
                   "Check the {}/geo_* folder").format(geom_dir))
            m_date12 = None

    # mli_looks in case both radar_2rlks.hgt and radar_8rlks.hgt exist.
    if m_date12:
        lookup_file = glob.glob(os.path.join(geom_dir, 'geo_{}/geomap*.trans'.format(m_date12)))[0]
        mli_looks = re.findall('_\d{1}rlks', lookup_file)[0]
        auto_dict[prefix+'demFile'] = 'radar{}.hgt'.format(mli_looks)

    for suffix in ['demFile', 'incAngleFile', 'azAngleFile', 'shadowMaskFile']:
        key = prefix+suffix
        if template[key] == 'auto':
            if m_date12 and auto_dict[key] != 'None':
                template[key] = os.path.join(project_dir,
                                             'PROCESS/DONE/*{}*'.format(m_date12),
                                             auto_dict[key])
            else:
                template[key] = None

    for suffix in ['lookupYFile', 'lookupXFile']:
        key = prefix+suffix
        if template[key] == 'auto':
            if m_date12:
                template[key] = os.path.join(geom_dir,
                                             '*{}*'.format(m_date12),
                                             auto_dict[key])
            else:
                template[key] = None

    return template


def get_auto_path4gamma(project_name, template=dict()):
    # default file pattern
    auto_dict = read_str2dict(gammaAutoPath, print_msg=False)
    for key, value in auto_dict.items():
        auto_dict[key] = os.path.basename(value)

    project_dir = os.path.join(os.getenv('SCRATCHDIR'), project_name)
    ifgram_dir = os.path.join(project_dir, 'PROCESS', 'DONE', 'IFG*')
    geom_dir = os.path.join(project_dir, 'PROCESS', 'SIM')

    # ifgramStack
    for suffix in ['unwFile', 'corFile', 'connCompFile', 'intFile']:
        key = prefix+suffix
        if template[key] == 'auto':
            if auto_dict[key] != 'None':
                template[key] = os.path.join(ifgram_dir, auto_dict[key])
            else:
                template[key] = None

    # geometry
    m_date12 = None
    masterIfgramTxtFile = os.path.join(
        project_dir, 'PROCESS', 'master_ifgram.txt')
    if os.path.isfile(masterIfgramTxtFile):
        m_date12 = str(np.loadtxt(masterIfgramTxtFile, dtype=bytes).astype(str))
    else:
        try:
            m_date12 = os.walk(geom_dir).next()[1][0].split('sim_')[1]
        except:
            print("No master interferogram found! Check the {} folder".format(
                os.path.join(geom_dir, 'sim_')))
            m_date12 = None

    for suffix in ['demFile', 'lookupYFile', 'lookupXFile',
                   'incAngleFile', 'azAngleFile', 'shadowMaskFile']:
        key = prefix+suffix
        if template[key] == 'auto':
            if m_date12 and auto_dict[key] != 'None':
                template[key] = os.path.join(geom_dir, '*'+m_date12+'*', auto_dict[key])
            else:
                template[key] = None

    return template
