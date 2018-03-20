#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2018, Zhang Yunjun                          #
# Author:  Zhang Yunjun, 2018 Mar                          #
############################################################


import os, sys
import numpy as np
from pysar.utils import utils as ut


def default_path(inps, template={}):
    '''Auto File Path Setting for Geodesy Lab - University of Miami'''
    print('Use auto path setting in University of Miami.'+\
          '(To turn it off, change auto_path_miami value to False in pysar/__init__.py)')
    # PYSAR working directory
    if not inps.timeseries_dir:
        inps.timeseries_dir = os.getenv('SCRATCHDIR')+'/'+inps.project_name+'/PYSAR'

    ##### .unw/.cor/.int files
    process_dir = os.getenv('SCRATCHDIR')+'/'+inps.project_name+'/PROCESS'
    print("PROCESS directory: "+process_dir)
    if inps.insarProcessor == 'roipac':
        if not inps.unw or inps.unw == 'auto':   inps.unw = process_dir+'/DONE/IFGRAM*/filt_*.unw'
        if not inps.cor or inps.cor == 'auto':   inps.cor = process_dir+'/DONE/IFGRAM*/filt_*rlks.cor'
        #if not inps.int or inps.int == 'auto':   inps.int = process_dir+'/DONE/IFGRAM*/filt_*rlks.int'
    elif inps.insarProcessor == 'gamma':
        if not inps.unw or inps.unw == 'auto':   inps.unw = process_dir+'/DONE/IFGRAM*/diff_*rlks.unw'
        if not inps.cor or inps.cor == 'auto':   inps.cor = process_dir+'/DONE/IFGRAM*/filt_*rlks.cor'
        #if not inps.int or inps.int == 'auto':   inps.int = process_dir+'/DONE/IFGRAM*/diff_*rlks.int'
    elif inps.insarProcessor == 'isce':
        process_dir = os.getenv('SCRATCHDIR')+'/'+inps.project_name
        if not inps.unw or inps.unw == 'auto':   inps.unw = process_dir+'/merged/interferograms/*/filt*.unw'
        if not inps.cor or inps.cor == 'auto':   inps.cor = process_dir+'/merged/interferograms/*/filt*.cor'
        if not inps.lut or inps.lut == 'auto':   inps.lut = process_dir+'/merged/geom_master/l*.rdr'
        if not inps.dem_radar or inps.dem_radar == 'auto':   inps.dem_radar = process_dir+'/merged/geom_master/hgt.rdr'
        if not inps.dem_geo or inps.dem_geo == 'auto':   inps.dem_geo = ""
        #if not inps.int or inps.int == 'auto':   inps.int = process_dir+'/DONE/IFGRAM*/diff_*rlks.int'

    ##### master interferogram for lookup table and DEM in radar coord
    if all(fname and fname != 'auto' for fname in [inps.lut, inps.dem_radar, inps.dem_geo]):
        return inps

    try:     m_date12 = str(np.loadtxt(process_dir+'/master_ifgram.txt', dtype=bytes).astype(str))
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
    if inps.dem_geo is None or inps.dem_geo == 'auto':
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








