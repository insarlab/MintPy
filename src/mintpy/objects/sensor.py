"""Utilities for SAR sensor parameters."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2016                               #
############################################################
# Recommend import:
#   from mintpy.objects import sensor


import os

SENSOR_NAME_VARIATION = {
    'alos'  : ['alos', 'alos1', 'palsar', 'palsar1'],
    'alos2' : ['alos2', 'palsar2'],
    'alos4' : ['alos4', 'palsar3'],
    'bio'   : ['bio', 'biomass'],
    'csk'   : ['csk', 'csk1', 'csk2', 'csk3', 'csk4', 'cos', 'cosmo', 'cosmoskymed'],
    'env'   : ['env', 'envisat', 'asar'],
    'ers'   : ['ers', 'ers1', 'ers2', 'ers12'],
    'gf3'   : ['gfen3', 'gaofen3', 'g3', 'gaofen'],
    'hj1c'  : ['hj1', 'huanjing1c'],
    'jers'  : ['jers', 'jers1'],
    'ksat5' : ['ksat5', 'kompsat5', 'kompsat', 'kmps5'],
    'lt1'   : ['lt1', 'lt', 'lutan', 'lutan1'],
    'ni'    : ['ni', 'nisar'],
    'rs1'   : ['rs1', 'rsat', 'rsat1', 'radarsat', 'radarsat1'],
    'rs2'   : ['rs2', 'rsat2', 'radarsat2'],
    'rcm'   : ['rcm', 'rsatc', 'radarsat-constellation', 'radarsat-constellation-mission'],
    'sen'   : ['sen', 's1', 's1a', 's1b', 'sent1', 'sentinel1', 'sentinel1a', 'sentinel1b'],
    'tsx'   : ['tsx', 'terra', 'terrasar', 'terrasarx', 'tdx', 'tandemx'],
    'uav'   : ['uav', 'uavsar'],
}

# duplicated in mintpy.cli.prep_gamma
SENSOR_NAMES = list(SENSOR_NAME_VARIATION.keys())


###########################  Util functions  ###################################

def standardize_sensor_name(sensor_name):
    """"""
    # decode if encoded
    try:
        sensor_name = sensor_name.decode()
    except (UnicodeDecodeError, AttributeError):
        pass

    # remove -_ and user lower case before standardize sensor names
    sensor_name = sensor_name.replace('-', '').replace('_', '').lower()

    if sensor_name in SENSOR_NAMES:
        # if input name is already standardized, do nothing
        pass
    else:
        # otherwise, check all the possible variations
        for key, values in SENSOR_NAME_VARIATION.items():
            if sensor_name in values:
                sensor_name = key
            else:
                pass
    return sensor_name


def project_name2sensor_name(proj_names):
    """Get sensor name from project_name or path
    Parameters: proj_names : str or list of str, name or path of template file containing project name
    Returns:    sensor     : str, SAR sensor name
                proj_name  : str, project name
    Examples:   ('Sen',
                 'AlcedoSenDT128') = project_name2sensor_name('AlcedoSenDT128')
                ('Env',
                 'GalapagosEnvA2T061') = project_name2sensor_name('~/insarlab/GalapagosEnvA2T061/'+
                                                                  'mintpy/GalapagosEnvA2T061.template')
    """
    sensor = None
    proj_name = None

    if isinstance(proj_names, str):
        proj_names = [proj_names]
    proj_names = [p for p in proj_names if p is not None]

    # get proj_name if input is the path of template file
    for proj_path in proj_names:
        if any(s in proj_path.lower() for s in SENSOR_NAMES) and proj_name is None:
            # exclude known words overlapping with sensor names
            # ers: users
            proj_path_segs = [p for p in proj_path.split(os.sep)
                              if p.lower() not in ['users']]
            for proj_path_seg in proj_path_segs:
                if any(s.capitalize() in proj_path_seg for s in SENSOR_NAMES):
                    proj_name = os.path.splitext(proj_path_seg)[0]

    # proj_name --> sensor
    if proj_name:
        # potetial sensor names
        # priority: captialized_case > all_upper_case > all_lower_case
        sensors = [s.capitalize() for s in SENSOR_NAMES if s.capitalize() in proj_name]
        if len(sensors) == 0:
            sensors = [s.capitalize() for s in SENSOR_NAMES if s.upper() in proj_name]
        if len(sensors) == 0:
            sensors = [s.capitalize() for s in SENSOR_NAMES if s.lower() in proj_name]

        if len(sensors) > 0:
            # if more than one, i.e. ['Alos','Alos2'], use the last one
            sensor = sorted(sensors)[-1]
        else:
            msg = f"No sensor name found in project_name: {proj_name}\n"
            msg += f"Available sensor names: {SENSOR_NAMES}"
            raise ValueError(msg)

    return sensor, proj_name


def get_unavco_mission_name(meta_dict):
    """Get mission name in UNAVCO InSAR Archive format from attribute mission/PLATFORM
    Parameters: meta_dict : dict, attributes
    Returns:    mission   : str, mission name in standard UNAVCO format.
    """
    mission_name = None

    if 'mission' in meta_dict.keys():
        value = meta_dict['mission'].lower()
    elif 'PLATFORM' in meta_dict.keys():
        value = meta_dict['PLATFORM'].lower()
    else:
        print('No PLATFORM nor mission attribute found, can not identify mission name.')
        print('return None')
        return mission_name

    # Convert to UNAVCO Mission name
    ## ERS, ENV, S1, RS1, RS2, CSK, TSX, JERS, ALOS, ALOS2
    if value.startswith(('alos', 'palsar')):
        if value.endswith('2'):
            mission_name = 'ALOS2'
        else:
            mission_name = 'ALOS'

    elif value.startswith(('csk', 'cos')):
        mission_name = 'CSK'

    elif value.startswith(('env', 'asar')):
        mission_name = 'ENV'

    elif value.startswith('ers'):
        mission_name = 'ERS'

    elif value.startswith('jers'):
        mission_name = 'JERS'

    elif value.startswith(('rs', 'rsat', 'radarsat')):
        mission_name = 'RS'
        if value.endswith('1'):
            mission_name += '1'
        else:
            mission_name += '2'

    elif value.startswith(('s1', 'sen')):
        mission_name = 'S1'

    elif value.startswith(('tsx', 'tdx', 'terra', 'tandem')):
        mission_name = 'TSX'

    elif value.startswith('uav'):
        mission_name = 'UAV'

    else:
        print('Un-recognized PLATFORM attribute:', value)
        print('return None')
    return mission_name



###########################  Hardwired parameters  #############################
## References:
# Giudici, D., A. Monti Guarnieri, and J. P. Cuesta Gonzalez (2017), Pre-Flight
#   SAOCOM-1A SAR Performance Assessment by Outdoor Campaign, Remote Sensing, 9(7),
#   729, doi:10.3390/rs9070729.
# Guarnieri, A. M. (2013), Introduction to RADAR, Politecnico di Milano Dipartimento
#   di Elettronica e Informazione, Milano.
# Kim, Y., and R. L. Jordan (2006), Spaceborne SAR Antennas for Earth Science,
#   in Spaceborne Antennas for Planetary Exploration, edited by W. A. Imbriale,
#   pp. 305-340, doi:10.1002/0470052783.ch6.
# Jung, H.-S., W.-J. Lee, and L. Zhang (2014), Theoretical accuracy of along-track
#   displacement measurements from multiple-aperture interferometry (MAI), Sensors (Basel),
#   14(9), 17703-17724, doi:10.3390/s140917703.
# Zebker, H. A., C. L. Werner, P. A. Rosen, and S. Hensley (1994), Accuracy of topographic maps
#   derived from ERS-1 interferometric radar, Geoscience and Remote Sensing, IEEE Transactions on,
#   32(4), 823-836, doi:10.1109/36.298010.


## Signal-to-noise-ratio (SNR) and Noise Equivalent Sigma Zero (NESZ):
# ERS        : SNR = 11.7 dB from Table 2 in Zebker et al. (1994)
# Envisat    : SNR = 19.5 dB from Table 3.3 in Guarnieri (2013)
# Sentinel-1 : NESZ = -22 dB from Table 1 in https://sentinels.copernicus.eu/web/sentinel/
#     user-guides/sentinel-1-sar/acquisition-modes/interferometric-wide-swath
# JERS       : SNR = 14 dB from https://directory.eoportal.org/web/eoportal/satellite-missions/j/jers-1


##--------------------  Ka-band  -------------------##
# SWOT
# launch date: 2022-12-16
# end    date: operational
# KaRIn (Ka-band Radar Interferometer) On-Board Processor (OBP) ATBD
SWOT = {
    # orbit
    'altitude'                   : 906e3,     # m
    'orbit_inclination'          : 78,        # deg
    # sar / antenna
    'carrier_frequency'          : 35.75e9,   # Hz
    'pulse_repetition_frequency' : 4420,      # Hz
    'chirp_bandwidth'            : 200e6,     # Hz
    'sampling_frequency'         : 300e6,     # Hz
}



##--------------------  X-band  --------------------##
# TerraSAR-X stripmap mode in single polarization
# launch date: 2007-06-15 (TSX), 2010-06-21 (TDX)
# end    date: operational
# from Table 1 in Jung et al. (2014)
# https://www.eoportal.org/satellite-missions/terrasar-x
TSX = {
    # orbit
    'altitude'                   : 514.8e3,   # m, mean value, 505-533 km
    'orbit_inclination'          : 97.44,     # deg
    'repeat_cycle'               : 11,        # day
    # sar / antenna
    'carrier_frequency'          : 9.65e9,    # Hz
    'antenna_length'             : 4.8,       # m
    'antenna_width'              : 0.8,       # m
    'doppler_bandwidth'          : 2770,      # Hz
    'pulse_repetition_frequency' : 3800,      # Hz
    'chirp_bandwidth'            : 100e6,     # Hz
    'sampling_frequency'         : 109.89e6,  # Hz
    'azimuth_pixel_size'         : 2.0,       # m
    'ground_range_pixel_size'    : 2.2,       # m
}

# COSMO-SkyMed stripmap HIMAGE mode
# launch date: 2007-06-08 (CSK1), 2007-12-09 (CSK2), 2008-10-25 (CSK3), 2010-11-06 (CSK4)
# end    date: operational
# from Table 1 in Jung et al. (2014)
# https://www.eoportal.org/satellite-missions/cosmo-skymed
CSK = {
    # orbit
    'altitude'                   : 619.6e3,   # m, mean value
    'orbit_inclination'          : 97.86,     # deg
    'repeat_cycle'               : 16,        # day, single satellite
    # sar / antenna
    'carrier_frequency'          : 9.6e9,     # Hz
    'antenna_length'             : 5.7,       # m
    'antenna_width'              : 1.4,       # m
    'doppler_bandwidth'          : 2670,      # Hz
    'pulse_repetition_frequency' : 3000,      # Hz
    'chirp_bandwidth'            : 117e6,     # Hz
    'sampling_frequency'         : 146.25e6,  # Hz
    'azimuth_pixel_size'         : 2.4,       # m
    'ground_range_pixel_size'    : 1.6,       # m
}

# Kompsat-5 (Korea Multi-Purpose Satellite-5) stripmap mode
# launch date: 2013-08-22
# end    date: operational
# from Table 1 in Jung et al. (2014)
# https://www.eoportal.org/satellite-missions/kompsat-5
KSAT5 = {
    # orbit
    'altitude'                   : 550e3,     # m, mean value
    'orbit_inclination'          : 98.1,      # deg
    'repeat_cycle'               : 28,        # day
    # sar / antenna
    'carrier_frequency'          : 9.66e9,    # Hz
    'antenna_length'             : 4.48,      # m
    'doppler_bandwidth'          : 3110,      # Hz
    'pulse_repetition_frequency' : 3530,      # Hz
    'chirp_bandwidth'            : 73.24e6,   # Hz
    'sampling_frequency'         : 88.125e6,  # Hz
    'azimuth_pixel_size'         : 2.1,       # m
    'ground_range_pixel_size'    : 2.7,       # m
}

# ICEYE SAR constellation
# from Table 2.1 in ICEYE SAR Product Guide at:
# https://earth.esa.int/eogateway/documents/20142/37627/ICEYE-SAR-Product-Guide-V4.pdf
# https://www.eoportal.org/satellite-missions/iceye-constellation
ICEYE = {
    # orbit
    'altitude'                   : 570e3,           # m, mean value, 560-580 km
    'orbit_inclination'          : 97.7,            # deg
    'repeat_cycle'               : 17,              # day, single satellite
    # sar / antenna
    'carrier_frequency'          : 9.65e9,          # Hz
    'antenna_length'             : 3.2,             # m
    'antenna_width'              : 0.4,             # m
    'pulse_repetition_frequency' : [2e3, 10e3],     # Hz
    'chirp_bandwidth'            : [37.6e6, 299e6], # Hz
}


##--------------------  C-band  --------------------##

# ERS-1/2
# launch date: 1991-07-17 (ERS-1), 1995-04-21 (ERS-2)
# end    date: 2000-03-10 (ERS-1), 2011-09-05 (ERS-2)
# from Table 2 in Jung et al. (2014)
# from Imaging Radar class by Howard Zebker, 2021.
# https://www.esa.int/esapub/bulletin/bullet83/duc83.htm
# https://www.eoportal.org/satellite-missions/ers-1
ERS = {
    # orbit
    'altitude'                   : 785e3,       # m, mean value, 782-785 km
    'orbit_inclination'          : 98.52,       # deg
    'repeat_cycle'               : 35,          # day
    # sar / antenna
    'carrier_frequency'          : 5.30e9,      # Hz
    'antenna_length'             : 10.0,        # m
    'antenna_width'              : 1.0,         # m
    'doppler_bandwidth'          : 1500,        # Hz
    'pulse_repetition_frequency' : 1679.9,      # Hz
    'pulse_length'               : 37.12e-6,    # s
    'chirp_bandwidth'            : 15.55e6,     # Hz
    'chirp_slope'                : 4.189166e11, # Hz
    'sampling_frequency'         : 18.96e6,     # Hz
    'azimuth_pixel_size'         : 4.2,         # m
    'ground_range_pixel_size'    : 20.2,        # m
}

# Envisat
# lauch date: 2002-03-01
# end   date: 2012-04-08
# from Table 2 in Jung et al. (2014)
# https://earth.esa.int/eogateway/missions/envisat/description
# https://www.eoportal.org/satellite-missions/envisat#asar-advanced-sar
ENV = {
    # orbit
    'altitude'                   : 799.8e3,   # m, mean value, 780-820 km
    'orbit_inclination'          : 98.55,     # deg
    'repeat_cycle'               : 35,        # day
    # sar / antenna
    'carrier_frequency'          : 5.331e9,   # Hz
    'antenna_length'             : 10.0,      # m
    'antenna_width'              : 1.3,       # m
    'doppler_bandwidth'          : 1500,      # Hz
    'pulse_repetition_frequency' : 1650,      # Hz
    'chirp_bandwidth'            : 16.00e6,   # Hz
    'sampling_frequency'         : 18.00e6,   # Hz
    'azimuth_pixel_size'         : 4.3,       # m
    'ground_range_pixel_size'    : 21.3,      # m
}

# Radarsat-1
# launch date: 1995-11-04
# end    date: 2013-03-29
# https://www.asc-csa.gc.ca/eng/satellites/radarsat/technical-features/radarsat-comparison.asp
RSAT1 = {
    # orbit
    'altitude'                   : 807e3,     # m, mean value, 793-821 km
    'orbit_inclination'          : 98.6,      # deg
    'repeat_cycle'               : 14,        # day, 14.29 orbits per day (14 7/24)
    # sar / antenna
    'carrier_frequency'          : 5.3e9,     # Hz
    'antenna_length'             : 15,        # m
    'antenna_width'              : 1.5,       # m
}

# Radarsat-2 stripmap ultra-fine mode
# launch date: 2007-12-14
# end    date: operational
# from Table 2 in Jung et al. (2014)
RSAT2 = {
    # orbit
    'altitude'                   : 798e3,     # m
    'orbit_inclination'          : 98.6,      # deg
    'repeat_cycle'               : 14,        # day, 14.29 orbits per day (14 7/24)
    # sar / antenna
    'carrier_frequency'          : 5.405e9,   # Hz
    'antenna_length'             : 6.55,      # m
    'doppler_bandwidth'          : 2308,      # Hz
    'pulse_repetition_frequency' : 3637,      # Hz
    'chirp_bandwidth'            : 78.16e6,   # Hz
    'sampling_frequency'         : 112.68e6,  # Hz
    'azimuth_pixel_size'         : 2.2,       # m
    'ground_range_pixel_size'    : 2.1,       # m
}

# Radarsat Constellation Mission
# launch date: 2019-06-12
# end    date: operational
# https://www.asc-csa.gc.ca/eng/satellites/radarsat/technical-features/radarsat-comparison.asp
# Cote et al. (2021) at https://ieeexplore.ieee.org/document/9472534
RCM = {
    # orbit
    'altitude'                   : 600.5e3,   # m, mean value, 586-615 km
    'orbit_inclination'          : 97.74,     # deg
    'repeat_cycle'               : 14,        # day, 14.92 orbits per day (14 11/12)
    # sar / antenna
    'carrier_frequency'          : 5.405e9,   # Hz
    'antenna_length'             : 6.75,      # m
    'antenna_width'              : 1.38,      # m
    'chirp_bandwidth'            : 100e6,     # Hz
    'noise_equivalent_sigma_zero': -19,       # dB, [-17, -25] for all modes
}

# GaoFen-3
# launch date: 2016-08-10 (GF3-01), 2021-11-23 (GF3-02), 2022-04-07 (GF3-03)
# end    date: operational
# Table 2 & 6 in https://directory.eoportal.org/web/eoportal/satellite-missions/g/gaofen-3
# https://www.eoportal.org/satellite-missions/gaofen-3
# Li et al. (2018, RS) at https://doi.org/10.3390/rs10121929
# Table I in Yang et al. (2023, IEEE-TGRS) at https://doi.org/10.1109/TGRS.2023.3238707
GF3 = {
    # orbit
    'altitude'                   : 755e3,     # m
    'orbit_inclination'          : 98.41,     # deg
    'repeat_cycle'               : 29,        # day
    # sar / antenna
    'carrier_frequency'          : 5.4e9,     # Hz
    'antenna_length'             : 15,        # m
    'antenna_width'              : 1.232,     # m
    'pulse_repetition_frequency' : 1412.18,   # Hz
    'chirp_bandwidth'            : 60.00e6,   # Hz
    'sampling_frequency'         : 533.33e6,  # Hz, IF sampling
    'azimuth_pixel_size'         : 4.77,      # m, FSII mode
    'range_pixel_size'           : 2.25,      # m, FSII mode
    'noise_equivalent_sigma_zero': -19.5,     # dB, -19.5/-21.3 dB for P,UF,F,WF,Q,WV/S,NS,WS,G,WQ,E modes
}

# Sentinel-1 Interferometric Wide (IW / TOPS) swath mode
# launch date: 2014-04-03  (S1A), 2016-04-25 (S1B), 2024-12-05  (S1C)
# end    date: operational (S1A), 2021-12-23 (S1B), operational (S1C)
# Typical value:
# azfact = azResolution / azPixelSize = 1.46
# rgfact = rgResolution / rgPixelSize = 1.33
# reference:
#   1. Table 2 & Fig. 5d in Jung et al. (2014)
#   2. Table 3-1 & 7-5 in https://sentinel.esa.int/documents/247904/1877131/Sentinel-1-Product-Definition
#   3. Fig. 4 in Miranda et al. (2016) at https://doi.org/10.1016/j.procs.2016.09.247
#   4. https://www.eoportal.org/satellite-missions/copernicus-sentinel-1
SEN = {
    # orbit
    'altitude'                   : 705e3,     # m, mean value
    'orbit_inclination'          : 98.18,     # deg
    'repeat_cycle'               : 12,        # day, single satellite
    # sar / antenna
    'carrier_frequency'          : 5.405e9,   # Hz, 5.55 cm wavelength
    'antenna_length'             : 12.3,      # m
    'antenna_width'              : 0.821,      # m
    'doppler_bandwidth'          : 380,       # Hz
    'pulse_repetition_frequency' : 1717.13,   # Hz, based on real data; 1000-3000 (programmable)
    'chirp_bandwidth'            : 56.50e6,   # Hz
    'sampling_frequency'         : 64.35e6,   # Hz
    'azimuth_pixel_size'         : 14.1,      # m, this is the ground azimuth pixel spacing, NOT on orbits!
    'range_pixel_size'           : 2.3,       # m
    'ground_range_pixel_size'    : 4.1,       # m
    'IW1' : {'range_resolution' : 2.7, 'azimuth_resolution': 22.5},
    'IW2' : {'range_resolution' : 3.1, 'azimuth_resolution': 22.7},
    'IW3' : {'range_resolution' : 3.5, 'azimuth_resolution': 22.6},
    'noise_equivalent_sigma_zero': -22,       # dB
}


##--------------------  S-band  --------------------##

# HJ-1C (Huan Jing-1C)
# launch date: 2012-11-19
# end    date: operational
# https://www.eoportal.org/satellite-missions/hj-1
# Liu et al. (2014, J Radar), doi: 10.3724/SP.J.1300.2013.13050
# Zhang et al. (2014, J Radar), doi: https://doi.org/10.3724/SP.J.1300.2014.13135
# spatial resolution: 10 m (4 looks)
# swath width: 100 km
HJ1C = {
    # orbit
    'altitude'                   : 502e3,     # m
    'orbit_inclination'          : 97.3,      # deg
    'repeat_cycle'               : 31,        # day
    # sar / antenna
    'carrier_frequency'          : 3.13e9,    # Hz
    'pulse_repetition_frequency' : 2600,      # Hz, 2600-3700
    'chirp_bandwidth'            : 60.0e6,    # Hz
    'noise_equivalent_sigma_zero': -22,       # dB
}

# NISAR S-band
# https://nisar.jpl.nasa.gov/system/documents/files/26_NISAR_FINAL_9-6-19.pdf
NISAR_S = {
    # orbit
    'altitude'                   : 747e3,     # m
    'orbit_inclination'          : 98.4,      # degree
    'repeat_cycle'               : 12,        # day
    # sar / antenna
    'carrier_frequency'          : 3.2e9,     # Hz
    'antenna_radius'             : 12,        # m
    'pulse_repetition_frequency' : 2200,      # Hz
    'chirp_bandwidth'            : 75e6,      # Hz
    'azimuth_resolution'         : 7,         # m
    'noise_equivalent_sigma_zero': -20,       # dB (threshold), -25 dB (baseline)
}


##--------------------  L-band  --------------------##

# Seasat
# launch date: 1978-06-27
# end    date: 1978-10-10
# from Table 6-1 in Kim and Jordan (2006)
# https://www.eoportal.org/satellite-missions/seasat
SEASAT = {
    # orbit
    'altitude'                   : 787e3,     # m, mean value, 775-799 km
    'orbit_inclination'          : 108.0,     # deg
    'repeat_cycle'               : 17,        # day
    # sar / antenna
    'carrier_frequency'          : 1.275e9,   # Hz
    'antenna_length'             : 10.74,     # m
    'antenna_width'              : 2.16,      # m
    'pulse_repetition_frequency' : 1555,      # Hz, 1463-1647
    'chirp_bandwidth'            : 19e6,      # Hz
}

# JERS-1
# launch date: 1992-02-11
# end    date: 1998-10-12
# from Table 3 in Jung et al. (2014)
# https://www.eoportal.org/satellite-missions/jers-1
# https://www.eorc.jaxa.jp/ALOS/en/jers-1/sensor/sar_e.htm
# swath_width: 75 km
JERS = {
    # orbit
    'altitude'                   : 568e3,     # m, mean value
    'orbit_inclination'          : 97.7,      # deg
    'repeat_cycle'               : 44,        # day
    # sar / antenna
    'carrier_frequency'          : 1.275e9,   # Hz
    'antenna_length'             : 11.92,     # m
    'antenna_width'              : 2.2,       # m
    'doppler_bandwidth'          : 1157,      # Hz
    'pulse_repetition_frequency' : 1600,      # Hz, 1505-1606
    'chirp_bandwidth'            : 15.00e6,   # Hz
    'sampling_frequency'         : 17.10e6,   # Hz
    'azimuth_pixel_size'         : 4.3,       # m
    'ground_range_pixel_size'    : 13.9,      # m
}

# ALOS PALSAR FBS (fine beam single polarization) mode
# launch date: 2006-01-24
# end    date: 2011-04-22
# from Table 3 in Jung et al. (2014)
# https://www.eorc.jaxa.jp/ALOS/en/alos/a1_about_e.htm
# https://www.eorc.jaxa.jp/ALOS/en/alos/sensor/palsar_e.htm
ALOS = {
    # orbit
    'altitude'                   : 691.65e3,  # m, mean value
    'orbit_inclination'          : 98.16,     # deg
    'repeat_cycle'               : 46,        # day, sub cycle: 2 days
    # sar / antenna
    'carrier_frequency'          : 1.270e9,   # Hz
    'antenna_length'             : 8.9,       # m
    'antenna_width'              : 3.1,       # m
    'doppler_bandwidth'          : 1700,      # Hz
    'pulse_repetition_frequency' : 2160,      # Hz
    'chirp_bandwidth'            : 28.00e6,   # Hz, 14/28
    'sampling_frequency'         : 32.00e6,   # Hz
    'azimuth_pixel_size'         : 3.5,       # m
    'ground_range_pixel_size'    : 7.4,       # m
    'range_pixel_size' : {
        'stripmap_FBD' : 9.37,                # m
        'stripmap_FBS' : 4.68,                # m
    },
    'noise_equivalent_sigma_zero': -23,       # dB, -23/-25 dB for swath width 70/60km
}

# ALOS-2 PALSAR-2 stripmap ultra-fine single polarization mode
# launch date: 2014-05-24
# end    date: operational
# from Table 3 in Jung et al. (2014) and eoPortal Table 10-11.
# eoPortal: https://www.eoportal.org/satellite-missions/alos-2
#   Parameter       Spotlight               Stripmap            ScanSAR
#                               ultrafine / high-     / fine
#                                           sensitive
#   Bandwidth (MHz)     84          84          42      28      14
#   Ground reso (m)     3 x 1       3           6       10      100
#   Swath (km)          25 x 25     50          50      70      350
#                                            (FP:30)  (FP:30) (5 looks)
# Notes: FP for full polarization
ALOS2 = {
    # orbit
    'altitude'                   : 628e3,     # m, mean value
    'orbit_inclination'          : 97.9,      # deg
    'repeat_cycle'               : 14,        # day
    # sar / antenna
    'carrier_frequency'          : 1.2575e9,  # Hz
    'antenna_length'             : 9.9,       # m
    'antenna_width'              : 2.9,       # m
    'doppler_bandwidth'          : 1515,      # Hz
    'pulse_repetition_frequency' : 2000,      # Hz
    'chirp_bandwidth'            : 84.0e6,    # Hz
    'sampling_frequency'         : 100.0e6,   # Hz
    'azimuth_pixel_size'         : 3.8,       # m
    'ground_range_pixel_size'    : 2.4,       # m
    'range_pixel_size' : {
        'stripmap_ultrafine'     : 1.43,      # m
        'stripmap_highsensitive' : 2.86,      # m
        'scansar_normal'         : 8.58,      # m
    },
    'noise_equivalent_sigma_zero': -24,       # dB, -24/-26/-28
}

# ALOS-4 PALSAR-3
# launch date: 2024-07-01
# end    date: operational
# https://www.eorc.jaxa.jp/ALOS/en/alos-4/a4_about_e.htm
# https://www.eorc.jaxa.jp/ALOS/en/alos-4/a4_sensor_e.htm
# using stripmap 200km@3m mode as reference
ALOS4 = {
    # orbit (same as ALOS-2)
    'altitude'                   : 628e3,     # m, mean value
    'orbit_inclination'          : 97.9,      # deg
    'repeat_cycle'               : 14,        # day, (15-3/14 rev/day)
    # sar / antenna
    'carrier_frequency'          : 1257.5e6,  # Hz (spotlight, 3m SM), 1236.5/1257.5/1278.5 MHz
    'chirp_bandwidth'            : 84e6,      # Hz, 84/42/28
    'range_resolution'           : 3,         # m
    'noise_equivalent_sigma_zero': -20,       # dB, -20/-24/-28
}

# SAOCOM-1A/B stripmap
# launch date: 2018-10-08 (SAOCOM-1A), 2020-08-30 (SAOCOM-1B)
# end    date: operational
# from Giudici et al. (2017) and
# https://directory.eoportal.org/web/eoportal/satellite-missions/s/saocom
# https://www.eoportal.org/satellite-missions/saocom
SAOCOM = {
    # orbit
    'altitude'                   : 619.6e3,   # m, mean value
    'orbit_inclination'          : 97.86,     # deg
    'repeat_cycle'               : 16,        # day, single satellite
    # sar / antenna
    'carrrier_frequency'         : 1.27414e9, # Hz
    'antenna_length'             : 10,        # m
    'pulse_repetition_frequency' : 4545,      # Hz
    'sampling_frequency'         : 50.0e6,    # Hz
    'noise_equivalent_sigma_zero': -28,       # dB (single/dual pol modes), -34 dB for quad pol modes
}

# LuTan-1 (stripmap mode)
# launch date: 2022-01-26 (LT1A), 2022-02-27 (LT1B)
# end    date: operational
# Table 1 from Wang et al. (2024, GRSM) at https://doi.org/10.1109/MGRS.2024.3478761
# Table 1 from Liu et al. (2022, EUSAR) at https://ieeexplore.ieee.org/document/9944327
# preliminary version: the azimuth bandwidth/frequency/pixelsize might change
LT1 = {
    # orbit
    'altitude'                   : 607e3,     # m, mean value
    'orbit_inclination'          : 97.8,      # deg
    'repeat_cycle'               : 8,         # day, single satellite
    # sar / antenna
    'carrier_frequency'          : 1.26e9,    # Hz
    'antenna_length'             : 9.8,       # m
    'antenna_width'              : 3.4,       # m
    'doppler_bandwidth'          : 2544,      # Hz
    'pulse_repetition_frequency' : 2934,      # Hz
    'chirp_bandwidth'            : 80.0e6,    # Hz
    'azimuth_pixel_size'         : 2.35,      # m
    'azimuth_resolution'         : 7.15,      # m
    'range_pixel_size'           : 1.67,      # m
    'range_resolution'           : 2.50,      # m
    'noise_equivalent_sigma_zero': -28,       # dB
}

# UAVSAR-L
# Reference:
#   https://earth.jpl.nasa.gov/estd-missions/airborne/uavsar/
#   https://airbornescience.nasa.gov/instrument/Uninhabited_Aerial_Vehicle_Synthetic_Aperture_Radar
#   https://www.eoportal.org/other-space-activities/uavsar
#   Fore et al. (2015) at https://doi.org/10.1109/TGRS.2014.2377637
#   Hu et al. (2020) at https://doi.org/10.1038/s41467-020-16617-7
# Parameters:
#   swath width = 16e3 m
#   fly path accuracy <= 10 m
#   instrument power = 3.1e3 W
#   noise equivalent sigma zero < -50 dB
#   operating altitude range = 2e3 - 18e3 m
#   ground speed range = 100 - 250 m/s
UAV_L = {
    # orbit
    'altitude'                   : 13.8e3,    # m, 2 - 18 km
    # sar / antenna
    'carrier_frequency'          : 1.2575e9,  # Hz, 23.70 cm wavelength
    'antenna_length'             : 1.6,       # m
    'antenna_width'              : 0.5,       # m
    'chirp_bandwidth'            : 80e6,      # Hz
    'range_resolution'           : 1.8,       # m
    'range_pixel_size'           : 1.67,      # m
    'azimuth_resolution'         : 0.8,       # m
    'azimuth_pixel_size'         : 0.6,       # m
    'noise_equivalent_sigma_zero': -50,       # dB
}

# NISAR
# https://nisar.jpl.nasa.gov/system/documents/files/26_NISAR_FINAL_9-6-19.pdf
NISAR_L = {
    # orbit
    'altitude'                   : 747e3,     # m, mean value
    'orbit_inclination'          : 98.4,      # deg
    'repeat_cycle'               : 12,        # day
    # sar / antenna
    'carrier_frequency'          : 1.257e9,   # Hz
    'antenna_length'             : 12,        # m
    'pulse_repetition_frequency' : 1650,      # Hz
    'chirp_bandwidth'            : 80.0e6,    # Hz
    'range_pixel_size' : {
        '24MHz'                  : 6.25,      # m
        '44MHz'                  : 3.41,      # m
        '80MHz'                  : 1.87,      # m
    },
    'noise_equivalent_sigma_zero': -25,       # dB
}


##--------------------  P-band  --------------------##

# Biomass
# launch date: 2025-04-29
# end    date: operational
# https://www.eoportal.org/satellite-missions/biomass
# Zhu et al. (2024) at https://doi.org/10.13203/j.whugis20240220
# swath width ~= 50e3  # m
BIOMASS = {
    # orbit
    'altitude'                   : 666e3,     # m, mean value
    'orbit_inclination'          : 98,        # deg
    'repeat_cycle'               : 3,         # day
    # sar / antenna
    'carrier_frequency'          : 435e6,     # Hz, 70 cm wavelength
    'antenna_length'             : 12,        # m
    'pulse_repetition_frequency' : 2000,      # Hz, 2000~4000
    'chirp_bandwidth'            : 6e6,       # Hz
    'range_pixel_size'           : 25,        # m
    'azimuth_pixel_size'         : 8,         # m
    'noise_equivalent_sigma_zero': -27,       # dB (threshold), -30 dB (goal)
}



SENSOR_DICT = {
    # X-band
    'tsx'   : TSX,
    'csk'   : CSK,
    'ksat5' : KSAT5,
    'iceye' : ICEYE,
    # C-band
    'ers'   : ERS,
    'env'   : ENV,
    'sen'   : SEN,
    'rsat1' : RSAT1,
    'rsat2' : RSAT2,
    'rcm'   : RCM,
    'gf3'   : GF3,
    # S-band
    'hj1c'  : HJ1C,
    # L-band
    'jers'  : JERS,
    'alos'  : ALOS,
    'alos2' : ALOS2,
    'alos4' : ALOS4,
    'lt1'   : LT1,
    'uav'   : UAV_L,
    'ni'    : NISAR_L,
    # P-band
    'bio'   : BIOMASS,
}
