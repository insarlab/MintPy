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
    'csk'   : ['csk', 'csk1', 'csk2', 'csk3', 'csk4', 'cos', 'cosmo', 'cosmoskymed'],
    'env'   : ['env', 'envisat', 'asar'],
    'ers'   : ['ers', 'ers1', 'ers2', 'ers12'],
    'gf3'   : ['gfen3', 'gaofen3', 'g3', 'gaofen'],
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


##--------------------  X-band  --------------------##
# TerraSAR-X stripmap mode in single polarization
# from Table 1 in Jung et al. (2014)
TSX = {
    'carrier_frequency'          : 9.65e9,    # Hz
    'altitude'                   : 516e3,     # m, mean value
    'antenna_length'             : 4.8,       # m
    'doppler_bandwidth'          : 2770,      # Hz
    'pulse_repetition_frequency' : 3800,      # Hz
    'chirp_bandwidth'            : 100e6,     # Hz
    'sampling_frequency'         : 109.89e6,  # Hz
    'azimuth_pixel_size'         : 2.0,       # m
    'ground_range_pixel_size'    : 2.2,       # m
}

# COSMO-SkyMed stripmap HIMAGE mode
# from Table 1 in Jung et al. (2014)
CSK = {
    'carrier_frequency'          : 9.6e9,     # Hz
    'altitude'                   : 619e3,     # m, mean value
    'antenna_length'             : 5.7,       # m
    'doppler_bandwidth'          : 2670,      # Hz
    'pulse_repetition_frequency' : 3000,      # Hz
    'chirp_bandwidth'            : 117e6,     # Hz
    'sampling_frequency'         : 146.25e6,  # Hz
    'azimuth_pixel_size'         : 2.4,       # m
    'ground_range_pixel_size'    : 1.6,       # m
}

# Kompsat-5 stripmap mode
# from Table 1 in Jung et al. (2014)
KSAT5 = {
    'carrier_frequency'          : 9.66e9,    # Hz
    'altitude'                   : 550e3,     # m, mean value
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
ICEYE = {
    'carrier_frequency'          : 9.65e9,          # Hz
    'altitude'                   : 570e3,           # m
    'antenna_length'             : 3.2,             # m
    'antenna_width'              : 0.4,             # m
    'pulse_repetition_frequency' : [2e3, 10e3],     # Hz
    'chirp_bandwidth'            : [37.6e6, 299e6], # Hz
}


##--------------------  C-band  --------------------##

# ERS-1/2
# from Table 2 in Jung et al. (2014)
# from Imaging Radar class by Howard Zebker, 2021.
ERS = {
    'carrier_frequency'          : 5.300e9,     # Hz
    'altitude'                   : 783e3,       # m, mean value
    'antenna_length'             : 10.0,        # m
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
# from Table 2 in Jung et al. (2014)
ENV = {
    'carrier_frequency'          : 5.331e9,   # Hz
    'altitude'                   : 800e3,     # m, mean value
    'antenna_length'             : 10.0,      # m
    'doppler_bandwidth'          : 1500,      # Hz
    'pulse_repetition_frequency' : 1650,      # Hz
    'chirp_bandwidth'            : 16.00e6,   # Hz
    'sampling_frequency'         : 18.00e6,   # Hz
    'azimuth_pixel_size'         : 4.3,       # m
    'ground_range_pixel_size'    : 21.3,      # m
}

# Radarsat-2 stripmap ultra-fine mode
# from Table 2 in Jung et al. (2014)
RSAT2 = {
    'carrier_frequency'          : 5.405e9,   # Hz
    'altitude'                   : 798e3,     # m, mean value
    'antenna_length'             : 6.55,      # m
    'doppler_bandwidth'          : 2308,      # Hz
    'pulse_repetition_frequency' : 3637,      # Hz
    'chirp_bandwidth'            : 78.16e6,   # Hz
    'sampling_frequency'         : 112.68e6,  # Hz
    'azimuth_pixel_size'         : 2.2,       # m
    'ground_range_pixel_size'    : 2.1,       # m
}

# GaoFen-3
# Table 2 & 6 in https://directory.eoportal.org/web/eoportal/satellite-missions/g/gaofen-3
GF3 = {
    'carrier_frequency'          : 5.4e9,     # Hz
    'altitude'                   : 755e3,     # m
    'antenna_length'             : 15,        # m
    'sampling_frequency'         : 533.33e6,  # Hz
}

# Sentinel-1 Interferometric Wide (IW / TOPS) swath mode
# Typical value:
# azfact = azResolution / azPixelSize = 1.46
# rgfact = rgResolution / rgPixelSize = 1.33
# reference:
#   1. Table 2 & Fig. 5d in Jung et al. (2014)
#   2. Table 3-1 & 7-5 in https://sentinel.esa.int/documents/247904/1877131/Sentinel-1-Product-Definition
SEN = {
    'carrier_frequency'          : 5.405e9,   # Hz
    'altitude'                   : 705e3,     # m, mean value
    'antenna_length'             : 12.3,      # m
    'antenna_width'              : 0.82,      # m
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
}


##--------------------  L-band  --------------------##

# Seasat
# from Table 6-1 in Kim and Jordan (2006)
SEASAT = {
    'carrier_frequency'          : 1.275e9,   # Hz
    'altitude'                   : 787e3,     # m, mean value
    'antenna_length'             : 10.74,     # m
    'antenna_width'              : 2.16,      # m
    'pulse_repetition_frequency' : 1555,      # Hz, 1463-1647
    'chirp_bandwidth'            : 19e6,      # Hz
}

# JERS-1
# from Table 3 in Jung et al. (2014)
JERS = {
    'carrier_frequency'          : 1.275e9,   # Hz
    'altitude'                   : 568e3,     # m, mean value
    'antenna_length'             : 11.92,     # m
    'doppler_bandwidth'          : 1157,      # Hz
    'pulse_repetition_frequency' : 1600,      # Hz, 1505-1606
    'chirp_bandwidth'            : 15.00e6,   # Hz
    'sampling_frequency'         : 17.10e6,   # Hz
    'azimuth_pixel_size'         : 4.3,       # m
    'ground_range_pixel_size'    : 13.9,      # m
}

# ALOS PALSAR FBS (fine beam single polarization) mode
# from Table 3 in Jung et al. (2014)
ALOS = {
    'carrier_frequency'          : 1.270e9,   # Hz
    'altitude'                   : 691.65e3,  # m, mean value
    'antenna_length'             : 8.9,       # m
    'doppler_bandwidth'          : 1700,      # Hz
    'pulse_repetition_frequency' : 2160,      # Hz
    'chirp_bandwidth'            : 28.00e6,   # Hz
    'sampling_frequency'         : 32.00e6,   # Hz
    'azimuth_pixel_size'         : 3.5,       # m
    'ground_range_pixel_size'    : 7.4,       # m
    'range_pixel_size' : {
        'stripmap_FBD' : 9.37,                # m
        'stripmap_FBS' : 4.68,                # m
    }
}

# ALOS-2 PALSAR-2 stripmap ultra-fine single polarization mode
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
    'carrier_frequency'          : 1.2575e9,  # Hz
    'altitude'                   : 628e3,     # m, mean value
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
    }
}

# SAOCOM-1A/B stripmap
# from Giudici et al. (2017) and
# https://directory.eoportal.org/web/eoportal/satellite-missions/s/saocom
SAOCOM = {
    'carrrier_frequency'         : 1.27414e9, # Hz
    'altitude'                   : 619.6e3,   # m, mean value
    'antenna_length'             : 10,        # m
    'pulse_repetition_frequency' : 4545,      # Hz
    'sampling_frequency'         : 50.0e6,    # Hz
}

# LuTan-1 (stripmap mode)
# preliminary version: the azimuth bandwidth/frequency/pixelsize might change
LT1 = {
    'carrier_frequency'          : 1.26e9,    # Hz
    'altitude'                   : 607e3,     # m, mean value
    'antenna_length'             : 9.8,       # m
    'antenna_width'              : 3.4,       # m
    'doppler_bandwidth'          : 2544,      # Hz
    'pulse_repetition_frequency' : 2934,      # Hz
    'chirp_bandwidth'            : 60.0e6,    # Hz
    'azimuth_pixel_size'         : 2.35,      # m
    'range_pixel_size'           : 1.67,      # m
    'azimuth_resolution'         : 7.15,      # m
    'range_resolution'           : 2.50,      # m
}

# NISAR
# https://nisar.jpl.nasa.gov/system/documents/files/26_NISAR_FINAL_9-6-19.pdf
NISAR_L = {
    'carrier_frequency'          : 1.257e9,   # Hz
    'altitude'                   : 747e3,     # m, mean value
    'antenna_length'             : 12,        # m
    'pulse_repetition_frequency' : 1650,      # Hz
    'chirp_bandwidth'            : 80.0e6,    # Hz
    'range_pixel_size' : {
        '24MHz'                  : 6.25,      # m
        '44MHz'                  : 3.41,      # m
        '80MHz'                  : 1.87,      # m
    }
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
    'rsat2' : RSAT2,
    # L-band
    'jers'  : JERS,
    'alos'  : ALOS,
    'alos2' : ALOS2,
    'lt1'   : LT1,
    'ni'    : NISAR_L,
}
