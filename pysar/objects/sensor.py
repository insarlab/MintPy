############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2016-2018, Zhang Yunjun                     #
# Author:  Zhang Yunjun                                    #
############################################################
# Recommend import:
#   from pysar.objects import sensor


import os

sensorNames = ['ers', 'env', 'sen', 'rsat', 'rsat2', 'ksat5', 'gaofen3',
               'jers', 'alos', 'alos2', 'nisar',
               'tsx', 'csk'
              ]

# remove -_ and user lower case before standardize sensor names
standardedSensorNames = {'ers1': 'ers', 'ers2': 'ers', 'ers12': 'ers',
                         'envisat': 'env', 'asar': 'env',
                         'sentinel1': 'sen', 's1': 'sen', 's1a': 'sen', 's1b': 'sen', 's1ab': 'sen',
                         'radarsat': 'rsat', 'radarsat1': 'rsat', 'rsat1': 'rsat',
                         'radarsat2': 'rsat2',
                         'kompsat5': 'ksat5', 'kompsat': 'ksat5', 'kmps5': 'ksat5',
                         'g3': 'gaofen3', 'gaofen': 'gaofen3',
                         'jers1': 'jers',
                         'alos1': 'alos', 'palsar': 'alos', 'palsar1': 'alos',
                         'palsar2': 'alos2',
                         'terra': 'tsx', 'terrasar': 'tsx', 'terrasarx': 'tsx', 'tdx': 'tsx', 'tandemx': 'tsx',
                         'cosmo': 'csk', 'cosmosky': 'csk', 'cosmoskymed': 'csk',
                         'csk1': 'csk', 'csk2': 'csk', 'csk3': 'csk', 'csk4': 'csk'
                         }


class JERS(object):

    def __init__(self, file=None):
        self.file = file
        self.center_frequency = 1.275e9  # Hz
        self.bandwidth = 15e6            # Hz
        self.swath_width = 75e3          # meter
        self.prf = 1505.8

##################################################################


##################################################################
def azimuth_bandwidth(sensor):
    """Find the hardwired azimuth bandwidth in hertz for the given satellite"""
    if    sensor == 'Ers'  :  bandwidth =  1300.0
    elif  sensor == 'Env'  :  bandwidth =  1340.0
    elif  sensor == 'Sen'  :  bandwidth =   327.0   #IW1-327; IW2-313; IW3-314 (Yague-Martinez et al., 2016)
    elif  sensor == 'Rsat' :  bandwidth =   900.0
    elif  sensor == 'Rsat2':  bandwidth =   900.0

    elif  sensor == 'Jers' :  bandwidth =   900.0   # FA 11/2015 just copied, need to research
    elif  sensor == 'Alos' :  bandwidth =  1720.0

    elif  sensor == 'Tsx'  :  bandwidth = 15000.0
    elif  sensor == 'Csk'  :  bandwidth = 15000.0   # FA 9/2015  shoud be checked
    elif  sensor == 'Kmps5':  bandwidth = 15000.0   # shong 08/2016 sould be checked
    else: print('satellite not found'); bandwidth = None
    return bandwidth


def range_bandwidth(sensor):
    # Range Bandwidth in Hz for the given satellite
    if    sensor == 'Ers' :  bandwidth = 15.55e6
    elif  sensor == 'Env' :  bandwidth = 16.00e6
    elif  sensor == 'Sen' :  bandwidth = 56.5e6    #IW1-56.5; IW2-48.3; IW3-42.79

    elif  sensor == 'Jers':  bandwidth = 15e6      # Jers only has HH pol
    elif  sensor == 'Alos':  bandwidth = 14e6      # for FBD, 28MHz for FBS

    elif  sensor == 'Tsx' :  bandwidth = 150e6
    return bandwidth


def wavelength(sensor):
    if    sensor == 'Ers'  :  center_frequency = 5.300e9
    elif  sensor == 'Env'  :  center_frequency = 5.331e9
    elif  sensor == 'Sen'  :  center_frequency = 5.405e9
    elif  sensor == 'Rsat' :  center_frequency = 5.300e9
    elif  sensor == 'Rsat2':  center_frequency = 5.405e9

    elif  sensor == 'Jers' :  center_frequency = 1.275e9
    elif  sensor == 'Alos' :  center_frequency = 1.270e9
    elif  sensor == 'Alos2':  center_frequency = 1.270e9

    elif  sensor == 'Tsx'  :  center_frequency = 9.65e9
    elif  sensor == 'Csk'  :  center_frequency = 9.60e9
    elif  sensor == 'Kmps5':  center_frequency = 9.66e9

    c = 299792458   # m/s, speed of light
    wavelength = c / center_frequency
    return wavelength


def incidence_angle(sensor, inc_angle=None):
    if not inc_angle:
        if   sensor == 'Ers' :  inc_angle = 34.3
        elif sensor == 'Env' :  inc_angle = 34.3
        elif sensor == 'Sen' :  inc_angle = 32.9     #IW1 - 32.9; IW2 - 38.3; IW3 - 43.1 (Yague-Martinez et al., 2016)

        elif sensor == 'Jers':  inc_angle = 35.21
        elif sensor == 'Alos':  inc_angle = 34.3     # degree, for ALOS PALSAR Fine mode

        elif sensor == 'Tsx' :  inc_angle = 39.23    # Yunjun 5/2016, for TaizhouTsx, not sure it's for all cases.
    return inc_angle


def signal2noise_ratio(sensor):
    """Fine the Signal to Noise Ratio in dB for the given satellite
    Reference:
        ERS - Zebker et al., 1994, TGRS
        Envisat - Guarnieri, A.M., 2013. Introduction to RADAR. POLIMI DEI, Milano.
        JERS - https://directory.eoportal.org/web/eoportal/satellite-missions/j/jers-1
        Sentinel-1 - https://sentinels.copernicus.eu/web/sentinel/user-guides/sentinel-1-sar/acquisition-modes/interferometric-wide-swath
    """
    if   sensor.startswith('Ers') :  SNR = 11.7
    elif sensor.startswith('Env') :  SNR = 19.5
    elif sensor.startswith('Jers'):  SNR = 14
    elif sensor.startswith('Sen'):     SNR = 22
    else: print('satellite not found'); SNR = None
    return SNR


def project_name2sensor_name(project_names):
    """Get sensor name from project_name or path
    Parameters: project_names : str or list of str, name or path of template file containing project name
    Returns:    sensor : str, SAR sensor name
                project_name : str, project name
    Examples:   ('Sen',
                 'AlcedoSenDT128') = project_name2sensor_name('AlcedoSenDT128')
                ('Env',
                 'GalapagosEnvA2T061') = project_name2sensor_name('/Users/yunjunz/insarlab/GalapagosEnvA2T061/'+
                                                                  'PYSAR/GalapagosEnvA2T061.template')
    """
    sensor = None
    project_name = None

    if isinstance(project_names, str):
        project_names = [project_names]
    project_names = [p for p in project_names if p != None]        

    # get project_name if input is the path of template file
    for p in project_names:
        if any(s in p.lower() for s in sensorNames) and project_name is None:
            for i in p.split('/'):
                if any(s in i.lower() for s in sensorNames):
                    project_name = os.path.splitext(i)[0]

    # project_name --> sensor
    if project_name:
        sensor = [s.capitalize() for s in sensorNames
                  if s.lower() in project_name.lower()]
        if len(sensor) > 1:
            sensor = [s for s in sensor if s in project_name][0]
        elif len(sensor) == 1:
            sensor = sensor[0]
        else:
            msg = "No sensor name found in project_name: {}\n".format(project_name)
            msg += "Available sensor names: {}".format(sensorNames)
            raise ValueError(msg)

    return sensor, project_name



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
    if value.startswith('ers'):
        mission_name = 'ERS'
    elif value.startswith(('env', 'asar')):
        mission_name = 'ENV'
    elif value.startswith(('s1', 'sen')):
        mission_name = 'S1'
    elif value.startswith(('rs', 'rsat', 'radarsat')):
        mission_name = 'RS'
        if value.endswith('1'):
            mission_name += '1'
        else:
            mission_name += '2'
    elif value.startswith(('csk', 'cos')):
        mission_name = 'CSK'
    elif value.startswith(('tsx', 'tdx', 'terra', 'tandem')):
        mission_name = 'TSX'
    elif value.startswith('jers'):
        mission_name = 'JERS'
    elif value.startswith(('alos', 'palsar')):
        if value.endswith('2'):
            mission_name = 'ALOS2'
        else:
            mission_name = 'ALOS'
    else:
        print('Un-recognized PLATFORM attribute:', value)
        print('return None')
    return mission_name




