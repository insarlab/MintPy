############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2016, Yunjun Zhang                          #
# Author:  Yunjun Zhang                                    #
############################################################
# Recommended Usage:
#   import pysar.objects.sensor as sensor
#

import os

sensors = ['ers', 'env', 'sen', 'rsat', 'rsat2', 'ksat5', 'gaofen3',
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


def project_name2sensor(project_names):
    """Get sensor from project_name or path
    Parameters: project_name_in : str or list of str, name or path of template file containing project name
    Returns:    sensor : str, SAR sensor name
                project_name : str, project name
    Examples:   'Sen', 'AlcedoSenDT128' = project_name2sensor('AlcedoSenDT128')
                'Env', 'GalapagosEnvA2T061' = project_name2sensor('/Users/yunjunz/insarlab/Galapagos/'+
                                                                  'GalapagosEnvA2T061/PYSAR/GalapagosEnvA2T061.template')
    """
    sensor = None
    project_name = None

    if isinstance(project_names, str):
        project_names = [project_names]
    project_names = [p for p in project_names if p != None]        

    # get project_name if input is the path of template file
    for p in project_names:
        if any(s in p.lower() for s in sensors) and project_name is None:
            for i in p.split('/'):
                if any(s in i.lower() for s in sensors):
                    project_name = os.path.splitext(i)[0]

    import pdb; pdb.set_trace()

    # project_name --> sensor
    if project_name:
        sensor = [s for s in sensors if s in project_name.lower()][0].capitalize()

    return sensor, project_name
