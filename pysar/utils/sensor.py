############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2016, Yunjun Zhang                          #
# Author:  Yunjun Zhang                                    #
############################################################
# Recommended Usage:
#   import pysar.utils.sensor as sensor
#

sensors=['ers','env','sen','rsat','rsat2','ksat5','gaofen3',
         'jers','alos','alos2','nisar',
         'tsx','csk'
        ]

#remove -_ and user lower case before standardize sensor names
standardedSensorNames={'ers1':'ers','ers2':'ers','ers12':'ers',
                       'envisat':'env','asar':'env',
                       'sentinel1':'sen','s1':'sen','s1a':'sen','s1b':'sen','s1ab':'sen',
                       'radarsat':'rsat','radarsat1':'rsat','rsat1':'rsat',
                       'radarsat2':'rsat2',
                       'kompsat5':'ksat5','kompsat':'ksat5','kmps5':'ksat5',
                       'g3':'gaofen3','gaofen':'gaofen3',
                       'jers1':'jers',
                       'alos1':'alos','palsar':'alos','palsar1':'alos',
                       'palsar2':'alos2',
                       'terra':'tsx','terrasar':'tsx','terrasarx':'tsx','tdx':'tsx','tandemx':'tsx',
                       'cosmo':'csk','cosmosky':'csk','cosmoskymed':'csk',
                       'csk1':'csk','csk2':'csk','csk3':'csk','csk4':'csk'
                      }


class JERS(object):

    def __init__(self, file=None):
        self.file = file
        self.center_frequency = 1.275e9  # Hz
        self.bandwidth = 15e6            # Hz
        self.swath_width = 75e3          # meter
        self.prf         = 1505.8
        
    
