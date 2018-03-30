#!/usr/bin/env python3

##------------------ Variables ---------------------##
timeseriesDatasetNames = ['raw',
                          'troposphericDelay',
                          'topographicResidual',
                          'ramp',
                          'displacement']

geometryDatasetNames = ['height',
                        'latitude',
                        'longitude',
                        'rangeCoord',
                        'azimuthCoord',
                        'incidenceAngle',
                        'headingAngle',
                        'slantRangeDistance',
                        'shadowMask',
                        'waterMask',
                        'commonMask',
                        'bperp']

ifgramDatasetNames = ['unwrapPhase',
                      'coherence',
                      'connectComponent',
                      'wrapPhase',
                      'iono',
                      'rangeOffset',
                      'azimuthOffset']


##------------------ Classes ---------------------##
from .pysarobj import *



