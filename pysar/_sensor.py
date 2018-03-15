#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2016, Yunjun Zhang                          #
# Author:  Yunjun Zhang                                    #
############################################################
# Recommended Usage:
#   import pysar._sensor as sensor
#


class JERS(object):

    def __init__(self, file=None):
        self.file = file
        self.center_frequency = 1.275e9  # Hz
        self.bandwidth = 15e6            # Hz
        self.swath_width = 75e3          # meter
        self.prf         = 1505.8
        
    
