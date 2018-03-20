# Author: Zhang Yunjun, Heresh Fattahi

import os, sys, glob
import h5py
import numpy as np
from pysar.utils import readfile, writefile


class geometry:
    '''
    Geometry object for Lat, Lon, Heigt, Incidence, Heading, Bperp, ... from the same platform and track.
    Attributes are saved in the root level.
    It contains a "geometry" group and some datasets: date12, bperp, unwrapIfgram, coherence, ...
    
    /                    Root level
    Attributes           Dictionary for metadata
    geometry             group name
        /height          2D array of float32 in size of (l, w) in meter.
        /incidenceAngle  2D array of float32 in size of (l, w) in degree.
        /latitude        2D array of float32 in size of (l, w) in degree.
        /longitude       2D array of float32 in size of (l, w) in degree.
        /rangeCoord      2D array of float32 in size of (l, w).
        /azimuthCoord    2D array of float32 in size of (l, w).
        /headingAngle    2D array of float32 in size of (l, w) in degree. (optional) 
        /shadowMask      2D array of bool in size of (l, w). (optional)
        /waterMask       2D array of bool in size of (l, w). (optional)
        /bperp           3D array of float32 in size of (n, lc, wc) in meter (optional)
        ...
    '''
    def __init__(self, file=None):
        self.file = file
    


    def close(self):
        print('close geometry: {}'.format(os.path.basename(self.file)))
        self.h5.close()
