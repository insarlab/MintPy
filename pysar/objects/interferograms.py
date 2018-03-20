# Author: Zhang Yunjun, Heresh Fattahi

import os, sys, glob
import h5py
import numpy as np
from pysar.utils import readfile, writefile


class ifgramStack:
    '''
    Stack of Interferograms object for a set of interferograms and coherence from the same platform and track.
    Attributes are saved in the root level.
    It contains a "interferograms" group and some datasets: date12, bperp, unwrapIfgram, coherence, ...
    
    /                    Root level
    Attributes           Dictionary for metadata
    interferograms       group name
        /date12          1D array of string in size of (m,) in YYYYMMDD_YYYYMMDD format
        /bperp           1D array of float32 in size of (m,) in meter.
        /unwrapIfgram    3D array of float32 in size of (m, l, w) in radian.
        /coherence       3D array of float32 in size of (m, l, w).
        /connComp        3D array of int16 in size of (m, l, w). (optional)
        /wrapIfgram      3D array of float32 in size of (m, l, w) in radian. (optional)
        /rangeOffset     3D array of float32 in size of (m, l, w). (optional)
        /azimuthOffset   3D array of float32 in size of (m, l, w). (optional)
        ...
    '''
    def __init__(self, file=None):
        self.file = file
    


    def close(self):
        print('close interferograms stack: {}'.format(os.path.basename(self.file)))
        self.h5.close()
