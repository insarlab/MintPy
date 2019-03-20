from __future__ import print_function
import sys
import os


pysar_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(1, pysar_path)
sys.path.insert(1, os.path.join(pysar_path, 'defaults'))
sys.path.insert(1, os.path.join(pysar_path, 'objects'))
sys.path.insert(1, os.path.join(pysar_path, 'simulation'))
sys.path.insert(1, os.path.join(pysar_path, 'utils'))

from pysar.version import *
__version__ = release_version

try:
    os.environ['PYSAR_HOME']
except KeyError:
    print('Using default PySAR Path: %s' % (pysar_path))
    os.environ['PYSAR_HOME'] = pysar_path


## Modules dependency graph
#/pysar
#    # Level 0
#    /simulation
#        forward_model
#        fractal
#    /objects
#        giantobj
#        pysarobj
#        ramp
#        sensor
#    /utils
#        /solvers
#            l1
#            l1regls
#            lstl1
#        ptime
#        utils0
#
#    # Level 1 (depend on level 0)
#    /simulation
#        variance  (utils/ptime)
#    /objects
#        conncomp  (objects/ramp)
#    /utils
#        variance  (utils/ptime)
#        readfile  (objects/*)
#        writefile (objects/*, utils/readfile)
#        network   (objects/*, utils/readfile)
#        utils1    (objects/*, utils/writefile)
#
#    # Level 2 (depends on level 0 + 1)
#    /objects
#        resample     (utils/readfile)
#        coord        (utils/utils1)
#    /utils
#        utils        (objects/*, utils/coord)
#        gps          (objects/*, utils/utils)
#        plot         (objects/*, utils/utils)
#        insarobj     (objects/*, utils/utils)
#
#    # Level 3 (depends on level 0 + 1 + 2)
#    /objects
#        insar_vs_gps (objects/*, utils/{gps, plot}, simulation/*)
#    /simulation
#        simulation   (objects/*, utils/*, simulation/*, ./ifgram_inversion)
#
