from __future__ import print_function
import sys
import os


mintpy_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(1, mintpy_path)
sys.path.insert(1, os.path.join(mintpy_path, 'defaults'))
sys.path.insert(1, os.path.join(mintpy_path, 'objects'))
sys.path.insert(1, os.path.join(mintpy_path, 'simulation'))
sys.path.insert(1, os.path.join(mintpy_path, 'utils'))

from mintpy.version import *
__version__ = release_version

try:
    os.environ['MINTPY_HOME']
except KeyError:
    print('Using default MintPy Path: %s' % (mintpy_path))
    os.environ['MINTPY_HOME'] = mintpy_path


## Modules dependency graph
#/mintpy
#    # Level 0
#    /simulation
#        forward_model
#        fractal
#    /objects
#        giant
#        stack
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
#        stackDict    (objects/*, utils/utils)
#
#    # Level 3 (depends on level 0 + 1 + 2)
#    /objects
#        insar_vs_gps (objects/*, utils/{gps, plot}, simulation/*)
#    /simulation
#        simulation   (objects/*, utils/*, simulation/*, ./ifgram_inversion)
#
