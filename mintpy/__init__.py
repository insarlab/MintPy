from __future__ import print_function
import sys
import os

# get version info
from mintpy.version import release_version, logo
__version__ = release_version
__logo__ = logo

# check environmental variable
mintpy_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(1, mintpy_path)
sys.path.insert(1, os.path.join(mintpy_path, 'defaults'))
sys.path.insert(1, os.path.join(mintpy_path, 'objects'))
sys.path.insert(1, os.path.join(mintpy_path, 'simulation'))
sys.path.insert(1, os.path.join(mintpy_path, 'utils'))

try:
    os.environ['MINTPY_HOME']
except KeyError:
    print('Using default MintPy Path: %s' % (mintpy_path))
    os.environ['MINTPY_HOME'] = mintpy_path


module_dependency_graph = """
# level N depends on level N-1, N-2, ..., 0
/mintpy
------------------ level 0 --------------------
    /simulation
        fractal
    /objects
        colors
        giant
        ramp
        sensor
        stack
    /utils
        /solvers
            l1
            l1regls
            lstl1
        ptime
        utils0

------------------ level 1 --------------------
    /simulation
        decorrelation (utils/ptime)
        variance      (utils/ptime)
        defo_model    (utils/utils0)
    /objects
        conncomp      (objects/ramp)
    /utils
        variance      (utils/ptime)
        readfile      (objects/*)
        writefile     (objects/*, utils/readfile)
        network       (objects/*, utils/readfile)
        utils1        (objects/*, utils/writefile)

------------------ level 2 --------------------
    /objects
        resample      (utils/readfile)
        coord         (utils/utils1)
    /utils
        utils         (objects/*, utils/coord)
        gps           (objects/*, utils/utils)
        plot          (objects/*, utils/utils)
        stackDict     (objects/*, utils/utils)

------------------ level 3 --------------------
    /objects
        insar_vs_gps  (objects/*, utils/{gps, plot}, simulation/*)
    /simulation
        simulation    (objects/*, utils/*, simulation/*, ./ifgram_inversion)

"""
