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

module_dependency_graph = """# level N depends on level N-1, N-2, ..., 0
/mintpy
------------------ level 0 --------------------
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
    /simulation
        fractal

------------------ level 1 --------------------
    /objects
        conncomp      (objects/ramp)
    /simulation
        decorrelation (utils/ptime)
        defo_model    (utils/utils0)
        variance      (utils/ptime)
    /utils
        readfile      (objects/{stack, giant})

------------------ level 2 --------------------
    /objects
        resample      (utils/{readfile, utils0, ptime})
        coord         (utils/{readfile, utils0, utils1})
    /utils
        writefile     (objects/{stack},         utils/{readfile})
        network       (objects/{stack, sensor}, utils/{readfile})

------------------ level 3 --------------------
    /objects
        gps           (objects/{stack, coord}, utils/{ptime, utils0, readfile})
        stackDict     (objects/{stack},        utils/{ptime, utils0, readfile})
    /simulation
        simulation    (objects/{stack},        utils/{ptime, network}, simulation/{fractal, decorrelation, defo_model})
    /utils
        utils1        (objects/{stack, ramp},  utils/{ptime, utils0, readfile, writefile})

------------------ level 4 --------------------
    /utils
        plot          (objects/{stack, coord, colors}, utils/{ptime, utils0, utils1, readfile, network})
        utils         (objects/{stack, coord},         utils/{ptime, utils0, utils1, readfile})

------------------ level 5 --------------------
    /objects
        insar_vs_gps  (objects/{stack, giant}, utils/{readfile, gps, plot, utils})

"""
