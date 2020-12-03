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
