############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2013, Zhang Yunjun, Heresh Fattahi          #
# Author:  Zhang Yunjun, Heresh Fattahi, 2018 Mar          #
############################################################


from __future__ import print_function
from .version import release_version, release_date
__version__ = release_version


import sys, os
pysar_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#sys.path.insert(1,pysar_path)
#sys.path.insert(1,os.path.join(pysar_path,'pysar'))
#sys.path.insert(1,os.path.join(pysar_path,'pysar/pysarobj'))
#sys.path.insert(1,os.path.join(pysar_path,'library'))

try:
    os.environ['PYSAR_HOME']
except KeyError:
    print('Using default PySAR Path: %s'%(pysar_path))
    os.environ['PYSAR_HOME'] = pysar_path


##### PySAR modules listed by relative dependecies:
## 0. Independent modules:
# pysar.objects.pysarobj
# pysar.defaults.auto_path
# pysar.utils.writefile
# pysar.utils.datetime
# pysar.utils.sensor
# 
## Level 1 dependent modules (depends on Level 0):
# pysar.utils.readfile
# pysar.utils.network
# pysar.utils.deramp
#
## Level 2 dependent modules (depends on Level 0,1):
# pysar.utils.utils
#
## Level 3 dependent modules (depends on Level 0,1,2):
# pysar.objects.insarobj
# pysar.utils.plot
# 

