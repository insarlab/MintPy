#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Copyright(c) 2013, Zhang Yunjun, Heresh Fattahi          
# Author:  Zhang Yunjun, Heresh Fattahi                    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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


########################################################################
miami_path = True              # Package-wide variable, Auto setting for file structure of Univ. of Miami
                               # change it to False if you are not using it.
parallel_num = 8               # max core number used in parallel processing
figsize_single_min = 6.0       # default min size in inch, for single plot
figsize_single_max = 12.0      # default min size in inch, for single plot
figsize_multi = [15.0, 8.0]    # default size in inch, for multiple subplots


