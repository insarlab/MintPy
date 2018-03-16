#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Purpose: Python Module for InSAR Time Series Analysis    #
# Copyright(c) 2013, Heresh Fattahi, Zhang Yunjun          #
# Author:  Heresh Fattahi, Zhang Yunjun                    #
############################################################

########################################################################
miami_path = True              # Package-wide variable, Auto setting for file structure of Univ. of Miami
                               # change it to False if you are not using it.
parallel_num = 8               # max core number used in parallel processing
figsize_single_min = 6.0       # default min size in inch, for single plot
figsize_single_max = 12.0      # default min size in inch, for single plot
figsize_multi = [15.0, 8.0]    # default size in inch, for multiple subplots


###################### Do not change below this line ###################
#from . import _datetime as ptime
#from . import _readfile as readfile
#from . import _writefile as writefile
