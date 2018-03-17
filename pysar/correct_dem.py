#!/usr/bin/env python3
# Author: Emre Havazli

import os
import sys
import numpy as np
import h5py
from pysar.utils import readfile, writefile



def usage():
    print('''
*******************************************
  Usage: correct_dem.py demFile geo_demErrorFile

  Example:
      correct_dem.py $DEMDIR/Socorro-30/Socorro_30.dem geo_DEM_error.h5
      correct_dem.py $DEMDIR/Socorro-30/Socorro_30.dem geo_DEM_error.h5

*******************************************         
    ''')


def main(argv):
    try:
        dem_file = argv[1]
        dem_error = argv[2]
    except: usage();  sys.exit(1)

    dem, demrsc = readfile.read_real_int16(dem_file)
    g = h5py.File(dem_error,'r')
    dset  = g['dem'].get('dem')
    dem_error = dset[0:dset.shape[0]]

    print('Correcting the DEM')
    sum = dem + dem_error
    print('Creating the new DEM')
    writefile.write_real_int16(sum,'DEM_w_error.dem')

    rsc_file = open('DEM_w_error.dem.rsc','w')
    for k in list(demrsc.keys()):
        rsc_file.write(k+'	'+demrsc[k]+'\n')
    rsc_file.close()
          
    date12_file=open('111111-222222_baseline.rsc','w')
    date12_file.write('P_BASELINE_TOP_ODR'+'     '+ '000')
    date12_file.close()

##########
if __name__ == '__main__':
    main(sys.argv[:])
