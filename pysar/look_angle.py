#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import sys
import numpy as np
import getopt
import h5py


############################################################
def usage():
    print('''usage: look_angle.py  file_radarCoord

Generating look angle for each pixel

example:
  look_angle.py  timeseries.h5
  look_angle.py  velocity.h5
    ''')
    return


############################################################
def main(argv):
    try:  opts, args = getopt.getopt(argv,"f:h")
    except getopt.GetoptError:  usage() ; sys.exit(1)
  
  
    if opts==[]:  usage() ; sys.exit(1)
    for opt,arg in opts:
        if opt in ("-h","--help"):  usage();  sys.exit()
        elif opt == '-f':  File = arg
    
        h5file=h5py.File(File,'r')
        look_n=float(h5file['velocity'].attrs['LOOK_REF1'])
        look_f=float(h5file['velocity'].attrs['LOOK_REF2'])
        length=float(h5file['velocity'].attrs['FILE_LENGTH'])
        width = float(h5file['velocity'].attrs['WIDTH'])
        look_step=(look_f-look_n)/width
        lookx=np.arange(look_n,look_f,look_step)
        
        look_angle = np.tile(lookx,(length,1))
        print(look_n)
        print((lookx[0]))
        print(look_f)
        print((lookx[-1]))
        print((lookx.shape)) 
        print(width)
        print(length)
        print((np.shape(look_angle)))
    
        h5file2 = h5py.File('look_angle.h5','w')
        group=h5file2.create_group('mask')
        dset = group.create_dataset('mask', data=look_angle, compression='gzip')
    
        for key, value in list(h5file['velocity'].attrs.items()):
              group.attrs[key] = value
    
        h5file.close()
        h5file2.close()
    return

        
############################################################  
if __name__ == '__main__':
    main(sys.argv[1:])
