#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import sys
import os
import numpy as np
import getopt
import h5py


def Usage():
    print '''
***************************************************************
***************************************************************

  Generating look angle for each pixel

  Usage:
       
       look_angke velocity.h5

***************************************************************
***************************************************************
    '''

def main(argv):
    try:  opts, args = getopt.getopt(argv,"f:h")
    except getopt.GetoptError:  Usage() ; sys.exit(1)
  
  
    if opts==[]:  Usage() ; sys.exit(1)
    for opt,arg in opts:
        if opt in ("-h","--help"):  Usage();  sys.exit()
        elif opt == '-f':  File = arg
    
        h5file=h5py.File(File,'r')
        look_n=float(h5file['velocity'].attrs['LOOK_REF1'])
        look_f=float(h5file['velocity'].attrs['LOOK_REF2'])
        length=float(h5file['velocity'].attrs['FILE_LENGTH'])
        width = float(h5file['velocity'].attrs['WIDTH'])
        look_step=(look_f-look_n)/width
        lookx=np.arange(look_n,look_f,look_step)
        
        look_angle = np.tile(lookx,(length,1))
        print look_n
        print lookx[0]
        print look_f
        print lookx[-1]
        print lookx.shape 
        print width
        print length
        print np.shape(look_angle)
    
        h5file2 = h5py.File('look_angle.h5','w')
        group=h5file2.create_group('mask')
        dset = group.create_dataset('mask', data=look_angle, compression='gzip')
    
        for key, value in h5file['velocity'].attrs.iteritems():
              group.attrs[key] = value
    
        h5file.close()
        h5file2.close()

        
        
if __name__ == '__main__':
    main(sys.argv[1:])





