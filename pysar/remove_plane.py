#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Jun 2016: Add template input option
#                   Add multiple files support
# Yunjun, Aug 2016: Support multiple surfaces


import os
import sys
import argparse

import h5py
import numpy as np

import pysar._pysar_utilities as ut
import pysar._remove_surface as rm
import pysar._readfile as readfile
import pysar._writefile as writefile


######################################
EXAMPLE='''example:
  remove_plane.py  timeseries.h5      -m Mask.h5
  remove_plane.py  timeseries.h5      -m Mask.h5         -s quadratic
  remove_plane.py  090214_101120.unw  -m Mask_tempCoh.h5 -s quadratic  -y 0,2400,2000,6843
'''


def cmdLineParse():
    parser = argparse.ArgumentParser(description='Remove phase ramp',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs='+', help='File(s) for ramp removal')
    parser.add_argument('-m','--mask', dest='mask_file', default='maskTempCoh.h5',\
                        help='mask for pixels used in ramp estimation\n'+\
                             'default - maskTempCoh.h5\n'+\
                             'no - use the whole area')
    parser.add_argument('-s', dest='surface_type', default='plane', \
                        choices={'plane','quadratic','plane_range','quadratic_range','plane_azimuth','quadratic_azimuth'},\
                        help='type of surface/ramp to remove, plane by default')
    parser.add_argument('-y', dest='ysub', type=int, nargs='*',\
                        help='subset in azimuth/row direction for multiple surface removal within one track, i.e.:\n'+\
                             '0,2400,2000,6843')
    parser.add_argument('-o','--outfile', help='Output file name. Disabled when more than 1 input files')
    parser.add_argument('--no-parallel',dest='parallel',action='store_false',default=True,\
                        help='Disable parallel processing. Diabled auto for 1 input file.')

    inps = parser.parse_args()
    if inps.ysub and not len(inps.ysub)%2 == 0:
        raise Exception('ERROR: -y input has to have even length!')
    return inps


######################################
def main(argv):
    
    inps = cmdLineParse()
    inps.file = ut.get_file_list(inps.file)
    print('input file(s): '+str(len(inps.file)))
    print(inps.file)
    
    #print '\n*************** Phase Ramp Removal ***********************'
    atr = readfile.read_attribute(inps.file[0])
    length = int(atr['FILE_LENGTH'])
    width = int(atr['WIDTH'])

    # Read mask file if inputed
    if inps.mask_file == 'no':  inps.mask_file = None
    if inps.mask_file:
        try:
            mask_atr = readfile.read_attribute(inps.mask_file)
        except:
            print('Can not open mask file: '+inps.mask_file)
            inps.mask_file = None

    # Update mask for multiple surfaces
    if inps.ysub:
        # Read mask
        if not inps.mask_file:
            Mask_temp = readfile.read(inps.mask_file)[0]
            Mask = np.zeros((length, width), dtype=np.float32)
            Mask[Mask_temp!=0] = 1
        else:
            Mask = np.ones((length, width))
        
        # Update mask for multiple surface from inps.ysub
        mask_multiSurface = np.zeros((length,width), dtype=np.float32)
        surfNum = len(inps.ysub)/2
        if surfNum == 1:
            mask_multiSurface = Mask
        else:
            i = 0
            mask_multiSurface[inps.ysub[2*i]:inps.ysub[2*i+1],:] = Mask[inps.ysub[2*i]:inps.ysub[2*i+1],:]
            for i in range(1,surfNum):
                if inps.ysub[2*i] < inps.ysub[2*i-1]:
                    mask_multiSurface[inps.ysub[2*i]:inps.ysub[2*i-1],:]  += Mask[inps.ysub[2*i]:inps.ysub[2*i-1],:]*(i+1)
                    mask_multiSurface[inps.ysub[2*i]:inps.ysub[2*i-1],:]  /= 2
                    mask_multiSurface[inps.ysub[2*i-1]:inps.ysub[2*i+1],:] = Mask[inps.ysub[2*i-1]:inps.ysub[2*i+1],:]*(i+1)
                else:
                    mask_multiSurface[inps.ysub[2*i]:inps.ysub[2*i+1],:]   = Mask[inps.ysub[2*i]:inps.ysub[2*i+1],:]*(i+1)
         
        # Write updated mask for multiple surfaces into file
        outFile = 'mask_'+str(surfNum)+inps.surface_type+'.h5'
        atr['FILE_TYPE'] = 'mask'
        writefile.write(mask_multiSurface, atr, outFile)
        print('saved mask to '+outFile)

    ############################## Removing Phase Ramp #######################################
    # check outfile and parallel option
    if inps.parallel:
        num_cores, inps.parallel, Parallel, delayed = ut.check_parallel(len(inps.file))

    if len(inps.file) == 1:
        rm.remove_surface(inps.file[0], inps.surface_type, inps.mask_file, inps.outfile, inps.ysub)

    elif inps.parallel:
        #num_cores = min(multiprocessing.cpu_count(), len(inps.file))
        #print 'parallel processing using %d cores ...'%(num_cores)
        Parallel(n_jobs=num_cores)(delayed(rm.remove_surface)(file, inps.surface_type, inps.mask_file, ysub=inps.ysub)\
                                   for file in inps.file)

    else:
        for File in inps.file:
            print('------------------------------------------')
            rm.remove_surface(File, inps.surface_type, inps.mask_file, ysub=inps.ysub)
    
    print('Done.')
    return


###########################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])

