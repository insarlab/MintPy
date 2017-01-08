#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Jun 2016: Add template input option
#                   Add multiple files support
# Yunjun, Aug 2016: Support multiple surfaces


import os
import sys
import getopt
import glob

import numpy as np
import h5py

import pysar._remove_surface as rm
import pysar._readfile as readfile
import pysar._writefile as writefile


######################################
def usage():
    print '''
**************************************************************************
  Remove phase ramp

  Usage:
      remove_plane.py file method [Maskfile]

      -f : input file (list) that need to remove ramp
      -s : quadratic, plane, quardatic_range, quadratic_azimiuth, plane_range, plane_azimuth
      -m : (optional) a mask file with 0 values for those pixels which are not considered in
             plane estimation.
      -t : template file
      -o : output name
      -y : subset in azimuth/row direction for multiple surface removal within one track
      --save-mask : plot mask file and save it to file.

  example:
      remove_plane.py  timeseries.h5 plane
      remove_plane.py  timeseries.h5 plane             Mask.h5
      remove_plane.py  LoadedData.h5 quadratic_range   Mask.h5

      remove_plane.py  -f timeseries.h5 -t KyushuT424F640AlosA.template

      remove_plane.py  -f 'geo_100102_*.unw'  -s plane     -m Mask_tempCoh.h5
      remove_plane.py  -f 090214_101120.unw   -s quadratic -m Mask_tempCoh.h5 -y 0,2400,2000,6843

**************************************************************************
    '''

######################################
def main(argv):

    ########################## Check Inputs ################################################
    ## Default value
    Masking   = 'no'
    save_mask = 'no'
  
    if len(sys.argv) > 4:
        try: opts, args = getopt.getopt(argv,'h:f:m:o:s:t:y:',['help','save-mask'])
        except getopt.GetoptError:  print 'Error while getting args!\n';  usage(); sys.exit(1)
  
        for opt,arg in opts:
            if   opt in ['-h','--help']:    usage(); sys.exit()
            elif opt in '-f':    File     = arg
            elif opt in '-m':    maskFile = arg
            elif opt in '-o':    outName  = arg
            elif opt in '-s':    surfType = arg.lower()
            elif opt in '-t':    templateFile = arg
            elif opt in '-y':    ysub = [int(i) for i in arg.split(',')]
            elif opt in '--save-mask'  :    save_mask = 'yes'
  
    elif len(sys.argv) in [3,4]:
        File          = argv[0]
        surfType      = argv[1].lower()
        try: maskFile = argv[2]
        except: pass
    else: usage(); sys.exit(1)
  
    print '\n*************** Phase Ramp Removal ***********************'
    ## Multiple Surfaces
    try:
        ysub
        if not len(ysub)%2 == 0:
            print 'ERROR: -y input has to have even length!'
            sys.exit(1)
        surfNum = len(ysub)/2
    except:
        surfNum = 1
    print 'phase ramp number: '+str(surfNum)
  
    ## Tempate File
    try:
        templateFile
        templateContents = readfile.read_template(templateFile)
    except: pass

    try:        surfType
    except:
        try:    surfType = templateContents['pysar.orbitError.method']
        except: surfType = 'plane'; print 'No ramp type input, use plane as default'
    print 'phase ramp type  : '+surfType
  
    ## Input File(s)
    fileList = glob.glob(File)
    fileList = sorted(fileList)
    print 'input file(s): '+str(len(fileList))
    print fileList
  
    atr = readfile.read_attribute(fileList[0])
    length = int(atr['FILE_LENGTH'])
    width  = int(atr['WIDTH'])

    ## Output File(s)
    if   len(fileList) >  1:    outName = ''
    elif len(fileList) == 0:    print 'ERROR: Cannot find input file(s)!';  sys.exit(1)
    else:    ## Customized output name only works for single file input
        try:
            outName
        except:
            ext     = os.path.splitext(fileList[0])[1].lower()
            outName = os.path.basename(fileList[0]).split(ext)[0]+'_'+surfType+ext
  
    ##### Read Mask File 
    ## Priority:
    ## Input mask file > pysar.mask.file > existed Modified_Mask.h5 > existed Mask.h5
    try:      maskFile
    except:
        try:  maskFile = templateContents['pysar.mask.file']
        except: pass
  
    try:
        Mask,Matr = readfile.read(maskFile)
        print 'mask file: '+maskFile
        Masking = 'yes'
    except:
        print 'No mask. Use the whole area for ramp estimation.'
        Masking = 'no'
        Mask=np.ones((length,width))

    ## Plot mask
    if save_mask == 'yes':
        mask_dis = np.zeros((length,width))
        if surfNum == 1:
            mask_dis = Mask
        else:
            i = 0
            mask_dis[ysub[2*i]:ysub[2*i+1],:] = Mask[ysub[2*i]:ysub[2*i+1],:]
            for i in range(1,surfNum):
                if ysub[2*i] < ysub[2*i-1]:
                    mask_dis[ysub[2*i]:ysub[2*i-1],:]  += Mask[ysub[2*i]:ysub[2*i-1],:]*(i+1)
                    mask_dis[ysub[2*i]:ysub[2*i-1],:]  /= 2
                    mask_dis[ysub[2*i-1]:ysub[2*i+1],:] = Mask[ysub[2*i-1]:ysub[2*i+1],:]*(i+1)
                else:
                    mask_dis[ysub[2*i]:ysub[2*i+1],:]   = Mask[ysub[2*i]:ysub[2*i+1],:]*(i+1)
        maskOutName = 'mask_'+str(surfNum)+surfType+'.h5'
        writefile.write(mask_dis,Matr,maskOutName)
        print 'save mask to mask_'+str(surfNum)+surfType+'.h5'

    ############################## Removing Phase Ramp #######################################
    for file in fileList:
        print '------------------------------------------'
        print 'input file : '+file
        if surfNum == 1:
            rm.remove_surface(file,surfType,Mask,outName)
        else:
            rm.remove_multiple_surface(file,surfType,Mask,ysub,outName)


###########################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])

#
