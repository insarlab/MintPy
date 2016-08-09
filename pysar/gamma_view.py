#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2015, Yunjun Zhang                          #
# Author:  Yunjun Zhang                                    #
############################################################

import sys
import os
import getopt
import pysar._readfile as readfile

def Usage():
    print '''
**************************************************************
**************************************************************
  Display GAMMA products:
     support files: .mli, .slc 

  Usage:
           gamma_view.py FILE
           gamma_view.py -f FILE -x row_subset -y col_subset -r mli_rg -a mli_az
           gamma_view.py -f FILE -l lat_subset -L lon_subset (not implemented)

           -f: (input) SLC/intensity image (FLOAT or SCOMPLEX data type)
           -x: subset in x directioin
           -y: subset in y direction
           -l: subset in latitude
           -L: subset in longitude
           -r: multilook number in range/x direction [default: 1]
           -a: multilook number in azimuth/y direction [default: 1 for .mli, 2 for .slc]
           -P: display ras file: yes or no [default: no]

  Example:
           gamma_view.py 101016.mli
           gamma_view.py -f 101016.mli -x 760:1060  -y 620:960
           gamma_view.py -f 101016.slc -x 3200:4000 -y 5500:7100 -P yes
           
**************************************************************
**************************************************************
    '''

###########################################################
def main(argv):

    disRas = 'no'
  
    if len(sys.argv)>2:
  
        try:  opts, args = getopt.getopt(argv,"h:f:x:y:l:L:r:a:o:P:")
        except getopt.GetoptError:  Usage(); sys.exit(1)
    
        for opt, arg in opts:
            if opt in ("-h","--help"):    Usage(); sys.exit()
            elif opt == '-f':   file = arg
            elif opt == '-x':   xsub = [int(i) for i in arg.split(':')];      xsub.sort()
            elif opt == '-y':   ysub = [int(i) for i in arg.split(':')];      ysub.sort()
            elif opt == '-l':   latsub = [float(i) for i in arg.split(':')];  latsub.sort()
            elif opt == '-L':   lonsub = [float(i) for i in arg.split(':')];  lonsub.sort()
            elif opt == '-r':   mli_rg = int(arg)
            elif opt == '-a':   mli_az = int(arg)
            elif opt == '-o':   outname = arg
            elif opt == '-P':   disRas = arg
    
        try:     file
        except:  Usage(); sys.exit(1)
  
    elif len(sys.argv)==2:   file = argv[0]
    else:                    Usage(); sys.exit(1)
        
    ############################################################

    ext = os.path.splitext(file)[1]
    outname='subset_'+file
  
    try:
        parContents = readfile.read_roipac_rsc(file + '.rsc')
        width  = int(parContents['WIDTH'])
        length = int(parContents['FILE_LENGTH'])
    except:
        parContents = readfile.read_par_file(file + '.par')
        width  = int(parContents['range_samples:'])
        length = int(parContents['azimuth_lines:'])
  
    # subset
    try:
        ysub
        if ysub[1] > length: ysub[1]=length;   print 'ysub[1] > length! Set ysub[1]=length='+str(length)
    except:
        ysub=[0,length]
        print 'no subset in y direction'
    try:
        xsub
        if xsub[1] > width:  xsub[1]=width;  print 'xsub[1] > width! Set xsub[1]=width='+str(width)
    except:
        xsub=[0,width]
        print 'no subset in x direction'
  
    if (ysub[1]-ysub[0])*(xsub[1]-xsub[0]) < length*width:
        subsetCmd='subset.py -f '+file+' -x '+str(xsub[0])+':'+str(xsub[1])+' -y '+str(ysub[0])+':'+str(ysub[1])+' -o '+outname
        print subsetCmd
        os.system(subsetCmd)
    else:
        outname = file
        print 'No subset.'
  
    # generate .ras file
    if ext == '.mli':
        try:    mli_rg
        except: mli_rg=1
        try:    mli_az
        except: mli_az=1
        rasCmd='raspwr '+outname+' '+str(xsub[1]-xsub[0])+' 1 0 '+str(mli_rg)+' '+str(mli_az)+' 1. .35 1 - 0'
        print rasCmd
        os.system(rasCmd)
    elif ext in ('.slc','.SLC'):
        try:    mli_rg
        except: mli_rg=1
        try:    mli_az
        except: mli_az=2
        rasCmd='rasSLC '+outname+' '+str(xsub[1]-xsub[0])+' 1 0 '+str(mli_rg)+' '+str(mli_az)+' 1. .35 1 1'
        print rasCmd
        os.system(rasCmd)
    else:
        print 'Not recognized file extension!'
        Usage(); sys.exit(1)

    # display .ras file
    if disRas in ('yes','Yes','Y','y','YES'):
        disCmd = 'display '+outname+'.ras'
        print disCmd
        os.system(disCmd)

        
############################################################
if __name__ == '__main__':
    main(sys.argv[1:])



