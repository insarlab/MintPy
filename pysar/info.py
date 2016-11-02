#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# Yunjun, Jul 2015: add 'coherence','wrapped' option, eNum
# Yunjun, Oct 2015: merge 'interferograms','coherence','wrapped' into one
# Yunjun, Oct 2016: add --tree option to check hdf5 tree/structure


import sys
import os
import getopt
import time

import h5py


############################################################
## By andrewcollette at https://github.com/h5py/h5py/issues/406
def print_attrs(name, obj):
    print name
    for key, val in obj.attrs.iteritems():
        print "    %s: %s" % (key, val)


############################################################
def Usage():
    print '''
***************************************************************
Displayes the general information of the PySAR product h5 file.

  Usage:
      info.py hdf5File  [eNum]

      file : HDF5 file, support all .h5 files
      eNum : number of interferogram/coherence in the group
             (1 as the first)
      --struct/structure/tree : show the structure tree

  Example:

      info.py timeseries.h5
      info.py velocity_demCor_masked.h5
      info.py LoadedData.h5
      info.py LoadedData.h5    3

      info.py timeseries.h5 --tree

***************************************************************
    '''

############################################################
def main(argv):

    try:    File = argv[0]
    except: Usage();sys.exit(1)
    h5file=h5py.File(File,'r')

    ## Print Structure Tree of Input HDF5 File
    try:
        argv[1] in ['--struct','--structure','--tree']
        h5file.visititems(print_attrs)
        return
    except: pass

    k=h5file.keys()
    if 'interferograms' in k: k[0] = 'interferograms'
    elif 'coherence'    in k: k[0] = 'coherence'
    elif 'timeseries'   in k: k[0] = 'timeseries'
    print '******************************************'
    print '******************************************'
    print 'PySAR'
    print '**********************'
    print 'File contains: '+ k[0]
    print '**********************'

 
    if 'timeseries' in k:
        import datetime
        from numpy import std
   
        try:
            h5file[k[0]].attrs['X_FIRST']
            print 'coordinates : GEO'
        except:
            print 'coordinates : radar'
        print '**********************'
        dateList = h5file['timeseries'].keys()
        print 'Number of epochs: '+str(len(dateList))
        print 'Start Date: '+dateList[0]
        print 'End Date: '+dateList[-1]
        print '**********************'
        print 'List of the dates:'
        print dateList
        print '**********************'
        print 'List of the dates in years'
        t=[]
        for i in range(len(dateList)):
            ti=(datetime.datetime(*time.strptime(dateList[i],"%Y%m%d")[0:5]))
            tt = ti.timetuple()
            ty=tt.tm_year + (tt.tm_mon-1)/12.0 + tt.tm_mday/365.0
            t.append(ty)
        print t
        
        print '*****************************************'
        print 'Standard deviation of aquisition times :'
        print str(std(t)) + ' years'
        print '**********************'
        print 'Attributes:'
        print''
        for key,value in h5file['timeseries'].attrs.iteritems():
            print key + ' : ' + str(value)
        print '*****************************************'
        print 'All groups in this file:'
        print h5file.keys()
        #print k


    elif k[0] in ['interferograms', 'coherence', 'wrapped']:
        ifgramList = h5file[k[0]].keys()
   
        try:
            h5file[k[0]][ifgramList[0]].attrs['X_FIRST']
            print 'coordinates : GEO'
        except:
            print 'coordinates : radar'
        print '**********************'
   
        try: 
            igramNumber=int(argv[1])
            print ifgramList[igramNumber-1]
            print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
            for key, value in h5file[k[0]][ifgramList[igramNumber-1]].attrs.iteritems():
                print key + ' : ' + value
            print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
            print ifgramList[igramNumber-1]
    
        except: 
            print 'Number of '+k[0]+': '+str(len(ifgramList)) 
            print '**********************'
            print 'List of the '+k[0]+':     eNum'
            eNum=1
            for ifg in ifgramList:
                print ifg + '    ' + str(eNum)
                eNum += 1
            print '**********************'
            print 'File contains: '+ k[0]
            print 'Number of '+k[0]+': '+str(len(ifgramList))
            print 'All groups in this file:'
            print h5file.keys()
            #print k
   
        for key, value in h5file[k[0]].attrs.iteritems():
            print key, value
  
  
    elif len(k)==1:
        try:
            h5file[k[0]].attrs['X_FIRST']
            print 'coordinates : GEO'
        except:
            print 'coordinates : radar'
        print '**********************'
        print 'Attributes:'
        print''
        for key , value in h5file[k[0]].attrs.iteritems():
            print key + ' : ' + str(value)
  
    else: print 'Unrecognized file: '+File
  
  
    print '******************************************'
    print '******************************************'
  
    h5file.close()

############################################################
if __name__ == '__main__':
    main(sys.argv[1:])



