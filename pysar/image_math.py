#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2015, Yunjun Zhang                          #
# Author:  Yunjun Zhang                                    #
############################################################


import sys
import os

import h5py
import numpy as np

import pysar._readfile as readfile
import pysar._writefile as writefile


########################  Sub Functions  #########################
#####################  Operation  ####################
def operation(data,operator,operand):
    if   operator == 'plus':     data2 = data + operand;
    elif operator == 'minus':    data2 = data - operand;
    elif operator == 'multiply': data2 = data * operand;
    elif operator == 'divide':   data2 = data / operand;
    elif operator == 'exp':      data2 = data ^ operand;      ## not working, come back later
  
    return data2

#####################  Image Add  ####################
def add(data1,data2):
    data = data1 + data2;
    data[np.isnan(data1)] = data2[np.isnan(data1)];
    data[np.isnan(data2)] = data1[np.isnan(data2)];
  
    return data

#####################  Image Diff  ####################
def diff(data1,data2):
    data = data1 - data2;
    data[np.isnan(data2)] = data1[np.isnan(data2)];
  
    return data


#######################  Usage  ######################
def usage():
    print '''
***************************************************************
  Basic Mathmatic Operations of file and value
 
  Usage:  image_math.py file operator operand [outputFileName]

     file     : input file. Support all PySAR HDF5 and ROI_PAC files
                  PySAR HDF5 files: velocity, timeseries, interferograms, ...
                  ROI_PAC    files: .unw .cor .int .hgt .dem .trans
     operator : mathmatic operator, including: + - * / ^, other names:
                  +, plus, add, addition
                  -, minus, substract, substraction
                  *, multiply, multiplication, times
                  /, divide, division, obelus
                  ^, exp, exponential [not implemented yet]
     operand  : input value
     outName  : output file name (optional, default is file_operatorOperand)
  
  Example:
      image_math.py velocity.h5   '+' 0.5
      image_math.py velocity.h5   '-' 0.5
      image_math.py velocity.h5   '*' 1.5
      image_math.py velocity.h5   '/' 1.5
      image_math.py velocity.h5   add    0.5
      image_math.py velocity.h5   divide 0.5 velocity_divide0.5.h5

      image_math.py timeseries.h5 '+' 0.5
      image_math.py LoadedData.h5 '+' 0.5
      image_math.py temporal_coherence.h5 '+' 0.5

      image_math.py geo_080212_101120.unw '+' 0.5
      image_math.py geo_080212_101120.cor '+' 0.5

***************************************************************
    '''


########################  Main Functions  ########################
def main(argv):

    ########### Check Inputs #############
    try:
        file     = sys.argv[1]
        operator = sys.argv[2]
        operand  = float(sys.argv[3])
    except:
        usage();sys.exit(1)
  
    if   operator in ['+','plus',  'add',      'addition']:        operator = 'plus'
    elif operator in ['-','minus', 'substract','substraction']:    operator = 'minus'
    elif operator in ['*','times', 'multiply', 'multiplication']:  operator = 'multiply'
    elif operator in ['/','obelus','divide',   'division']:        operator = 'divide'
    elif operator in ['^','exp',   'exponential']:                 operator = 'exp'
    else:  print 'ERROR: Unrecognized operator: '+operator;  sys.exit(1)
    print '\n*************** Image Math ******************'
    print 'operation: '+operator+' '+str(operand)
  
    ext = os.path.splitext(file)[1]
    try:    outName = sys.argv[4]
    except: outName = file.split('.')[0]+'_'+operator+str(operand)+ext

    ########### Read - Calculate - Write  ###########
    ##### PySAR HDF5 files ######
    if ext in ['.h5','.he5']:
        try: h5file=h5py.File(file,'r')
        except: print 'ERROR: can not open file: '+file; sys.exit(1)
        k=h5file.keys()
        if 'interferograms' in k: k[0] = 'interferograms'
        elif 'coherence'    in k: k[0] = 'coherence'
        elif 'timeseries'   in k: k[0] = 'timeseries'
        print 'Input file is '+k[0]
      
        h5fileOut = h5py.File(outName,'w'); print 'writing >>> '+outName
        group = h5fileOut.create_group(k[0])
   
        if k[0] in ('velocity','temporal_coherence','rmse','mask','dem'):
            dset = h5file[k[0]].get(k[0])
            data = dset[0:dset.shape[0],0:dset.shape[1]]
       
            dataOut = operation(data,operator,operand)
       
            dset = group.create_dataset(k[0], data=dataOut, compression='gzip')
            for key , value in h5file[k[0]].attrs.iteritems():
                group.attrs[key]=value
   
        elif k[0] == 'timeseries':
            dateList = h5file[k[0]].keys();  print 'number of dates: '+str(len(dateList))
            for date in dateList:
                print date
                dset = h5file[k[0]].get(date)
                data = dset[0:dset.shape[0],0:dset.shape[1]]
       
                dataOut = operation(data,operator,operand)
       
                dset = group.create_dataset(date, data=dataOut, compression='gzip')
            for key,value in h5file[k[0]].attrs.iteritems():
                group.attrs[key] = value
   
        elif k[0] in ['interferograms','coherence','wrapped']:
            ifgramList = h5file[k[0]].keys();  print 'number of epochs: '+str(len(ifgramList))
            for igram in ifgramList:
                print igram
                dset = h5file[k[0]][igram].get(igram)
                data = dset[0:dset.shape[0],0:dset.shape[1]]
        
                dataOut = operation(data,operator,operand)
        
                group2 = group.create_group(igram)
                dset = group2.create_dataset(igram, data=dataOut, compression='gzip')
                for key, value in h5file[k[0]][igram].attrs.iteritems():
                    group2.attrs[key] = value
       
            try:
                mask = h5file['mask'].get('mask')
                gm = h5fileOut.create_group('mask')
                dset = gm.create_dataset('mask', data=mask, compression='gzip')
            except:  print 'No group for mask found in the file.'
       
            try:
                Cset = h5file['meanCoherence'].get('meanCoherence')
                gm = h5fileOut.create_group('meanCoherence')
                dset = gm.create_dataset('meanCoherence', data=Cset, compression='gzip')
            except:  print 'No group for meanCoherence found in the file'

        else: print 'ERROR: Unrecognized HDF5 file type: '+k[0]; sys.exit(1)
   
        h5file.close()
        h5fileOut.close()

    ##### ROI_PAC files #######
    elif ext in ['.unw','.cor','.hgt','.dem','.trans']:
        print 'Input file is '+ext+'\nwriting >>> '+outName
        if ext in ['.unw','.cor','.hgt']:
            a,p,r = readfile.read_float32(file)
            p2 = operation(p,operator,operand)
            writefile.write_float32(p2,outName)
        elif ext == '.dem':
            p,r = readfile.read_real_int16(file)
            p2 = operation(p,operator,operand)
            writefile.write_real_int16(p2,outName)
        elif ext == '.trans':
            a,p,r = readfile.read_float32(file)
            a2 = operation(a,operator,operand)
            p2 = operation(p,operator,operand)
            writefile.write_float32(a2,p2,outName)
   
        ## Write atrributes file
        f = open(outName+'.rsc','w')
        for k in r.keys():    f.write(k+'    '+r[k]+'\n')
        f.close()
  
  
    else: print 'ERROR: Unrecognized file extension: '+ext; sys.exit(1)
  
    print 'Done.'

#######################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])  

