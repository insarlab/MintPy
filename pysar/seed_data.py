#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Dec 2015: Add find_lat_lon(), add 'ref_lat','ref_lon' for geocoded file
# Yunjun, Jan 2016: Add input option for template file
# Yunjun, Apr 2016: Add maskFile input option
# Yunjun, Jun 2016: Add seed_attributes(), support to all file types
#                   Add reference file option


import os
import sys
import argparse

import h5py
import matplotlib.pyplot as plt
import numpy as np
import random
import multiprocessing
from joblib import Parallel, delayed

import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._pysar_utilities as ut
import pysar.subset as subset
from pysar._readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file


########################################## Sub Functions #############################################
###############################################################
###############################################################
def nearest(x, tbase,xstep):
    ## """ find nearest neighbour """
    dist = np.sqrt((tbase -x)**2)
    if min(dist) <= np.abs(xstep):
        indx=dist==min(dist)
    else:
        indx=[]
    return indx


###############################################################
def seed_file_reference_value(File, outName, refList, ref_y='', ref_x=''):
    ## Seed Input File with reference value in refList
    print 'Reference value: '
    print refList

    #####  IO Info
    atr = readfile.read_attribute(File)
    k = atr['FILE_TYPE']
    print 'file type: '+k

    ##### Multiple Dataset File
    if k in ['timeseries','interferograms','wrapped','coherence']:
        ##### Input File Info
        h5file = h5py.File(File,'r')
        epochList = sorted(h5file[k].keys())
        epochNum  = len(epochList)
        print 'number of epochs: '+str(epochNum)
        
        ##### Check Epoch Number
        if not epochNum == len(refList):
            print '\nERROR: Reference value has different epoch number'+\
                  'from input file.'
            print 'Reference List epoch number: '+str(refList)
            print 'Input file     epoch number: '+str(epochNum)
            sys.exit(1)
  
        ##### Output File Info
        h5out = h5py.File(outName,'w')
        group = h5out.create_group(k)
        print 'writing >>> '+outName

    ## Loop
    if k == 'timeseries':
        for i in range(epochNum):
            epoch = epochList[i]
            print epoch
            data = h5file[k].get(epoch)[:]
            
            data -= refList[i]
  
            dset = group.create_dataset(epoch, data=data, compression='gzip')

        atr  = seed_attributes(atr,ref_x,ref_y)
        for key,value in atr.iteritems():   group.attrs[key] = value

    elif k in ['interferograms','wrapped','coherence']:
        for i in range(epochNum):
            epoch = epochList[i]
            #print epoch
            data = h5file[k][epoch].get(epoch)[:]
            atr  = h5file[k][epoch].attrs

            data -= refList[i]
            atr  = seed_attributes(atr,ref_x,ref_y)

            gg = group.create_group(epoch)
            dset = gg.create_dataset(epoch, data=data, compression='gzip')
            for key, value in atr.iteritems():    gg.attrs[key] = value

            ut.printProgress(i+1,epochNum,'seeding:',epoch)
  
    ##### Single Dataset File
    else:
        data,atr = readfile.read(File)
        data -= refList
        atr  = seed_attributes(atr,ref_x,ref_y)
        writefile.write(data,atr,outName)
  
    ##### End & Cleaning
    try:
        h5file.close()
        h5out.close()
    except: pass

    return outName


def seed_file_inps(File, inps=None, outFile=None):
    '''Seed input file with option from input namespace
    Return output file name if succeed; otherwise, return None
    '''
    # Optional inputs
    if not outFile:  outFile = 'Seeded_'+os.path.basename(File)
    if not inps:  inps = cmdLineParse([''])
    print '----------------------------------------------------'
    print 'seeding file: '+File
    
    # Get stack and mask
    stack = ut.get_file_stack(File, inps.mask_file)
    mask = ~np.isnan(stack)
    if np.nansum(mask) == 0.0:
        print '\n*****************************************************'
        print   'ERROR:'
        print   'There is no pixel that has valid phase value in all datasets.' 
        print   'Check the file!'
        print   'Seeding failed'
        sys.exit(1)
    
    # 1. Reference using global average 
    if inps.method == 'global-average':
        print '\n---------------------------------------------------------'
        print 'Automatically Seeding using Global Spatial Average Value '
        print '---------------------------------------------------------'
        print 'Calculating the global spatial average value for each epoch'+\
              ' of all valid pixels ...'
        atr = readfile.read_attribute(File)
        width = int(atr['WIDTH'])
        length = int(atr['FILE_LENGTH'])
        box = (0,0,width,length)
        meanList = ut.spatial_average(File, mask, box)
        inps.ref_y = ''
        inps.ref_x = ''
        outFile = seed_file_reference_value(File, outFile, meanList, inps.ref_y, inps.ref_x)
        return outFile
    
    # 2. Reference using specific pixel
    # 2.1 Find reference y/x
    if not inps.ref_y or not inps.ref_x:
        if inps.coherence_file:
            inps.method = 'max-coherence'
            inps.ref_y, inps.ref_x = select_max_coherence_yx(inps.coherence_file, mask)
        elif inps.method == 'random':
            inps.ref_y, inps.ref_x = random_select_reference_yx(mask)
        elif inps.method == 'manual':
            inps = manual_select_reference_yx(stack, inps)

    # 2.2 Seeding file with reference y/x
    if inps.ref_y and inps.ref_x:
        print 'Seed file with input reference y/x ...'
        if mask[inps.ref_y, inps.ref_x]:
            print 'Referencing input file to pixel in y/x: (%d, %d)'%(inps.ref_y, inps.ref_x)
            box = (inps.ref_x, inps.ref_y, inps.ref_x+1, inps.ref_y+1)
            refList = ut.spatial_average(File, mask, box)
            outFile = seed_file_reference_value(File, outFile, refList, inps.ref_y, inps.ref_x)
        else:
            print '\nInput reference y/x has NaN value in file stacking, skip seeding.'
            return None
    else:
        sys.exit('ERROR: can not find reference y/x! Seeding FAILED.')
    
    return outFile


###############################################################
def seed_attributes(atr_in,x,y):
    atr = dict()
    for key, value in atr_in.iteritems():  atr[key] = str(value)
    
    atr['ref_y']=y
    atr['ref_x']=x
    try:
        atr['X_FIRST']
        lat = subset.coord_radar2geo(y,atr,'y')
        lon = subset.coord_radar2geo(x,atr,'x')
        atr['ref_lat']=lat
        atr['ref_lon']=lon
        geocoord='yes'
    except: geocoord='no'

    return atr


###############################################################
def random_select_reference_yx(data_mat):
    print '\n---------------------------------------------------------'
    print   'Random select reference point ...'
    print   '---------------------------------------------------------'
    
    nrow,ncol = np.shape(data_mat)
    y = random.choice(range(nrow))
    x = random.choice(range(ncol))
    while data_mat[y,x] == 0:
        y = random.choice(range(nrow))
        x = random.choice(range(ncol))
    return y,x


def manual_select_reference_yx(stack, inps):
    '''
    Input: 
        data4display : 2D np.array, stack of input file
        inps    : namespace, with key 'ref_x' and 'ref_y', which will be updated
    '''
    print '\n---------------------------------------------------------'
    print   'Manual select reference point ...'
    print   'Click on a pixel that you want to choose as the refernce '
    print   '    pixel in the time-series analysis;'
    print   'Then close the displayed window to continue.'
    print   '---------------------------------------------------------'

    ## Mutable object
    ## ref_url: http://stackoverflow.com/questions/15032638/how-to-return
    #           -a-value-from-button-press-event-matplotlib
    SeedingDone = {}
    SeedingDone['key'] = 'no'

    ##### Display
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.imshow(stack)

    ##### Selecting Point
    def onclick(event):
        if event.button==1:
            print 'click'
            x = int(event.xdata+0.5)
            y = int(event.ydata+0.5)

            if not np.isnan(stack[y][x]):
                print 'valid input reference y/x: '+str([y, x])
                inps.ref_y = y
                inps.ref_x = x
                #plt.close(fig) 
            else:
                print '\nWARNING:'
                print 'The selectd pixel has NaN value in data.'
                print 'Try a difference location please.'
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()

    return inps


def select_max_coherence_yx(corFile, mask=None):
    print '\n---------------------------------------------------------'
    print   'Searching pixel with max coherence ...'
    coh, coh_atr = readfile.read(corFile)
    if not mask is None:
        coh[mask==0] = 0.0
    y, x = np.unravel_index(np.argmax(coh), coh.shape)
    print   'y/x: '+str([y, x])
    print   '---------------------------------------------------------'

    return y, x


###############################################################
def print_warning(next_method):
    print '-----------------------------------------------------'
    print 'WARNING:'
    print 'Input file is not referenced to the same pixel yet!'
    print '-----------------------------------------------------'
    print 'Continue with default automatic seeding method: '+next_method+'\n'
    return


###############################################################
def read_seed_template2inps(template_file, inps=None):
    '''Read seed/reference info from template file and update input namespace'''
    if not inps:
        inps = cmdLineParse([''])
    
    template = readfile.read_template(template_file)
    templateKeyList = template.keys()
    
    if not inps.ref_y or not inps.ref_x:
        if 'pysar.seed.yx' in templateKeyList:
            inps.ref_y, inps.ref_x = [int(i) for i in template['pysar.seed.yx'].split(',')]
        elif 'pysar.reference.yx' in templateKeyList:
            inps.ref_y, inps.ref_x = [int(i) for i in template['pysar.reference.yx'].split(',')]
        else: print 'No y/x input from template'
    
    if not inps.ref_lat or not inps.ref_lon:
        if 'pysar.seed.lalo' in templateKeyList:
            inps.ref_lat, inps.ref_lon = [float(i) for i in template['pysar.seed.lalo'].split(',')]
        elif 'pysar.reference.lalo' in templateKeyList:
            inps.ref_lat, inps.ref_lon = [float(i) for i in template['pysar.reference.lalo'].split(',')]
        else: print 'No lat/lon input from template'
    
    return inps


def read_seed_reference2inps(reference_file, inps=None):
    '''Read seed/reference info from reference file and update input namespace'''
    if not inps:
        inps = cmdLineParse([''])
    atr_ref = readfile.read_attribute(inps.reference_file)
    atr_ref_key_list = atr_ref.keys()
    if (not inps.ref_y or not inps.ref_x) and 'ref_x' in atr_ref_key_list:
        inps.ref_y = int(atr_ref['ref_y'])
        inps.ref_x = int(atr_ref['ref_x'])
    if (not inps.ref_lat or not inps.ref_lon) and 'ref_lon' in atr_ref_key_list:
        inps.ref_lat = float(atr_ref['ref_lat'])
        inps.ref_lon = float(atr_ref['ref_lon'])
    return inps


#########################################  Usage  ##############################################
def usage():
    print '''
****************************************************************************************
  Referencing all interferograms to the same pixel.

  Usage:
      seed_data.py -f filename -t templateFile               [ -M MaskFile -o out_name]
      seed_data.py -f filename -y lineNumber -x pixelNumber  [ -M MaskFile]
      seed_data.py -f filename -l latitude   -L longitude    [ -M MaskFile]
      seed_data.py -f filename --method                      [ -M MaskFile]

      -f : the interferograms, time-series or velocity file saved in hdf5 format.
      -M : Mask file 
           ##### Priority:
               Input mask file > pysar.mask.file
      -o : output file name [Seeded_file by default]
  
      Input Method:
      -y : line number  of the reference pixel
      -x : pixel number of the reference pixel
      -l : latitude     of the reference pixel (for geocoded file)
      -L : longitude    of the reference pixel (for geocoded file)
      -r : reference file, use seeding info of this file to seed input file
      -t : template file with setting of reference point information
           Example: pysar.reference.yx   = 1160,300
                    pysar.reference.lalo = 33.1,130.0
  
           Priority:
           lat/lon > y/x
           Direct input coordinate (-y/x/l/L) > reference file (-r) > template file (-t)
  
      Non-Input Method: 
      *This is effective only when there is no valid reference value from input coordinate (options above)
      --manual         : display stack of input file and manually select reference point
      --max-coherence  : automatically select point with highest coherence value as reference point [default]
      --global-average : use spatial global average value as reference value for each epoch
      --random         : random select point as reference point
      
      -c : coherence file, used in automatic seeding based on max coherence.
  
      Reference value cannot be nan, thus, all selected reference point can not be:
          a. non zero in mask, if mask is given
          b. non nan  in data (stack)

    '''
    return

TEMPLATE='''
pysar.reference.yx   = 1160,300
pysar.reference.lalo = 33.1,130.0
'''

NOTE='''note: Reference value cannot be nan, thus, all selected reference point must be:
  a. non zero in mask, if mask is given
  b. non nan  in data (stack)
'''

EXAMPLE='''example:
  seed_data.py .h5 -t ShikokuT417F650_690AlosA.template
  seed_data.py 091120_100407.unw -y 257       -x 151             -m Mask.h5
  seed_data.py geo_velocity.h5   -l 34.45     -L -116.23         -m Mask.h5
  seed_data.py timeseries.h5     -r Seeded_velocity.h5
  
  seed_data.py unwrapIfgram.h5 -c average_spatial_coherence.h5
  seed_data.py unwrapIfgram.h5 --method manual
  seed_data.py unwrapIfgram.h5 --method random
  seed_data.py timeseries.h5   --method global-average 
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Reference to the same pixel in space.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=NOTE+'\n'+EXAMPLE)
    
    parser.add_argument('file', nargs='+', help='file(s) to be referenced.')
    parser.add_argument('-m','--mask', dest='mask_file', help='mask file')
    parser.add_argument('-o', '--outfile', help='output file name, disabled when more than 1 input files.')
    parser.add_argument('--no-parallel', dest='parallel', action='store_false',\
                        help='Disable parallel processing. Diabled auto for 1 input file.\n')
    
    coord_group = parser.add_argument_group('input coordinates')
    coord_group.add_argument('-y','--row', dest='ref_y', type=int, help='row/azimuth  number of reference pixel')
    coord_group.add_argument('-x','--col', dest='ref_x', type=int, help='column/range number of reference pixel')
    coord_group.add_argument('-l','--lat', dest='ref_lat', type=float, help='latitude  of reference pixel')
    coord_group.add_argument('-L','--lon', dest='ref_lon', type=float, help='longitude of reference pixel')
    
    coord_group.add_argument('-r','--reference', dest='reference_file', help='use reference/seed info of this file')
    coord_group.add_argument('-t','-template', dest='template_file',\
                             help='template with reference info as below:\n'+TEMPLATE)

    parser.add_argument('-c','--coherence', dest='coherence_file',\
                        help='use input coherence file to find the pixel with max coherence for reference pixel.')
    parser.add_argument('--method', default='random',\
                        choices=['input-coord','max-coherence','manual','random','global-average'], \
                        help='method to select reference pixel:\n\n'+\
                             'input-coord   : input specific coordinates, enabled when there are coordinates input\n'+\
                             'max-coherence : select pixel with highest coherence value as reference point\n'+\
                             '                enabled when there is --coherence option input\n'+\
                             'manual        : display stack of input file and manually select reference point\n'+\
                             'random        : random select pixel as reference point\n'+\
                             'global-average: for each dataset, use its spatial average value as reference value\n'+\
                             '                reference pixel is changing for different datasets\n')
    
    inps = parser.parse_args()
    return inps


#######################################  Main Function  ########################################
def main(argv):
    
    inps = cmdLineParse()
    inps.file = ut.get_file_list(inps.file)
    
    atr = readfile.read_attribute(inps.file[0])
    length = int(atr['FILE_LENGTH'])
    width  = int(atr['WIDTH'])

    # check outfile and parallel option
    if len(inps.file) > 1:
        inps.outfile = None
    elif len(inps.file) == 1 and inps.parallel:
        inps.parallel =  False
        print 'parallel processing is diabled for one input file'

    ##### Check Input Coordinates
    # Read ref_y/x/lat/lon from reference/template
    # priority: Direct Input > Reference File > Template File
    if inps.template_file:
        print 'reading reference info from template: '+inps.template_file
        inps = read_seed_template2inps(inps.template_file, inps)
    if inps.reference_file:
        print 'reading reference info from reference: '+inps.reference_file
        inps = read_seed_reference2inps(inps.reference_file, inps)
    
    # Do not use ref_lat/lon input for file in radar-coord
    if not 'X_FIRST' in atr.keys() and (inps.ref_lat or inps.ref_lon):
        print 'Lat/lon reference input is disabled for file in radar coord.'
        inps.ref_lat = None
        inps.ref_lon = None
    
    # Convert ref_lat/lon to ref_y/x
    if inps.ref_lat and inps.ref_lon:
        inps.ref_y = subset.coord_geo2radar(inps.ref_lat, atr, 'lat')
        inps.ref_x = subset.coord_geo2radar(inps.ref_lon, atr, 'lon')
        print 'Input reference point in lat/lon: '+str([inps.ref_lat, inps.ref_lon])
    print 'Input reference point in   y/x  : '+str([inps.ref_y, inps.ref_x])
    
    if inps.ref_y and inps.ref_x:
        # Do not use ref_y/x outside of data coverage
        if not (0<= inps.ref_y <= length and 0<= inps.ref_x <= width):
            inps.ref_y = None
            inps.ref_x = None
            print 'WARNING: input reference point is OUT of data coverage!'
            print 'Continue with other method to select reference point.'
        
        # Do not use ref_y/x in masked out area
        if inps.mask_file:
            print 'mask: '+inps.mask_file
            mask = readfile.read(inps.mask_file)[0]
            if mask[inps.ref_y, inps.ref_x] == 0:
                inps.ref_y = None
                inps.ref_x = None
                print 'WARNING: input reference point is in masked OUT area!'
                print 'Continue with other method to select reference point.'
    
    ##### Select method
    if inps.ref_y and inps.ref_x:
        inps.method = 'input-coord'
    elif inps.coherence_file:
        if os.path.isfile(inps.coherence_file):
            inps.method = 'max-coherence'
        else: 
            inps.coherence_file = None
    
    if inps.method == 'manual':
        inps.parallel = False
        print 'Parallel processing is disabled for manual seeding method.'

    ##### Seeding file by file
    if inps.parallel:
        num_cores = multiprocessing.cpu_count()
        print 'parallel processing using %d cores ...'%(num_cores)
        Parallel(n_jobs=num_cores)(delayed(seed_file_inps)(file, inps) for file in inps.file)
    else:
        for file in inps.file:
            print '-------------------------------------------'
            seed_file_inps(inps.file[0], inps, inps.outfile)

    print 'Done.'
    return


################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])


