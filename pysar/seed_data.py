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
import getopt

import numpy as np
import h5py

import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._pysar_utilities as ut
import pysar.subset as sub


########################################## Sub Functions #############################################
###############################################################
def random_selection(stack):
    import random
    nrow,ncol = np.shape(stack)
  
    y = random.choice(range(nrow))
    x = random.choice(range(ncol))
  
    while stack[y,x] == 0:
        y = random.choice(range(nrow))
        x = random.choice(range(ncol))

    return y,x

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
def seed_xy(File,x,y,outName=''):
    ## Seed Input File with reference on point (y,x)
    print 'Referencing input file to pixel: (%d, %d)'%(y,x)
    ##4-tuple defining the left, upper, right, and lower pixel coordinate [optional]
    box = (x,y,x+1,y+1)

    #####  IO Info
    atr = readfile.read_attributes(File)
    k = atr['FILE_TYPE']
    if outName == '':  outName = 'Seeded_'+os.path.basename(File)

    ##### Mask
    length = int(atr['FILE_LENGTH'])
    width  = int(atr['WIDTH'])
    mask = np.ones((length,width))
    
    ## Read refernce value
    refList = ut.spatial_mean(File,mask,box)

    ## Seeding
    seed_file(File,outName,refList,x,y)

    return 1


###############################################################
def seed_file(File,outName,refList,ref_x='',ref_y=''):
    ## Seed Input File with reference value in refList
    print 'Reference value: '
    print refList

    #####  IO Info
    atr = readfile.read_attributes(File)
    k = atr['FILE_TYPE']
    print 'file type: '+k

    ##### Multiple Dataset File
    if k in ['timeseries','interferograms','wrapped','coherence']:
        ##### Input File Info
        h5file = h5py.File(File,'r')
        epochList = h5file[k].keys()
        epochList = sorted(epochList)
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

    return 1

###############################################################
def seed_attributes(atr_in,x,y):
    atr = dict()
    for key, value in atr_in.iteritems():  atr[key] = str(value)
    
    atr['ref_y']=y
    atr['ref_x']=x
    try:
        atr['X_FIRST']
        lat = sub.coord_radar2geo(y,atr,'y')
        lon = sub.coord_radar2geo(x,atr,'x')
        atr['ref_lat']=lat
        atr['ref_lon']=lon
        geocoord='yes'
    except: geocoord='no'

    return atr


###############################################################
def seed_manual(File,stack,outName):
    import matplotlib.pyplot as plt

    print '\n---------------------------------------------------------'
    print   'Manual select reference point ...'
    print   'Click on a pixel that you want to choose as the refernce '
    print   '    pixel in the time-series analysis, and then close the'
    print   '    displayed window'
    print   '---------------------------------------------------------'

    ## Mutable object
    ## ref_url: http://stackoverflow.com/questions/15032638/how-to-return-a-value-from-button-press-event-matplotlib
    SeedingDone = {}
    SeedingDone['key'] = 'no'

    ##### Display
    fig = plt.figure();
    ax  = fig.add_subplot(111);
    ax.imshow(stack)

    ##### Selecting Point
    def onclick(event):
        if event.button==1:
            print 'click'
            x = int(event.xdata+0.5)
            y = int(event.ydata+0.5)

            if not stack[y][x] == 0:
                seed_xy(File,x,y,outName)
                SeedingDone['key'] = 'yes'
                plt.close(fig) 
            else:
                print '\nWARNING:'
                print 'The selectd pixel has NaN value in data.'
                print 'Try a difference location please.'
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()

    return SeedingDone['key']


###############################################################
def seed_max_coherence(File,mask,outFile,corFile=''):
    print '\n---------------------------------------------------------'
    print   'Automatically select reference point ...'
    print   '    Based on maximum coherence.'
    print   '    Input coherence file or meanCoherence group within   '
    print   '    the file is needed.'
    print   '---------------------------------------------------------'

    SeedingDone = 'no'
    
    ##### Read Coherence
    try:
        h5file = h5py.File(File,'r')
        coh = h5file['meanCoherence'].get('meanCoherence')[:]
    except:
        try:  coh, coh_atr = readfile.read(corFile)
        except: print '\nERROR: No coherence data is found!'

    try:
        coh *= mask
        print 'Searching the pixel with maximum avergae coherence'
        y,x = np.unravel_index(np.argmax(coh), coh.shape)
        seed_xy(File,x,y,outFile)
        SeedingDone = 'yes'
    except: pass

    return SeedingDone


###############################################################
def print_warning(next_method):
    print '-----------------------------------------------------'
    print 'WARNING:'
    print 'Input file is not referenced to the same pixel yet!'
    print '-----------------------------------------------------'
    print 'Continue with default automatic seeding method: '+next_method+'\n'


#########################################  Usage  ##############################################
def Usage():
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
           Example: pysar.seed.yx   = 1160,300
                    pysar.seed.lalo = 33.1,130.0
  
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

  Examples:
     seed_data.py -f LoadedData.h5 -t ShikokuT417F650_690AlosA.template
     seed_data.py -f LoadedData.h5 -t ShikokuT417F650_690AlosA.template  -m Mask.h5
     seed_data.py -f 091120_100407.h5    -y 257       -x 151             -m Mask.h5
     seed_data.py -f velocity.h5         -l 34.45     -L -116.23         -m Mask.h5
     seed_data.py -f timeseries.h5       -r Seeded_velocity.h5

     seed_data.py -f LoadedData.h5 --manual
     seed_data.py -f LoadedData.h5 --max-coherence -c average_spatial_coherence.h5
     seed_data.py -f timeseries.h5 --global-average

****************************************************************************************
    '''


#######################################  Main Function  ########################################
def main(argv):

    global method_default
    ##### Referencing methods
    method_default = 'max_coherence'
    #method = 'manual'
    #method = 'max_coherence'        ## Use phase on point with max coherence [default]
    #method = 'global_average'       ## Use Nan Mean of phase on all pixels
    #method = 'random'
    #maskFile = 'Mask.h5'

    global SeedingDone
    
    ############################## Check Inputs ##############################
    if len(sys.argv) > 2:
        try:  opts, args = getopt.getopt(argv,"h:c:f:m:y:x:l:L:t:o:r:",\
                                         ['manual','max-coherence','global-average','random'])
        except getopt.GetoptError:  Usage() ; sys.exit(1)

        for opt,arg in opts:
            if   opt in ("-h","--help"):   Usage();  sys.exit()
            elif opt == '-f':        File     = arg
            elif opt == '-m':        maskFile = arg
            elif opt == '-c':        corFile  = arg
            elif opt == '-o':        outFile  = arg

            elif opt == '-y':        ry       = int(arg)
            elif opt == '-x':        rx       = int(arg)
            elif opt == '-l':        rlat     = float(arg)
            elif opt == '-L':        rlon     = float(arg)
            elif opt == '-r':        refFile  = arg
            elif opt == '-t':        templateFile = arg

            elif opt == '--global-average' :  method = 'global_average'
            elif opt == '--manual'         :  method = 'manual'
            elif opt == '--max-coherence'  :  method = 'max_coherence'
            elif opt == '--random'         :  method = 'random'

    elif len(sys.argv)==2:
        if   argv[0]=='-h':            Usage(); sys.exit(1)
        elif os.path.isfile(argv[0]):  File = argv[0]
        else:  print 'Input file does not existed: '+argv[0];  sys.exit(1)
    elif len(sys.argv)<2:             Usage(); sys.exit(1)

    ##### Input File Info
    try:
        File
        atr = readfile.read_attributes(File)
        k = atr['FILE_TYPE']
        length = int(atr['FILE_LENGTH'])
        width  = int(atr['WIDTH'])
    except:  Usage() ; sys.exit(1)
    ext = os.path.splitext(File)[1].lower()

    try:    outFile
    except: outFile = 'Seeded_'+File
  
    ############################## Reference Point Input ####################
    try:
        refFile
        atr_ref = readfile.read_attributes(refFile)
    except: pass
  
    try:
        templateFile
        templateContents = readfile.read_template(templateFile)
    except: pass

    ### Priority
    ## lat/lon > y/x
    ## Direct Input > Reference File > Template File
    try:
        rlat
        rlon
    except:
        try:
            rlat = float(atr_ref['ref_lat'])
            rlon = float(atr_ref['ref_lon'])
        except:
            try: rlat,rlon = [float(i) for i in templateContents['pysar.seed.lalo'].split(',')]
            except: pass

    try:
        ry
        rx
    except:
        try:
            ry = int(atr_ref['ref_y'])
            rx = int(atr_ref['ref_x'])
        except:
            try: ry,rx       = [int(i)   for i in templateContents['pysar.seed.yx'].split(',')]
            except: pass

    ##### Check lalo / YX
    print '\n************** Reference Point ******************'
    try:
        rlat
        rlon
        y = sub.coord_geo2radar(rlat,atr,'lat')
        x = sub.coord_geo2radar(rlon,atr,'lon')
        0<= x <= width
        0<= y <= length
        rx = x
        ry = y
        print 'Reference point: lat = %.4f,   lon = %.4f'%(rlat,rlon)
        print '                 y   = %d,     x   = %d'%(ry,rx)
    except:
        print 'Skip input lat/lon reference point.'
        print 'Continue with the y/x reference point.'


    ######################### a. Read Mask File #########################
    ## Priority: Input mask file > pysar.mask.file 
    try:     maskFile
    except:
        try: maskFile = templateContents['pysar.mask.file']
        except:  print 'No mask found!';
    try:
        M,Matr = readfile.read(maskFile);
        print 'mask: '+maskFile
    except:
        print '---------------------------------------------------------'
        print 'WARNING: No mask, use the whole area as mask'
        print '---------------------------------------------------------'
        M = np.ones((length,width))

    ## Message
    try:
        rx
        ry
        0<= rx <= width
        0<= ry <= length
        if M[ry,rx] == 0:
            print 'Input point has 0 value in mask.'
    except: pass

    ######################### b. Stack ##################################
    stackFile = os.path.basename(File).split(ext)[0] + '_stack.h5'
    try:
        os.path.isfile(stackFile)
        stack,atrStack = readfile.read(stackFile)
        print 'read stack from file: '+stackFile
    except:
        stack = ut.stacking(File)
        atrStack = atr.copy()
        atrStack['FILE_TYPE'] = 'mask'
        writefile.write(stack,atrStack,stackFile)

    ## Message
    try:
        rx
        ry
        if stack[ry,rx] == 0:
            print 'Input point has nan value in data.'
    except: pass

    stack[M==0] = 0
    if np.nansum(M) == 0.0:
        print '\n*****************************************************'
        print   'ERROR:'
        print   'There is no pixel that has valid phase value in all datasets.' 
        print   'Check the file!'
        print   'Seeding failed'
        sys.exit(1)

    ######################### Check Method ##############################
    try:
        not stack[ry,rx] == 0
        method = 'input_coord'
    except:
        try:    method
        except: method = method_default
        print 'Skip input y/x reference point.'
        print 'Continue with '+method

    #h5file = h5py.File(File)

    ######################### Seeding ###################################
    ##### Sub-function
    def seed_method(method,File,stack,outFile,corFile=''):
        SeedingDone = 'no'
        next_method = method_default
        M = stack != 0

        if   method == 'manual':
            SeedingDone = seed_manual(File,stack,outFile)
            if SeedingDone == 'no':
                next_method = method_default
                print_warning(next_method)

        elif method == 'max_coherence':
            try:    SeedingDone = seed_max_coherence(File,M,outFile,corFile)
            except: SeedingDone = seed_max_coherence(File,M,outFile)
            if SeedingDone == 'no':
                next_method = 'random'
                print_warning(next_method)

        elif method == 'random':
            y,x = random_selection(stack)
            seed_xy(File,x,y,outFile)
            SeedingDone = 'yes'

        elif method == 'global_average':
            print '\n---------------------------------------------------------'
            print 'Automatically Seeding using Global Spatial Average Value '
            print '---------------------------------------------------------'
            print 'Calculating the global spatial average value for each epoch'+\
                  ' of all valid pixels ...'
            box = (0,0,width,length)
            meanList = ut.spatial_mean(File,M,box)
            seed_file(File,outFile,meanList,'','')
            SeedingDone = 'yes'

        return SeedingDone, next_method

    ##### Seeding
    SeedingDone = 'no'

    if method == 'input_coord':
        seed_xy(File,rx,ry,outFile)
        SeedingDone = 'yes'

    else:
        i = 0
        while SeedingDone == 'no' and i < 5:
            try:    SeedingDone,method = seed_method(method,File,stack,outFile,corFile)
            except: SeedingDone,method = seed_method(method,File,stack,outFile)
            i += 1
        if i >= 5:
            print 'ERROR: Seeding failed after more than '+str(i)+' times try ...'
            sys.exit(1)

################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])


