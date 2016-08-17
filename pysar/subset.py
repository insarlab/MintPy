#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Mar 2015: Add lat/lon option
# Yunjun, Jul 2015: Add 'coherence' option
# Yunjun, Sep 2015: Read .par file for jpg
#                   Add '.mli' and '.slc' option for Gamma product
#                   Make x/y/l/L option independent
#                   Add min/max check for input x/y/lat/lon
#                   Merge all files into PySAR, ROI_PAC, Image and GAMMA
#                   Merge 'interferograms','coherence','wrapped' into one
#                   Merge ROI_APC, Image and GAMMA into one
# Yunjun, Oct 2015: Support '.trans' file
# Yunjun, May 2016: Add -t and -- option, simplied code
#                   Add coord_geo2radar(),check_subset(),subset_attributes()
# Yunjun, Jun 2016: Add geo_box()
# Yunjun, Jul 2016: add parallel support
#                   add outlier fill option
# Yunjun, Aug 2016: add coord_geo2radar()


import os
import sys
import getopt

import h5py
import numpy as np

import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._pysar_utilities as ut


################################################################
def coord_geo2radar(geoCoord,atr,coordType):
    ## convert geo coordinates into radar coordinates
    ## for Geocoded file only
    ## Inputs:
    ##     geoCoord  : coordinate (list) in latitude/longitude in float
    ##     atr       : dictionary of file attributes
    ##     coordType : coordinate type: latitude, longitude
    ##
    ## Example:
    ##      300        = coord_radar2geo(32.104990,    atr,'lat')
    ##     [1000,1500] = coord_radar2geo([130.5,131.4],atr,'lon')

    try: atr['X_FIRST']
    except: print 'Support geocoded file only!'; sys.exit(1)

    ## Convert to List if input is String
    if isinstance(geoCoord,float):
        geoCoord = [geoCoord]

    radarCoord = []
    for i in range(len(geoCoord)):
        if   coordType.lower() in ['lat','latitude']:
            coord = int(np.floor((geoCoord[i]-float(atr['Y_FIRST']))/float(atr['Y_STEP'])))
        elif coordType.lower() in ['lon','longitude']:
            coord = int(np.floor((geoCoord[i]-float(atr['X_FIRST']))/float(atr['X_STEP'])))
        else: print 'Unrecognized coordinate type: '+coordType
        radarCoord.append(coord)
    radarCoord.sort()

    if len(radarCoord) == 1:
        radarCoord = radarCoord[0]

    return radarCoord


################################################################
def coord_radar2geo(radarCoord,atr,coordType):
    ## convert radar coordinates into geo coordinates
    ## for Geocoded file only
    ##
    ## Inputs:
    ##     radarCoord : coordinate (list) in row/col in int
    ##     atr        : dictionary of file attributes
    ##     coordType  : coordinate type: row, col, y, x
    ##
    ## Example:
    ##     32.104990     = coord_radar2geo(300,        atr,'y')
    ##     [130.5,131.4] = coord_radar2geo([1000,1500],atr,'x')

    try: atr['X_FIRST']
    except: print 'Support geocoded file only!'; sys.exit(1)
    
    ## Convert to List if input is String
    if isinstance(radarCoord,int):
        radarCoord = [radarCoord]

    geoCoord = []
    for i in range(len(radarCoord)):
        if   coordType.lower() in ['row','y']:
            coord = (radarCoord[i] + 0.5)*float(atr['Y_STEP']) + float(atr['Y_FIRST'])
        elif coordType.lower() in ['col','x','column']:
            coord = (radarCoord[i] + 0.5)*float(atr['X_STEP']) + float(atr['X_FIRST'])
        else: print 'Unrecognized coordinate type: '+coordType
        geoCoord.append(coord)
    geoCoord.sort()

    if len(geoCoord) == 1:
        geoCoord = geoCoord[0]

    return geoCoord

################################################################
def check_subset_range(sub_y,sub_x,atr):
    ## Check the subset range
    ## Inputs:
    ##     sub   : list of coordinates in row or column
    ##     atr   : dictionary of file attributes
    ##     type  : coordinate type: row, col
  
    width  = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])
    sub_y = sorted(sub_y)
    sub_x = sorted(sub_x)
  
    ## Check Y/Azimuth/Latitude subset range
    if not all(i>=0 and i<=length for i in sub_y) or not all(i>=0 and i<=width for i in sub_x):
        print 'WARNING: input index is out of data range!\nData range:'
        print     'range in y - 0:'+str(length)
        print     'range in x - 0:'+str(width)
        try:
            print 'range in latitude  - %.8f:%.8f'%(float(atr['Y_FIRST']),float(atr['Y_FIRST'])+float(atr['Y_STEP'])*length)
            print 'range in longitude - %.8f:%.8f'%(float(atr['X_FIRST']),float(atr['X_FIRST'])+float(atr['X_STEP'])*width)
        except: Geo=0
   
  
    if sub_y[0]<0:        sub_y[0]=0;      print 'WARNING: input y < min (0)! Set it to min.'
    if sub_y[1]>length:   sub_y[1]=length; print 'WARNING: input y > max ('+str(length)+')! Set it to max.'
    if sub_y[0]>length or sub_y[1]<0:
        print 'ERROR: input index is out of data range!'
        print     '              y - 0 : '+str(length)
        try:
            print 'range in geo: lat - %.8f:%.8f'%(float(atr['Y_FIRST']),float(atr['Y_FIRST'])+float(atr['Y_STEP'])*length)
        except: Geo=0
        sys.exit(1)
  
    ## Check X/Range/Longitude subset range
    if sub_x[0]<0:       sub_x[0]=0;      print 'WARNING: input x < min (0)! Set it to min.'
    if sub_x[1]>width:   sub_x[1]=width;  print 'WARNING: input x > max ('+str(width)+')! Set it to max x.'
    if sub_x[0]>width or sub_x[1]<0:
        print 'ERROR: input index is out of data range!'
        print     'range in rdr: x - 0 : '+str(width)
        try:
            print '              lon - %.8f:%.8f'%(float(atr['X_FIRST']),float(atr['X_FIRST'])+float(atr['X_STEP'])*width)
        except: Geo=0
        sys.exit(1)

    ##### Display subset range Info
    if not sub_y[1]-sub_y[0] == length:
        print 'subset in y direction - '+str(sub_y[0])+':'+str(sub_y[1])
        try:
            atr['Y_FIRST']
            sub_lat = [0]*2
            sub_lat[0] = float(atr['Y_FIRST']) + sub_y[0]*float(atr['Y_STEP'])
            sub_lat[1] = float(atr['Y_FIRST']) + sub_y[1]*float(atr['Y_STEP'])
            print 'subset in latitude  - %.8f:%.8f'%(sub_lat[0],sub_lat[1])
        except: pass
  
    if not sub_x[1]-sub_x[0] == width:
        print 'subset in x direction - '+str(sub_x[0])+':'+str(sub_x[1])
        try:
            atr['Y_FIRST']
            sub_lon = [0]*2
            sub_lon[0] = float(atr['X_FIRST']) + sub_x[0]*float(atr['X_STEP'])
            sub_lon[1] = float(atr['X_FIRST']) + sub_x[1]*float(atr['X_STEP'])
            print 'subset in longitude - %.8f:%.8f'%(sub_lon[0],sub_lon[1])
        except: pass
  
    return sub_y, sub_x

################################################################
def subset_attributes(atr_dict,sub_y,sub_x):
    #####
    atr = dict()
    for key, value in atr_dict.iteritems():  atr[key] = str(value)
  
    ##### Update attribute variable
    atr['FILE_LENGTH'] = str(sub_y[1]-sub_y[0])
    atr['WIDTH']       = str(sub_x[1]-sub_x[0])
    atr['YMAX']        = str(sub_y[1]-sub_y[0] - 1)
    atr['XMAX']        = str(sub_x[1]-sub_x[0] - 1)
  
    try:
        subset_y0_ori = int(atr['subset_y0'])
        atr['subset_y0'] = str(sub_y[0] + subset_y0_ori)
        atr['subset_y1'] = str(sub_y[1] + subset_y0_ori)
    except:
        atr['subset_y0'] = str(sub_y[0])
        atr['subset_y1'] = str(sub_y[1])
    try:
        subset_x0_ori = int(atr['subset_x0'])
        atr['subset_x0'] = str(sub_x[0] + subset_x0_ori)
        atr['subset_x1'] = str(sub_x[1] + subset_x0_ori)
    except:
        atr['subset_x0'] = str(sub_x[0])
        atr['subset_x1'] = str(sub_x[1])
    try:
        atr['Y_FIRST'] = str(float(atr['Y_FIRST'])+sub_y[0]*float(atr['Y_STEP']))
        atr['X_FIRST'] = str(float(atr['X_FIRST'])+sub_x[0]*float(atr['X_STEP']))
    except: pass
  
    try:
        atr['ref_y'] = str(int(atr['ref_y']) - sub_y[0])
        atr['ref_x'] = str(int(atr['ref_x']) - sub_x[0])
    except: pass
  
    return atr


###########################################################
def geo_box(atr):
    ## calculate coverage box in lalo
    ## box = (lon_min,lat_max,lon_max,lat_min)
    ## box = (UL_X,   UL_Y,   LR_X,   LR_Y)
  
    length = int(atr['FILE_LENGTH'])
    width  = int(atr['WIDTH'])
    lat_step = float(atr['Y_STEP'])
    lon_step = float(atr['X_STEP'])
    lat_max  = float(atr['Y_FIRST'])
    lon_min  = float(atr['X_FIRST'])
    lat_min  = lat_max + lat_step*length
    lon_max  = lon_min + lon_step*width
  
    box = (lon_min,lat_max,lon_max,lat_min)
  
    return box

def box_overlap_index(box1,box2):
    ## Calculate the overlap of two input boxes
    ##   and output the index box of the overlap in each's coord.
  
    x0 = max(box1[0],box2[0])
    y0 = max(box1[1],box2[1])
    x1 = min(box1[2],box2[2])
    y1 = min(box1[3],box2[3])
  
    if x0 >= x1 or y0 >= y1:
        print 'No overlap between two ranges!'
        print 'box 1:'
        print box1
        print 'box 2:'
        print box2
        sys.exit(1)
  
    box  = (x0,y0,x1,y1)
    idx1 = (box[0]-box1[0],box[1]-box1[1],box[2]-box1[0],box[3]-box1[1])
    idx2 = (box[0]-box2[0],box[1]-box2[1],box[2]-box2[0],box[3]-box2[1])
  
    return idx1, idx2


################################################################
def subset_file(File,sub_x,sub_y,outfill=np.nan,outName=''):

    ##### Overlap between subset and data range
    atr = readfile.read_attributes(File)
    width  = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])
    box1 = (0,0,width,length)
    box2 = (sub_x[0],sub_y[0],sub_x[1],sub_y[1])
    idx1,idx2 = box_overlap_index(box1,box2)
    print 'data   range:'
    print box1
    print 'subset range:'
    print box2
  
    ###########################  Data Read and Write  ######################
    k = atr['FILE_TYPE']
    print 'file type: '+k
    if outName == '':  outName = 'subset_'+os.path.basename(File)
  
    ##### Multiple Dataset File
    if k in ['timeseries','interferograms','wrapped','coherence']:
        ##### Input File Info
        h5file = h5py.File(File,'r')
        epochList = h5file[k].keys()
        epochList = sorted(epochList)
        print 'number of epochs: '+str(len(epochList))
  
        ##### Output File Info
        h5out = h5py.File(outName,'w')
        group = h5out.create_group(k)
        print 'writing >>> '+outName

    ## Loop
    if k == 'timeseries':
        for epoch in epochList:
            print epoch
            dset = h5file[k].get(epoch)
            data_overlap = dset[idx1[1]:idx1[3],idx1[0]:idx1[2]]
  
            data = np.ones((box2[3]-box2[1],box2[2]-box2[0]))*outfill
            data[idx2[1]:idx2[3],idx2[0]:idx2[2]] = data_overlap
  
            dset = group.create_dataset(epoch, data=data, compression='gzip')

        atr  = subset_attributes(atr,sub_y,sub_x)
        for key,value in atr.iteritems():   group.attrs[key] = value

    elif k in ['interferograms','wrapped','coherence']:
        for epoch in epochList:
            print epoch
            dset = h5file[k][epoch].get(epoch)
            atr  = h5file[k][epoch].attrs
            data_overlap = dset[idx1[1]:idx1[3],idx1[0]:idx1[2]]
  
            data = np.ones((box2[3]-box2[1],box2[2]-box2[0]))*outfill
            data[idx2[1]:idx2[3],idx2[0]:idx2[2]] = data_overlap
  
            atr  = subset_attributes(atr,sub_y,sub_x)
            gg = group.create_group(epoch)
            dset = gg.create_dataset(epoch, data=data, compression='gzip')
            for key, value in atr.iteritems():    gg.attrs[key] = value
  
    ##### Single Dataset File
    elif k == '.trans':
        rg_overlap,az_overlap,atr = readfile.read(File,idx1)
  
        rg = np.ones((box2[3]-box2[1],box2[2]-box2[0]))*outfill
        rg[idx2[1]:idx2[3],idx2[0]:idx2[2]] = rg_overlap
  
        az = np.ones((box2[3]-box2[1],box2[2]-box2[0]))*outfill
        az[idx2[1]:idx2[3],idx2[0]:idx2[2]] = az_overlap
  
        atr = subset_attributes(atr,sub_y,sub_x)
        writefile.write(rg,az,atr,outName)
    else:
        data_overlap,atr = readfile.read(File,idx1)
  
        data = np.ones((box2[3]-box2[1],box2[2]-box2[0]))*outfill
        data[idx2[1]:idx2[3],idx2[0]:idx2[2]] = data_overlap
  
        atr = subset_attributes(atr,sub_y,sub_x)
        writefile.write(data,atr,outName)
  
    ##### End Cleaning
    try:
        h5file.close()
        h5out.close()
    except: pass


################################################################
def Usage():
    print '''
****************************************************************
  Generate a subset of the dataset

  Usage:
      subset.py    file    template_file
      subset.py -f file -t template_file
      subset.py -f file -y subset_row      -x subset_column    -o output_name
      subset.py -f file -l subset_latitude -L subset_longitude -o output_name 

      -f : PySAR h5 files [interferograms, coherence, timeseries, velocity,
                          temporal_coherence, rmse, mask] 
           roi_pac files  [.unw, .cor, .dem, .hgt]
           image files    [jpeg, jpg, png, bmp]
           GAMMA files    [.mli, .slc]
      -o : output file name

      -t : template file with subset setting. Priority: subset range input > template
           pysar.subset.yx    = 300:800,1000:3500
           pysar.subset.lalo  = 30.2:30.5,130.1:131.3
      -x/--col : subset range in column
      -y/--row : subset range in row
      -l/--lat : subset range in latitude
      -L/--lon : subset range in longitude
      -r       : reference file, subset to the same lalo as reference file

      --parallel     : enable parallel computing
      --outfill-nan  : fill outside area with numpy.nan
      --outfill-zero : fill outside area with zero
      --outfill      : fill with input value if subset is outside of input data range
                       --outfill 0

  Example:
      subset.py velocity.h5 SinabungT495F50AlosA.template

      subset.py -f LoadedData.h5      -y 400:1500       -x 200:600
      subset.py -f geo_velocity.h5    -l 31.88:31.95    -L 130.85:130.92
      subset.py -f geo_timeseries.h5  --lat=30.5:30.8   --lon=130.3:130.9
      subset.py -f 030405_090801.unw  -t SinabungT495F50AlosA.template
      subset.py -f geo_incidence.h5   -r subset_geo_velocity.h5

      subset.py -f '*velocity*.h5,timeseries*.h5'  -y 400:1500  -x 200:600  --parallel
      subset.py -f geo_velocity.h5  -l 32.2:33.5  --outfill-nan
      subset.py -f Mask.h5          -x 500:3500   --outfill 0

****************************************************************
    '''

################################################################
def main(argv):

    global outName
    parallel = 'no'
  
    ############### Check Inputs ###############
    if len(sys.argv)>3:
        try:
            opts, args = getopt.getopt(argv,'f:l:L:o:t:x:y:r:',['lat=','lon=','row=','col=',\
                                            'parallel','outfill=','outfill-nan','outfill-zero'])
        except getopt.GetoptError:
            print 'Error while getting args'
            Usage() ; sys.exit(1)
  
        for opt,arg in opts:
            if   opt == '-f':   File         = arg.split(',')
            elif opt == '-o':   outName      = arg
            elif opt == '-t':   templateFile = arg
            elif opt == '-r':   refFile      = arg
            elif opt in ['-x','--col']  :   sub_x   = [int(i)   for i in arg.split(':')];    sub_x.sort()
            elif opt in ['-y','--row']  :   sub_y   = [int(i)   for i in arg.split(':')];    sub_y.sort()
            elif opt in ['-l','--lat']  :   sub_lat = [float(i) for i in arg.split(':')];  sub_lat.sort()
            elif opt in ['-L','--lon']  :   sub_lon = [float(i) for i in arg.split(':')];  sub_lon.sort()
            elif opt in '--parallel'    :   parallel = 'yes'
            elif opt in '--outfill'     :   out_fill = float(arg)
            elif opt in '--outfill-nan' :   out_fill = np.nan
            elif opt in '--outfill-zero':   out_fill = 0.0
  
    elif len(sys.argv)==3:
        File         = argv[0].split(',')
        templateFile = argv[1]
    elif len(sys.argv)==2:
        if argv[0] in ['-h','--help']:  Usage(); sys.exit()
        else: print '\nERROR: A minimum of 3 inputs is needed.\n'; Usage(); sys.exit()
    else: Usage(); sys.exit(1)
  
    ##### Check Input file Info
    print '\n**************** Subset *********************'
    fileList = ut.get_file_list(File)
    print 'number of file: '+str(len(fileList))
    print fileList
    atr = readfile.read_attributes(fileList[0])
  
    if len(fileList) == 1 and parallel == 'yes':
        print 'parallel is disabled for one input file.'
        parallel = 'no'

    ################## Subset Setting ###########
    try:
        atr['X_FIRST']
        print 'geo coordinate'
    except:
        print 'radar coordinate'
    ## Read Subset Inputs
    try:
        templateFile
        template = readfile.read_template(templateFile)
    except: pass
  
    try:
        refFile
        atr_ref = readfile.read_attributes(refFile)
        box_ref = geo_box(atr_ref)
        lat_ref = [box_ref[3],box_ref[1]]
        lon_ref = [box_ref[0],box_ref[2]]
    except: pass
  
    try:
        sub_lat
        sub_lon
    except:
        try:
            sub_lat = lat_ref
            sub_lon = lon_ref
        except:
            try:
                sub = template['pysar.subset.lalo'].split(',')
                sub_lat = [float(i) for i in sub[0].split(':')];  sub_lat.sort()
                sub_lon = [float(i) for i in sub[1].split(':')];  sub_lon.sort()
            except: pass; #print 'No pysar.subset.lalo option found in template file!'
  
    try:
        sub_y
        sub_x
    except:
        try:
            sub = template['pysar.subset.yx'].split(',')
            sub_y = [int(i) for i in sub[0].split(':')];  sub_y.sort()
            sub_x = [int(i) for i in sub[1].split(':')];  sub_x.sort()
        except: pass; #print 'No pysar.subset.yx option found in template file!'

    ## Check Subset Inputs Existed or not
    try:     sub_y
    except:
        try: sub_x
        except:
            try: sub_lat
            except:
                try: sub_lon
                except: print 'ERROR: no subset is setted.'; Usage(); sys.exit(1)

    ##### Subset range radar to geo
    width  = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])
    print 'input file length: '+str(length)
    print 'input file width : '+str(width)
  
    try: sub_y = coord_geo2radar(sub_lat,atr,'latitude')
    except:
        try:    sub_y
        except: sub_y = [0,length]
    try: sub_x = coord_geo2radar(sub_lon,atr,'longitude')
    except:
        try:    sub_x
        except: sub_x = [0,width]

    ##### Check subset range
    try:
        out_fill
    except:
        sub_y,sub_x = check_subset_range(sub_y,sub_x,atr)
        out_fill = np.nan
        if sub_y[1]-sub_y[0] == length and sub_x[1]-sub_x[0] == width:
            print 'Input subset range == data size, no need to subset.'
            sys.exit(0)

    ################### Subset #######################
    if parallel == 'no':
        for file in fileList:
            print '-------------------------------------------'
            print 'subseting  : '+file
            try:    subset_file(file,sub_x,sub_y,out_fill,outName)
            except: subset_file(file,sub_x,sub_y,out_fill)
  
    else:
        print '-------------------------'
        print 'parallel subseting ...'
        print '-------------------------'
        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        Parallel(n_jobs=num_cores)(delayed(subset_file)(file,sub_x,sub_y,out_fill) for file in fileList)
  
    print 'Done.'


###########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])



