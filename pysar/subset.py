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
#                   Add coord_geo2radar(),check_subset(),subset_attribute()
# Yunjun, Jun 2016: Add geo_box()
# Yunjun, Jul 2016: add parallel support
#                   add outlier fill option
# Yunjun, Aug 2016: add coord_geo2radar()
# Yunjun, Dec 2016: add cmdLineParse(), --tight option


import os
import sys
import argparse

import h5py
import numpy as np
import multiprocessing
from joblib import Parallel, delayed

import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._pysar_utilities as ut
from pysar._readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file


################################################################
def coord_geo2radar(geoCoord,atr,coordType):
    ## convert geo coordinates into radar coordinates (round to nearest integer)
    ## for Geocoded file only
    ## Inputs:
    ##     geoCoord  : coordinate (list) in latitude/longitude in float
    ##     atr       : dictionary of file attributes
    ##     coordType : coordinate type: latitude, longitude
    ##
    ## Example:
    ##      300        = coord_geo2radar(32.104990,    atr,'lat')
    ##     [1000,1500] = coord_geo2radar([130.5,131.4],atr,'lon')

    try: atr['X_FIRST']
    except: print 'Support geocoded file only!'; sys.exit(1)

    ## Convert to List if input is String
    if isinstance(geoCoord,float):
        geoCoord = [geoCoord]

    radarCoord = []
    coordType = coordType.lower()
    for i in range(len(geoCoord)):
        if   coordType in ['lat','latitude' ]:  coord = np.rint((geoCoord[i]-float(atr['Y_FIRST']))/float(atr['Y_STEP']))
        elif coordType in ['lon','longitude']:  coord = np.rint((geoCoord[i]-float(atr['X_FIRST']))/float(atr['X_STEP']))
        else: print 'Unrecognized coordinate type: '+coordType
        radarCoord.append(int(coord))
    #radarCoord.sort()

    if len(radarCoord) == 1:
        radarCoord = radarCoord[0]

    return radarCoord


################################################################
def coord_radar2geo(radarCoord,atr,coordType):
    ## convert radar coordinates into geo coordinates (pixel UL corner)
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
    coordType = coordType.lower()
    for i in range(len(radarCoord)):
        if   coordType in ['row','y']:           coord = radarCoord[i]*float(atr['Y_STEP']) + float(atr['Y_FIRST'])
        elif coordType in ['col','x','column']:  coord = radarCoord[i]*float(atr['X_STEP']) + float(atr['X_FIRST'])
        else: print 'Unrecognized coordinate type: '+coordType
        geoCoord.append(coord)
    #geoCoord.sort()

    if len(geoCoord) == 1:
        geoCoord = geoCoord[0]

    return geoCoord


################################################################
def check_box_within_data_coverage(pixel_box, atr_dict):
    '''Check the subset box's conflict with data coverage
    Inputs:
        pixel_box : 4-tuple of int, indicating y/x coordinates of subset
        atr       : dictionary of file attributes
    '''

    width  = int(atr_dict['WIDTH'])
    length = int(atr_dict['FILE_LENGTH'])
    sub_x = [pixel_box[0], pixel_box[2]]
    sub_y = [pixel_box[1], pixel_box[3]]

    if sub_y[0]>=length or sub_y[1]<=0 or sub_x[0]>=width or sub_x[1]<=0:
        print 'ERROR: input index is out of data range!'
        data_box = (0,0,width,length)
        print '\tdata   range in y/x: '+str(data_box)
        print '\tsubset range in y/x: '+str(pixel_box)
        print '\tdata   range in lat/lon: '+str(box_pixel2geo(data_box, atr_dict))
        print '\tsubset range in lat/lon: '+str(box_pixel2geo(pixel_box, atr_dict))
        sys.exit(1)

    ## Check Y/Azimuth/Latitude subset range
    if sub_y[0]<0:        sub_y[0]=0;      print 'WARNING: input y < min (0)! Set it to min.'
    if sub_y[1]>length:   sub_y[1]=length; print 'WARNING: input y > max ('+str(length)+')! Set it to max.'

    ## Check X/Range/Longitude subset range
    if sub_x[0]<0:       sub_x[0]=0;      print 'WARNING: input x < min (0)! Set it to min.'
    if sub_x[1]>width:   sub_x[1]=width;  print 'WARNING: input x > max ('+str(width)+')! Set it to max x.'


    out_box = (sub_x[0], sub_y[0], sub_x[1], sub_y[1])
    return out_box


################################################################
def subset_attribute(atr_dict, subset_box):
    '''Update attributes dictionary due to subset
    Inputs:
        atr_dict   : dict, data attributes to update
        subset_box : 4-tuple of int, subset box defined in (x0, y0, x1, y1)
    Outputs:
        atr      : dict, updated data attributes
    '''

    sub_x = [subset_box[0], subset_box[2]]
    sub_y = [subset_box[1], subset_box[3]]
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
def get_coverage_box(atr):
    '''Get Coverage Box of data in geo and pixel coordinates
    Inputs: atr - dict, meta data dictionary
    Outputs:
        pix_box : 4-tuple of int, defining in (UL_X, UL_Y, LR_X, LR_Y)
        geo_box : 4-tuple of float in lat/lon
    '''

    length = int(atr['FILE_LENGTH'])
    width  = int(atr['WIDTH'])

    # Get geo box
    try:
        lat_step = float(atr['Y_STEP'])
        lon_step = float(atr['X_STEP'])
        ul_lat  = float(atr['Y_FIRST'])
        ul_lon  = float(atr['X_FIRST'])
        lr_lat  = ul_lat + lat_step*length
        lr_lon  = ul_lon + lon_step*width
        geo_box = (ul_lon, ul_lat, lr_lon, lr_lat)
    except:
        geo_box = None

    # Get pixel box
    try:    pix_box = (int(atr['subset_x0']), int(atr['subset_y0']), int(atr['subset_x1']), int(atr['subset_y1']))
    except: pix_box = None

    return pix_box, geo_box


def read_subset_template2box(templateFile):
    '''Read pysar.subset.lalo/yx option from template file into box type
    Return None if not specified.
    '''
    tmpl = readfile.read_template(templateFile)
    try:
        sub = [i.strip() for i in tmpl['pysar.subset.lalo'].split(',')]
        sub_lat = sorted([float(i.strip()) for i in sub[0].split(':')])
        sub_lon = sorted([float(i.strip()) for i in sub[1].split(':')])
        geo_box = (sub_lon[0], sub_lat[1], sub_lon[1], sub_lat[0])
    except:
        geo_box = None
    try:
        sub = [i.strip() for i in tmpl['pysar.subset.yx'].split(',')]
        sub_y = sorted([int(i.strip()) for i in sub[0].split(':')])
        sub_x = sorted([int(i.strip()) for i in sub[1].split(':')])
        pix_box = (sub_x[0], sub_y[0], sub_x[1], sub_y[1])    
    except:
        pix_box = None
    return pix_box, geo_box


def bbox_geo2radar(geo_box, atr_rdr=dict(), transFile='geomap*.trans'):
    '''Calculate bounding box in x/y for file in radar coord, based on input geo box.
    Inputs:
        geo_box   - tuple of 4 float, indicating the UL/LR lon/lat 
        atr_rdr   - dict, attributes of file in radar coord
        transFile - string, path of transformation file, i.e. geomap_4rlks.trans
    Output:
        pix_box - tuple of 4 int, indicating the UL/LR x/y of the bounding box in radar coord
                  for the corresponding lat/lon coverage.
    '''
    lat = np.array([geo_box[3],geo_box[3],geo_box[1],geo_box[1]])
    lon = np.array([geo_box[0],geo_box[2],geo_box[0],geo_box[2]])
    y, x, y_res, x_res = ut.glob2radar(lat, lon, transFile, atr_rdr)
    buf = 10*(np.max([x_res, y_res]))
    pix_box = (np.min(x)-buf, np.min(y)-buf, np.max(x)+buf, np.max(y)+buf)
    return pix_box


def bbox_radar2geo(pix_box, atr_rdr=dict(), transFile='geomap*.trans'):
    '''Calculate bounding box in lat/lon for file in geo coord, based on input radar/pixel box
    Inputs:
        pix_box   - tuple of 4 int, indicating the UL/LR x/y
        atr_rdr   - dict, attributes of file in radar coord
        transFile - string, path of transformation file, i.e. geomap_4rlks.trans
    Output:
        geo_box - tuple of 4 float, indicating the UL/LR lon/lat of the bounding box
    '''
    x = np.array([pix_box[0],pix_box[2],pix_box[0],pix_box[2]])
    y = np.array([pix_box[1],pix_box[1],pix_box[3],pix_box[3]])
    lat, lon, lat_res, lon_res = ut.radar2glob(y, x, transFile, atr_rdr)
    buf = 10*(np.max([lat_res,lon_res]))
    geo_box = (np.min(lon)-buf, np.max(lat)+buf, np.max(lon)+buf, np.min(lat)-buf)
    return geo_box


def subset_box2inps(inps, pix_box, geo_box):
    '''Update inps.subset_y/x/lat/lon from pixel_box and geo_box'''
    if geo_box:
        inps.subset_lon = [geo_box[0], geo_box[2]]
        inps.subset_lat = [geo_box[1], geo_box[3]]
    else:
        inps.subset_lon = None
        inps.subset_lat = None
    if pix_box:
        inps.subset_x = [pix_box[0], pix_box[2]]
        inps.subset_y = [pix_box[1], pix_box[3]]
    else:
        inps.subset_x = None
        inps.subset_y = None
    return inps


def get_box_overlap_index(box1,box2):
    '''Get index box overlap area of two input boxes
    
    Inputs:
        box1/2 : 4-tuple of int, indicating coverage of box1/2
                 defining in (x0, y0, x1, y1)
    Outputs:
        overlap_idx_box1/2 : 4-tuple of int, indicating index of overlap area in box1/2
                             defining in (idx_x0, idx_y0, idx_x1, idx_y1)
    '''
    ## Calculate the overlap of two input boxes
    ##   and output the index box of the overlap in each's coord.

    # Calculate Overlap Box
    x0 = max(box1[0], box2[0])
    y0 = max(box1[1], box2[1])
    x1 = min(box1[2], box2[2])
    y1 = min(box1[3], box2[3])
    if x0 >= x1 or y0 >= y1:
        print 'ERROR: No overlap between two input box range!'
        print 'box 1: '+str(box1)
        print 'box 2: '+str(box2)
        sys.exit(1)
    overlap_box  = (x0,y0,x1,y1)

    # Overlap index for box1
    overlap_idx_box1 = (overlap_box[0]-box1[0], overlap_box[1]-box1[1], overlap_box[2]-box1[0], overlap_box[3]-box1[1])
    # Overlap index for box2
    overlap_idx_box2 = (overlap_box[0]-box2[0], overlap_box[1]-box2[1], overlap_box[2]-box2[0], overlap_box[3]-box2[1])

    return overlap_idx_box1, overlap_idx_box2


################################################################
def subset_input_dict2box(subset_dict, meta_dict):
    '''Convert subset inputs dict into box in radar and/or geo coord.
    Inputs:
        subset_dict : dict, including the following 4 objects:
                      subset_x   : list of 2 int,   subset in x direction,   default=None
                      subset_y   : list of 2 int,   subset in y direction,   default=None
                      subset_lat : list of 2 float, subset in lat direction, default=None
                      subset_lon : list of 2 float, subset in lon direction, default=None
        meta_dict   : dict, including the following items:
                      'WIDTH'      : int
                      'FILE_LENGTH': int
                      'X_FIRST'    : float, optional
                      'Y_FIRST'    : float, optional
                      'X_STEP'     : float, optional
                      'Y_STEP'     : float, optional
    Outputs:
        # box defined by 4-tuple of number, defining (left, upper, right, lower) coordinate,
        #                                            (UL_X, UL_Y,  LR_X,  LR_Y )
        pixel_box   : 4-tuple of int, in pixel unit - 1
        geo_box     : 4-tuple of float, in  lat/lon unit - degree
                      None if file is in radar coordinate.
    example:
        subset_dict = {'subset_x': None, 'subset_y': None, 'subset_lat': [30.5, 31.0], 'subset_lon': [130.0, 131.0]}
        subset_dict = {'subset_x': [100, 1100], 'subset_y': [2050, 2550], 'subset_lat': None, 'subset_lon': None}
        pixel_box          = subset_input_dict2box(subset_dict, pysar_meta_dict)[0]
        pixel_box, geo_box = subset_input_dict2box(subset_dict, pysar_meta_dict)
    '''

    # Data Coverage
    width = int(float(meta_dict['WIDTH']))
    length = int(float(meta_dict['FILE_LENGTH']))

    # Use subset_lat/lon input if existed,  priority: lat/lon > y/x > len/wid
    if subset_dict['subset_lat']:
        sub_y = coord_geo2radar(subset_dict['subset_lat'],meta_dict,'latitude')
    elif subset_dict['subset_y']:
        sub_y = subset_dict['subset_y']
    else:
        sub_y = [0,length]

    if subset_dict['subset_lon']:
        sub_x = coord_geo2radar(subset_dict['subset_lon'],meta_dict,'longitude')
    elif subset_dict['subset_x']:
        sub_x = subset_dict['subset_x']
    else:
        sub_x = [0,width]

    # Get subset box in y/x
    sub_x = sorted(sub_x)
    sub_y = sorted(sub_y)
    pixel_box = (sub_x[0],sub_y[0],sub_x[1],sub_y[1])

    # Get subset box in lat/lon from subset box in y/x
    geo_box = box_pixel2geo(pixel_box, meta_dict)

    return pixel_box, geo_box


def box_pixel2geo(pixel_box, meta_dict):
    '''Convert pixel_box to geo_box'''
    try:
        lon_step = float(meta_dict['X_STEP'])
        lat_step = float(meta_dict['Y_STEP'])
        ul_lon = float(meta_dict['X_FIRST']) + pixel_box[0]*lon_step
        ul_lat = float(meta_dict['Y_FIRST']) + pixel_box[1]*lat_step
        lr_lon = ul_lon + lon_step*(pixel_box[2]-pixel_box[0])
        lr_lat = ul_lat + lat_step*(pixel_box[3]-pixel_box[1])
        geo_box = (ul_lon, ul_lat, lr_lon, lr_lat)
    except:
        geo_box = None
    return geo_box


def box_geo2pixel(geo_box, meta_dict):
    '''Convert geo_box to pixel_box'''
    try:
        y = coord_geo2radar([geo_box[1],geo_box[3]], meta_dict, 'latitude')
        x = coord_geo2radar([geo_box[0],geo_box[2]], meta_dict, 'longitude')
        pixel_box = (x[0],y[0],x[1],y[1])
    except:
        pixel_box = None
    return pixel_box


################################################################
def subset_file(File, subset_dict, outFile=None):
    '''Subset file with
    Inputs:
        File        : str, path/name of file
        outFile     : str, path/name of output file
        subset_dict : dict, subsut parameter, including the following items:
                      subset_x   : list of 2 int,   subset in x direction,   default=None
                      subset_y   : list of 2 int,   subset in y direction,   default=None
                      subset_lat : list of 2 float, subset in lat direction, default=None
                      subset_lon : list of 2 float, subset in lon direction, default=None
                      fill_value : float, optional. filled value for area outside of data coverage. default=None
                                   None/not-existed to subset within data coverage only.
                      tight  : bool, tight subset or not, for lookup table file, i.e. geomap*.trans
    Outputs:
        outFile :  str, path/name of output file; 
                   outFile = 'subset_'+File, if File is in current directory;
                   outFile = File, if File is not in the current directory.
    '''
    
    # Input File Info
    try:  atr_dict = readfile.read_attribute(File)
    except:  return None
    width = int(atr_dict['WIDTH'])
    length = int(atr_dict['FILE_LENGTH'])
    k = atr_dict['FILE_TYPE']
    print 'subset '+k+' file: '+File+' ...'

    # Read Subset Inputs into 4-tuple box in pixel and geo coord
    pix_box, geo_box = subset_input_dict2box(subset_dict, atr_dict)

    # if fill_value exists and not None, subset data and fill assigned value for area out of its coverage.
    # otherwise, re-check subset to make sure it's within data coverage and initialize the matrix with np.nan
    outfill = False
    if 'fill_value' in subset_dict.keys() and subset_dict['fill_value']:
        outfill = True
    else:
        outfill = False
    if not outfill:
        pix_box = check_box_within_data_coverage(pix_box, atr_dict)
        subset_dict['fill_value'] = np.nan

    geo_box = box_pixel2geo(pix_box, atr_dict)
    data_box = (0,0,width,length)
    print 'data   range in y/x: '+str(data_box)
    print 'subset range in y/x: '+str(pix_box)
    print 'data   range in lat/lon: '+str(box_pixel2geo(data_box, atr_dict))
    print 'subset range in lat/lon: '+str(geo_box)

    if pix_box == data_box:
        print 'Subset range == data coverage, no need to subset. Skip.'
        return File

    # Calculate Subset/Overlap Index
    pix_box4data, pix_box4subset = get_box_overlap_index(data_box, pix_box)

    ###########################  Data Read and Write  ######################
    # Output File Name
    if not outFile:
        if os.getcwd() == os.path.dirname(os.path.abspath(File)):
            if 'tight' in subset_dict.keys() and subset_dict['tight']:
                outFile = os.path.splitext(File)[0]+'_tight'+os.path.splitext(File)[1]
            else:
                outFile = 'subset_'+os.path.basename(File)
        else:
            outFile = os.path.basename(File)
    print 'writing >>> '+outFile

    ##### Multiple Dataset File
    if k in ['timeseries','interferograms','wrapped','coherence']:
        ##### Open Input File 
        h5file = h5py.File(File,'r')
        epochList = sorted(h5file[k].keys())
        epochNum = len(epochList)
        if k in multi_dataset_hdf5_file:
            print 'number of acquisitions: '+str(epochNum)
        else:
            print 'number of interferograms: '+str(epochNum)

        ##### Open Output File
        h5out = h5py.File(outFile)
        group = h5out.create_group(k)

    ## Loop
    if k == 'timeseries':
        for i in range(epochNum):
            epoch = epochList[i]
            ut.print_progress(i+1, epochNum, prefix='', suffix=epoch)
            
            dset = h5file[k].get(epoch)
            data_overlap = dset[pix_box4data[1]:pix_box4data[3],pix_box4data[0]:pix_box4data[2]]

            data = np.ones((pix_box[3]-pix_box[1], pix_box[2]-pix_box[0]))*subset_dict['fill_value']
            data[pix_box4subset[1]:pix_box4subset[3], pix_box4subset[0]:pix_box4subset[2]] = data_overlap

            dset = group.create_dataset(epoch, data=data, compression='gzip')

        atr_dict  = subset_attribute(atr_dict, pix_box)
        for key,value in atr_dict.iteritems():   group.attrs[key] = value

    elif k in ['interferograms','wrapped','coherence']:
        for i in range(epochNum):
            epoch = epochList[i]
            ut.print_progress(i+1, epochNum, prefix='', suffix=epoch)
            
            dset = h5file[k][epoch].get(epoch)
            atr_dict  = h5file[k][epoch].attrs
            data_overlap = dset[pix_box4data[1]:pix_box4data[3],pix_box4data[0]:pix_box4data[2]]

            data = np.ones((pix_box[3]-pix_box[1], pix_box[2]-pix_box[0]))*subset_dict['fill_value']
            data[pix_box4subset[1]:pix_box4subset[3], pix_box4subset[0]:pix_box4subset[2]] = data_overlap

            atr_dict  = subset_attribute(atr_dict, pix_box)
            gg = group.create_group(epoch)
            dset = gg.create_dataset(epoch, data=data, compression='gzip')
            for key, value in atr_dict.iteritems():    gg.attrs[key] = value

    ##### Single Dataset File
    elif k in ['.jpeg','.jpg','.png','.ras','.bmp']:
        data, atr_dict = readfile.read(File, pix_box)
        atr_dict = subset_attribute(atr_dict, pix_box)
        writefile.write(data,atr_dict,outFile)

    elif k == '.trans':
        rg_overlap,az_overlap,atr_dict = readfile.read(File, pix_box4data)

        rg = np.ones((pix_box[3]-pix_box[1], pix_box[2]-pix_box[0]))*subset_dict['fill_value']
        rg[pix_box4subset[1]:pix_box4subset[3], pix_box4subset[0]:pix_box4subset[2]] = rg_overlap

        az = np.ones((pix_box[3]-pix_box[1], pix_box[2]-pix_box[0]))*subset_dict['fill_value']
        az[pix_box4subset[1]:pix_box4subset[3], pix_box4subset[0]:pix_box4subset[2]] = az_overlap

        atr_dict = subset_attribute(atr_dict, pix_box)
        writefile.write(rg,az,atr_dict,outFile)
    else:
        data_overlap,atr_dict = readfile.read(File, pix_box4data)

        data = np.ones((pix_box[3]-pix_box[1], pix_box[2]-pix_box[0]))*subset_dict['fill_value']
        data[pix_box4subset[1]:pix_box4subset[3], pix_box4subset[0]:pix_box4subset[2]] = data_overlap

        atr_dict = subset_attribute(atr_dict, pix_box)
        writefile.write(data, atr_dict, outFile)

    ##### End Cleaning
    try:
        h5file.close()
        h5out.close()
    except: pass
    
    return outFile


def subset_file_list(fileList, inps):
    '''Subset file list'''
    # check outfile and parallel option
    if len(fileList) > 1:
        inps.outfile = None
    elif len(fileList) == 1 and inps.parallel:
        inps.parallel =  False
        print 'parallel processing is diabled for one input file'

    ##### Subset files
    if inps.parallel:
        num_cores = min(multiprocessing.cpu_count(), len(fileList))
        print 'parallel processing using %d cores ...'%(num_cores)
        Parallel(n_jobs=num_cores)(delayed(subset_file)(file, vars(inps)) for file in fileList)
    else:
        for File in fileList:
            print '----------------------------------------------------'
            subset_file(File, vars(inps), inps.outfile)
    return


###########################################################################################
EXAMPLE='''example:
  subset.py unwrapIfgram.h5    -y    400  1500   -x    200   600
  subset.py geo_velocity.h5    -l    30.5 30.8   -L    130.3 130.9
  subset.py geo_timeseries.h5  --lat 30.5 30.8   --lon 130.3 130.9
  subset.py 030405_090801.unw  -t SinabungT495F50AlosA.template
  subset.py geo_incidence.h5   -r subset_geo_velocity.h
  subset.py *velocity*.h5 timeseries*.h5  -y 400 1500  -x 200 600
  subset.py geo_velocity.h5    -l 32.2:33.5  --outfill-nan
  subset.py Mask.h5            -x 500:3500   --outfill 0
  subset.py geomap_4rlks.trans --tight
  
  subset.py unwrapIfgram.h5 coherence.h5 geomap*.trans  -l 33.10 33.50 -L 131.30 131.80 --bbox geomap_4rlks.trans
  subset.py *.unw *.cor *.trans *.dem  -y 50 450 -x 1300 1800 --bbox geomap_4rlks.trans
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Generate a subset from file/dataset',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)
    parser.add_argument('file', nargs='+', help='File(s) to subset/crop')

    parser.add_argument('-x', dest='subset_x', type=int, nargs=2, help='subset range in x/cross-track/column direction')
    parser.add_argument('-y', dest='subset_y', type=int, nargs=2, help='subset range in y/along-track/row direction')
    parser.add_argument('-l','--lat', dest='subset_lat', type=float, nargs=2, help='subset range in latitude')
    parser.add_argument('-L','--lon', dest='subset_lon', type=float, nargs=2, help='subset range in column\n\n')

    parser.add_argument('-t','--template', dest='template_file',\
                        help='template file with subset setting.  i.e. \n'
                             'pysar.subset.yx    = 300:800,1000:3500\n'
                             'pysar.subset.lalo  = 30.2:30.5,130.1:131.3')
    parser.add_argument('-r','--reference',\
                        help='reference file, subset to the same lalo as reference file')
    parser.add_argument('--tight', action='store_true',\
                        help='subset geomap_*.trans file based on non-zero values.\n'+\
                             'For geocoded file(s) only'
                             'A convenient way to get rid of extra wide space due to "too large" DEM.\n\n')

    parser.add_argument('--outfill', dest='fill_value', type=float,\
                        help="fill subset area out of data coverage with input value. i.e. \n"
                             "np.nan, 0, 1000, ... \n"
                             "By default, it's None for no-outfill.")
    parser.add_argument('--no-parallel',dest='parallel',action='store_false',default=True,\
                        help='Disable parallel processing. Diabled auto for 1 input file.\n\n')

    parser.add_argument('-o','--output', dest='outfile',\
                        help='output file name\n'+\
                             'add prefix "subset_" if input/output files are in the same directory;\n'+\
                             'same filename otherwise.')
    
    dset_group = parser.add_argument_group('Datasets',\
                                           'Create a subset of entire dataset in radar using y/x or lat/lon option\n'+\
                                           'Including *.trans and *.dem in geo coord.')
    dset_group.add_argument('--bbox', dest='trans_file',\
                            help='calculate bounding box in geo/radar coord from input radar/geo subset range\n'+\
                                 'using transformation file, i.e. geomap_4rlks.trans\n'+\
                                 'All input radar coord file should be same size/coverage; same for all geo coord files.')

    inps = parser.parse_args()
    return inps


###########################################################################################
def main(argv):
    inps = cmdLineParse()
    inps.file = ut.get_file_list(inps.file)

    #print '\n**************** Subset *********************'
    atr = readfile.read_attribute(inps.file[0])

    ##### Convert All Inputs into subset_y/x/lat/lon
    # Input Priority: subset_y/x/lat/lon > reference > template > tight
    if not inps.subset_x and not inps.subset_y and not inps.subset_lat and not inps.subset_lon:
        # 1. Read subset info from Reference File
        if inps.reference:
            ref_atr = readfile.read_attribute(inps.reference)
            pix_box, geo_box = get_coverage_box(ref_atr)
            print 'using subset info from '+inps.reference

        # 2. Read subset info from template options
        elif inps.template_file:
            pix_box, geo_box = read_subset_template2box(inps.template_file)
            print 'using subset info from '+inps.template_file

        # 3. Use subset from tight info
        elif inps.tight:
            if atr['FILE_TYPE']=='.trans':
                # Non-zero area in geomap_*.trans file, accurate
                trans_rg, trans_atr = readfile.read(inps.file[0], (), 'range')
                idx_row, idx_col = np.nonzero(trans_rg)
                pix_box = (np.min(idx_col)-10, np.min(idx_row)-10, np.max(idx_col)+10, np.max(idx_row)+10)
                geo_box = box_pixel2geo(pix_box, trans_atr)
            else:
                print 'ERROR: --tight option only works for geomap_*.trans file.\n'
                inps.tight = False
                sys.exit(1)

            ## from LAT/LON_REF*, which is not accurate
            #lats = [atr['LAT_REF1'], atr['LAT_REF3'], atr['LAT_REF4'], atr['LAT_REF2']]
            #lons = [atr['LON_REF1'], atr['LON_REF3'], atr['LON_REF4'], atr['LON_REF2']]
            #lats = [float(i) for i in lats]
            #lons = [float(i) for i in lons]
            #lalo_buff = min([max(lats)-min(lats), max(lons)-min(lons)]) * 0.05
            #geo_box = (min(lons)-lalo_buff, max(lats)+lalo_buff, max(lons)+lalo_buff, min(lats)-lalo_buff)
            #pix_box = None
            #if not inps.fill_value: inps.fill_value = np.nan
            #print 'using subset info from scene footprint - LAT/LON_REF1/2/3/4'
        else:
            raise Exception('No subset inputs found!')
        # Update subset_y/x/lat/lon
        inps = subset_box2inps(inps, pix_box, geo_box)

    ##### --bbox option
    if inps.trans_file:
        ## Seperate files in radar and geo coord
        rdrFileList = []
        geoFileList = []
        for File in inps.file:
            atr = readfile.read_attribute(File)
            if 'X_FIRST' in atr.keys():
                geoFileList.append(File)
            else:
                rdrFileList.append(File)
        
        ## Calculate bbox
        rdrFile = rdrFileList[0]
        atr_rdr = readfile.read_attribute(rdrFile)
        if inps.subset_lat and inps.subset_lon:
            print 'use subset input in lat/lon'
            print 'calculate corresponding bounding box in radar coordinate.'
            geo_box = (inps.subset_lon[0], inps.subset_lat[1], inps.subset_lon[1], inps.subset_lat[0])
            pix_box = bbox_geo2radar(geo_box, atr_rdr, inps.trans_file)
        else:
            print 'use subset input in y/x'
            print 'calculate corresponding bounding box in geo coordinate.'
            pix_box = (inps.subset_x[0], inps.subset_y[0], inps.subset_x[1], inps.subset_y[1])
            geo_box = bbox_radar2geo(pix_box, atr_rdr, inps.trans_file)
        print 'geo   box: '+str(geo_box)
        print 'pixel box: '+str(pix_box)
        
        ## Subset files
        inps.fill_value = 0
        print '--------------------------------------------'
        print 'subseting dataset in geo coord geo_box: '+str(geo_box)
        inps = subset_box2inps(inps, None, geo_box)
        subset_file_list(geoFileList, inps)
        print '--------------------------------------------'
        print 'subseting dataset in radar coord pix_box: '+str(pix_box)
        inps = subset_box2inps(inps, pix_box, None)
        subset_file_list(rdrFileList, inps)

    else:
        ##### Subset files
        subset_file_list(inps.file, inps)

    print 'Done.'
    return


###########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])



