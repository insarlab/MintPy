#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013, Zhang Yunjun, Heresh Fattahi          #
# Author:  Zhang Yunjun, Heresh Fattahi                    #
############################################################


import os
import sys
import argparse
import numpy as np
from pysar.utils import datetime as ptime, readfile, writefile, utils as ut


###########################################################################################
TEMPLATE = """template
## if both yx and lalo are specified, use lalo option unless a) no lookup file AND b) dataset is in radar coord
pysar.subset.yx       = auto    #[1800:2000,700:800 / no], auto for no
pysar.subset.lalo     = auto    #[31.5:32.5,130.5:131.0 / no], auto for no
"""
EXAMPLE = """example:
  subset.py INPUTS/ifgramStack.h5 -y 400  1500 -x 200   600
  subset.py geo_velocity.h5       -l 30.5 30.8 -L 130.3 130.9
  subset.py 030405_090801.unw     -t SinabungT495F50AlosA.template

  # subset to the same coverage as the reference file
  subset.py geo_incidence.h5 -r subset_geo_velocity.h5

  # multiple files input
  subset.py *velocity*.h5 timeseries*.h5  -y 400 1500  -x 200 600

  # crop to larger area with custom fill value 
  subset.py geo_velocity.h5 -l 32.2 33.5  --outfill-nan
  subset.py Mask.h5 -x 500 3500 --outfill 0

  # "tight" subset for geocoded lookup table larger than data file
  subset.py geomap_4rlks.trans --tight
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Generate a subset from file/dataset',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=TEMPLATE+'\n'+EXAMPLE)
    parser.add_argument('file', nargs='+', help='File(s) to subset/crop')

    parser.add_argument('-x', dest='subset_x', type=int, nargs=2,
                        help='subset range in x/cross-track/column direction')
    parser.add_argument('-y', dest='subset_y', type=int, nargs=2,
                        help='subset range in y/along-track/row direction')
    parser.add_argument('-l', '--lat', dest='subset_lat',
                        type=float, nargs=2, help='subset range in latitude')
    parser.add_argument('-L', '--lon', dest='subset_lon',
                        type=float, nargs=2, help='subset range in column\n\n')

    parser.add_argument('-t', '--template', dest='template_file',
                        help='template file with subset setting.  i.e. \n'
                             'pysar.subset.yx    = 300:800,1000:3500\n'
                             'pysar.subset.lalo  = 30.2:30.5,130.1:131.3')
    parser.add_argument('-r', '--reference',
                        help='reference file, subset to the same lalo as reference file')
    parser.add_argument('--tight', action='store_true',
                        help='subset geomap_*.trans file based on non-zero values.\n' +
                             'For geocoded file(s) only'
                             'A convenient way to get rid of extra wide space due to "too large" DEM.\n\n')

    parser.add_argument('--outfill', dest='fill_value', type=float,
                        help="fill subset area out of data coverage with input value. i.e. \n"
                             "np.nan, 0, 1000, ... \n"
                             "By default, it's None for no-outfill.")
    parser.add_argument('--no-parallel', dest='parallel', action='store_false', default=True,
                        help='Disable parallel processing. Diabled auto for 1 input file.\n\n')

    parser.add_argument('-o', '--output', dest='outfile',
                        help='output file name\n' +
                             'add prefix "subset_" if input/output files are in the same directory;\n' +
                             'same filename otherwise.')

    dset_group = parser.add_argument_group('Datasets',
                                           'Create a subset of entire dataset in radar using y/x or lat/lon option\n' +
                                           'Including *.trans and *.dem in geo coord.')
    dset_group.add_argument('--lookup', dest='lookup_file',
                            help='calculate bounding box in geo/radar coord from input radar/geo subset range\n' +
                                 'using transformation file, i.e. geomap_4rlks.trans\n' +
                                 'All input radar coord file should be same size/coverage; same for all geo coord files.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


################################################################
def check_box_within_data_coverage(pixel_box, metadata):
    """Check the subset box's conflict with data coverage
    Inputs:
        pixel_box : 4-tuple of int, indicating y/x coordinates of subset
        atr       : dictionary of file attributes
    """

    width = int(metadata['WIDTH'])
    length = int(metadata['LENGTH'])
    sub_x = [pixel_box[0], pixel_box[2]]
    sub_y = [pixel_box[1], pixel_box[3]]

    if sub_y[0] >= length or sub_y[1] <= 0 or sub_x[0] >= width or sub_x[1] <= 0:
        print('ERROR: input index is out of data range!')
        data_box = (0, 0, width, length)
        print('\tdata   range in x/y: '+str(data_box))
        print('\tsubset range in x/y: '+str(pixel_box))
        print('\tdata   range in lat/lon: '+str(box_pixel2geo(data_box, metadata)))
        print('\tsubset range in lat/lon: '+str(box_pixel2geo(pixel_box, metadata)))
        sys.exit(1)

    # Check Y/Azimuth/Latitude subset range
    if sub_y[0] < 0:
        sub_y[0] = 0
        print('WARNING: input y < min (0)! Set it to min.')
    if sub_y[1] > length:
        sub_y[1] = length
        print('WARNING: input y > max ('+str(length)+')! Set it to max.')

    # Check X/Range/Longitude subset range
    if sub_x[0] < 0:
        sub_x[0] = 0
        print('WARNING: input x < min (0)! Set it to min.')
    if sub_x[1] > width:
        sub_x[1] = width
        print('WARNING: input x > max ('+str(width)+')! Set it to max.')

    out_box = (sub_x[0], sub_y[0], sub_x[1], sub_y[1])
    return out_box


###########################################################
def get_coverage_box(atr):
    """Get Coverage Box of data in geo and pixel coordinates
    Inputs: atr - dict, meta data dictionary
    Outputs:
        pix_box : 4-tuple of int, defining in (UL_X, UL_Y, LR_X, LR_Y)
        geo_box : 4-tuple of float in lat/lon
    """

    length = int(atr['LENGTH'])
    width = int(atr['WIDTH'])

    # Get geo box
    try:
        lat_step = float(atr['Y_STEP'])
        lon_step = float(atr['X_STEP'])
        ul_lat = float(atr['Y_FIRST'])
        ul_lon = float(atr['X_FIRST'])
        lr_lat = ul_lat + lat_step*length
        lr_lon = ul_lon + lon_step*width
        geo_box = (ul_lon, ul_lat, lr_lon, lr_lat)
    except:
        geo_box = None

    # Get pixel box
    try:
        pix_box = (int(atr['SUBSET_XMIN']),
                   int(atr['SUBSET_YMIN']),
                   int(atr['SUBSET_XMAX']),
                   int(atr['SUBSET_YMAX']))
    except:
        pix_box = None

    return pix_box, geo_box


def read_subset_template2box(template_file):
    """Read pysar.subset.lalo/yx option from template file into box type
    Return None if not specified.
    """
    tmpl = readfile.read_template(template_file)
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


def bbox_geo2radar(geo_box, atr_rdr=dict(), lookup_file=None, print_msg=False):
    """Calculate bounding box in x/y for file in radar coord, based on input geo box.
    Inputs:
        geo_box    - tuple of 4 float, indicating the UL/LR lon/lat 
        atr_rdr    - dict, attributes of file in radar coord
        lookup_file - string / list of string, path of transformation file, i.e. geomap_4rlks.trans
    Output:
        pix_box - tuple of 4 int, indicating the UL/LR x/y of the bounding box in radar coord
                  for the corresponding lat/lon coverage.
    """
    lat = np.array([geo_box[3], geo_box[3], geo_box[1], geo_box[1]])
    lon = np.array([geo_box[0], geo_box[2], geo_box[0], geo_box[2]])
    if 'Y_FIRST' in atr_rdr.keys():
        y = ut.coord_geo2radar(lat, atr_rdr, 'lat')
        x = ut.coord_geo2radar(lon, atr_rdr, 'lon')
        pix_box = (x[0], y[2], x[1], y[0])
    else:
        y, x, y_res, x_res = ut.glob2radar(lat, lon,
                                           lookup_file,
                                           atr_rdr,
                                           print_msg=print_msg)
        buf = 2*(np.max(np.abs([x_res, y_res])))
        pix_box = (np.min(x)-buf, np.min(y)-buf,
                   np.max(x)+buf, np.max(y)+buf)
    return pix_box


def bbox_radar2geo(pix_box, atr_rdr=dict(), lookup_file=None, print_msg=False):
    """Calculate bounding box in lat/lon for file in geo coord, based on input radar/pixel box
    Inputs:
        pix_box    - tuple of 4 int, indicating the UL/LR x/y
        atr_rdr    - dict, attributes of file in radar coord
        lookup_file - string / list of string, path of transformation file, i.e. geomap_4rlks.trans
    Output:
        geo_box - tuple of 4 float, indicating the UL/LR lon/lat of the bounding box
    """
    x = np.array([pix_box[0], pix_box[2], pix_box[0], pix_box[2]])
    y = np.array([pix_box[1], pix_box[1], pix_box[3], pix_box[3]])
    if 'Y_FIRST' in atr_rdr.keys():
        lat = ut.coord_radar2geo(y, atr_rdr, 'y')
        lon = ut.coord_radar2geo(x, atr_rdr, 'x')
        geo_box = (lon[0], lat[0], lon[1], lat[2])
    else:
        lat, lon, lat_res, lon_res = ut.radar2glob(y, x,
                                                   lookup_file,
                                                   atr_rdr,
                                                   print_msg=print_msg)
        buf = 2*(np.max(np.abs([lat_res, lon_res])))
        geo_box = (np.min(lon)-buf, np.max(lat)+buf,
                   np.max(lon)+buf, np.min(lat)-buf)
    return geo_box


def subset_box2inps(inps, pix_box, geo_box):
    """Update inps.subset_y/x/lat/lon from pixel_box and geo_box"""
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


def get_box_overlap_index(box1, box2):
    """Get index box overlap area of two input boxes

    Inputs:
        box1/2 : 4-tuple of int, indicating coverage of box1/2
                 defining in (x0, y0, x1, y1)
    Outputs:
        overlap_idx_box1/2 : 4-tuple of int, indicating index of overlap area in box1/2
                             defining in (idx_x0, idx_y0, idx_x1, idx_y1)
    """
    # Calculate the overlap of two input boxes
    # and output the index box of the overlap in each's coord.

    # Calculate Overlap Box
    x0 = max(box1[0], box2[0])
    y0 = max(box1[1], box2[1])
    x1 = min(box1[2], box2[2])
    y1 = min(box1[3], box2[3])
    if x0 >= x1 or y0 >= y1:
        print('ERROR: No overlap between two input box range!')
        print('box 1: '+str(box1))
        print('box 2: '+str(box2))
        sys.exit(1)
    overlap_box = (x0, y0, x1, y1)

    # Overlap index for box1
    overlap_idx_box1 = (overlap_box[0] - box1[0],
                        overlap_box[1] - box1[1],
                        overlap_box[2] - box1[0],
                        overlap_box[3] - box1[1])
    # Overlap index for box2
    overlap_idx_box2 = (overlap_box[0] - box2[0],
                        overlap_box[1] - box2[1],
                        overlap_box[2] - box2[0],
                        overlap_box[3] - box2[1])

    return overlap_idx_box1, overlap_idx_box2


################################################################
def subset_input_dict2box(subset_dict, meta_dict):
    """Convert subset inputs dict into box in radar and/or geo coord.
    Inputs:
        subset_dict : dict, including the following 4 objects:
                      subset_x   : list of 2 int,   subset in x direction,   default=None
                      subset_y   : list of 2 int,   subset in y direction,   default=None
                      subset_lat : list of 2 float, subset in lat direction, default=None
                      subset_lon : list of 2 float, subset in lon direction, default=None
        meta_dict   : dict, including the following items:
                      'WIDTH'      : int
                      'LENGTH': int
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
    """

    # Data Coverage
    width = int(float(meta_dict['WIDTH']))
    length = int(float(meta_dict['LENGTH']))

    # Use subset_lat/lon input if existed,  priority: lat/lon > y/x > len/wid
    if subset_dict['subset_lat']:
        sub_y = ut.coord_geo2radar(subset_dict['subset_lat'], meta_dict, 'latitude')
    elif subset_dict['subset_y']:
        sub_y = subset_dict['subset_y']
    else:
        sub_y = [0, length]

    if subset_dict['subset_lon']:
        sub_x = ut.coord_geo2radar(subset_dict['subset_lon'], meta_dict, 'longitude')
    elif subset_dict['subset_x']:
        sub_x = subset_dict['subset_x']
    else:
        sub_x = [0, width]

    # Get subset box in y/x
    sub_x = sorted(sub_x)
    sub_y = sorted(sub_y)
    pixel_box = (sub_x[0], sub_y[0], sub_x[1], sub_y[1])

    # Get subset box in lat/lon from subset box in y/x
    geo_box = box_pixel2geo(pixel_box, meta_dict)

    return pixel_box, geo_box


def box_pixel2geo(pixel_box, meta_dict):
    """Convert pixel_box to geo_box"""
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
    """Convert geo_box to pixel_box"""
    try:
        y = ut.coord_geo2radar([geo_box[1], geo_box[3]], meta_dict, 'latitude')
        x = ut.coord_geo2radar([geo_box[0], geo_box[2]], meta_dict, 'longitude')
        pixel_box = (x[0], y[0], x[1], y[1])
    except:
        pixel_box = None
    return pixel_box


################################################################
def subset_file(fname, subset_dict_input, out_file=None):
    """Subset file with
    Inputs:
        fname        : str, path/name of file
        out_file     : str, path/name of output file
        subset_dict : dict, subsut parameter, including the following items:
                      subset_x   : list of 2 int,   subset in x direction,   default=None
                      subset_y   : list of 2 int,   subset in y direction,   default=None
                      subset_lat : list of 2 float, subset in lat direction, default=None
                      subset_lon : list of 2 float, subset in lon direction, default=None
                      fill_value : float, optional. filled value for area outside of data coverage. default=None
                                   None/not-existed to subset within data coverage only.
                      tight  : bool, tight subset or not, for lookup table file, i.e. geomap*.trans
    Outputs:
        out_file :  str, path/name of output file; 
                   out_file = 'subset_'+fname, if fname is in current directory;
                   out_file = fname, if fname is not in the current directory.
    """

    # Input File Info
    try:
        atr = readfile.read_attribute(fname)
    except:
        return None

    width = int(atr['WIDTH'])
    length = int(atr['LENGTH'])
    k = atr['FILE_TYPE']
    print('subset '+k+' file: '+fname+' ...')

    subset_dict = subset_dict_input.copy()
    # Read Subset Inputs into 4-tuple box in pixel and geo coord
    pix_box, geo_box = subset_input_dict2box(subset_dict, atr)

    # if fill_value exists and not None, subset data and fill assigned value for area out of its coverage.
    # otherwise, re-check subset to make sure it's within data coverage and initialize the matrix with np.nan
    outfill = False
    if 'fill_value' in subset_dict.keys() and subset_dict['fill_value']:
        outfill = True
    else:
        outfill = False
    if not outfill:
        pix_box = check_box_within_data_coverage(pix_box, atr)
        subset_dict['fill_value'] = np.nan

    geo_box = box_pixel2geo(pix_box, atr)
    data_box = (0, 0, width, length)
    print('data   range in y/x: '+str(data_box))
    print('subset range in y/x: '+str(pix_box))
    print('data   range in lat/lon: '+str(box_pixel2geo(data_box, atr)))
    print('subset range in lat/lon: '+str(geo_box))

    if pix_box == data_box:
        print('Subset range == data coverage, no need to subset. Skip.')
        return fname

    # Calculate Subset/Overlap Index
    pix_box4data, pix_box4subset = get_box_overlap_index(data_box, pix_box)

    ###########################  Data Read and Write  ######################
    # Output File Name
    if not out_file:
        if os.getcwd() == os.path.dirname(os.path.abspath(fname)):
            if 'tight' in subset_dict.keys() and subset_dict['tight']:
                out_file = '{}_tight{}'.format(os.path.splitext(fname)[0],
                                               os.path.splitext(fname)[1])
            else:
                out_file = 'subset_'+os.path.basename(fname)
        else:
            out_file = os.path.basename(fname)
    print('writing >>> '+out_file)

    # subset datasets one by one
    dsNames = readfile.get_dataset_list(fname)
    maxDigit = max([len(i) for i in dsNames])
    dsDict = dict()
    for dsName in dsNames:
        print('subsetting {d:<{w}} from {f} ...'.format(
            d=dsName, w=maxDigit, f=os.path.basename(fname)))
        data = readfile.read(fname, datasetName=dsName, print_msg=False)[0]

        # subset 2D data
        if len(data.shape) == 2:
            data_overlap = data[pix_box4data[1]:pix_box4data[3],
                                pix_box4data[0]:pix_box4data[2]]
            data = np.ones((pix_box[3] - pix_box[1],
                            pix_box[2] - pix_box[0])) * subset_dict['fill_value']
            data[pix_box4subset[1]:pix_box4subset[3],
                 pix_box4subset[0]:pix_box4subset[2]] = data_overlap

        # subset 3D data
        elif len(data.shape) == 3:
            data_overlap = data[:,
                                pix_box4data[1]:pix_box4data[3],
                                pix_box4data[0]:pix_box4data[2]]
            data = np.ones((data.shape[0],
                            pix_box[3] - pix_box[1],
                            pix_box[2] - pix_box[0])) * subset_dict['fill_value']
            data[:,
                 pix_box4subset[1]:pix_box4subset[3],
                 pix_box4subset[0]:pix_box4subset[2]] = data_overlap

        dsDict[dsName] = data

    atr = ut.subset_attribute(atr, pix_box)
    writefile.write(dsDict, out_file=out_file, metadata=atr, ref_file=fname)
    return out_file


def subset_file_list(fileList, inps):
    """Subset file list"""
    # check outfile and parallel option
    if inps.parallel:
        num_cores, inps.parallel, Parallel, delayed = ut.check_parallel(len(fileList))

    # Subset files
    if len(fileList) == 1:
        subset_file(fileList[0], vars(inps), inps.outfile)

    elif inps.parallel:
        #num_cores = min(multiprocessing.cpu_count(), len(fileList))
        # print 'parallel processing using %d cores ...'%(num_cores)
        Parallel(n_jobs=num_cores)(delayed(subset_file)(file, vars(inps)) for file in fileList)
    else:
        for fname in fileList:
            print('----------------------------------------------------')
            subset_file(fname, vars(inps))
    return


def read_aux_subset2inps(inps):
    # Convert All Inputs into subset_y/x/lat/lon
    # Input Priority: subset_y/x/lat/lon > reference > template > tight
    if all(not i for i in [inps.subset_x,
                           inps.subset_y,
                           inps.subset_lat,
                           inps.subset_lon]):
        # 1. Read subset info from Reference File
        if inps.reference:
            ref_atr = readfile.read_attribute(inps.reference)
            pix_box, geo_box = get_coverage_box(ref_atr)
            print('using subset info from '+inps.reference)

        # 2. Read subset info from template options
        elif inps.template_file:
            pix_box, geo_box = read_subset_template2box(inps.template_file)
            print('using subset info from '+inps.template_file)

        # 3. Use subset from tight info
        elif inps.tight:
            inps.lookup_file = ut.get_lookup_file(inps.lookup_file)
            if not inps.lookup_file:
                raise Exception('No lookup file found! Can not use --tight option without it.')

            atr_lut = readfile.read_attribute(inps.lookup_file)
            if 'Y_FIRST' in atr_lut.keys():
                rg_lut = readfile.read(inps.lookup_file, datasetName='range')[0]
                rg_unique, rg_pos = np.unique(rg_lut, return_inverse=True)
                idx_row, idx_col = np.where(rg_lut != rg_unique[np.bincount(rg_pos).argmax()])
                pix_box = (np.min(idx_col)-10, np.min(idx_row)-10,
                           np.max(idx_col)+10, np.max(idx_row)+10)
                geo_box = box_pixel2geo(pix_box, atr_lut)
                del rg_lut
            else:
                lat = readfile.read(inps.lookup_file, datasetName='latitude')[0]
                lon = readfile.read(inps.lookup_file, datasetName='longitude')[0]
                geo_box = (np.nanmin(lon), np.nanmax(lat),
                           np.nanmax(lon), np.nanmin(lat))
                pix_box = None
                del lat, lon
        else:
            raise Exception('No subset inputs found!')

        # Update subset_y/x/lat/lon
        inps = subset_box2inps(inps, pix_box, geo_box)
    return inps


###########################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    inps.file = ut.get_file_list(inps.file)
    print('number of input files: ({})\n{}'.format(len(inps.file), inps.file))

    inps = read_aux_subset2inps(inps)

    subset_file_list(inps.file, inps)

    print('Done.')
    return


###########################################################################
if __name__ == '__main__':
    main()
