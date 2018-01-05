#! /usr/bin/env python2
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2017, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################


import os
import sys
import time
import argparse
import warnings

import h5py
import numpy as np
from scipy.interpolate import RegularGridInterpolator as RGI
import scipy.spatial.qhull as qhull


import pysar._datetime as ptime
import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._pysar_utilities as ut
from pysar._readfile import multi_group_hdf5_file, multi_dataset_hdf5_file, single_dataset_hdf5_file


############################ Geocoded with lut in geo coord #########################
def update_attribute_geo_lut(atr_rdr, atr_lut, print_msg=True):
    '''Get attributes in geo coord from atr_rdr dict and atr_lut dict
    Inputs:
        atr_rdr : dict, attributes of file in radar coord
        atr_lut : dict, attributes of mapping transformation file
        print_msg : bool, print out message or not
    Output:
        atr : dict, attributes of output file in geo coord.
    '''

    # copy atr_rdr
    atr = dict()
    for key, value in atr_rdr.iteritems():
        atr[key] = str(value)

    atr['FILE_LENGTH'] = atr_lut['FILE_LENGTH']
    atr['WIDTH']   = atr_lut['WIDTH']
    atr['Y_FIRST'] = atr_lut['Y_FIRST']
    atr['X_FIRST'] = atr_lut['X_FIRST']
    atr['Y_STEP']  = atr_lut['Y_STEP']
    atr['X_STEP']  = atr_lut['X_STEP']
    try:    atr['Y_UNIT'] = atr_lut['Y_UNIT']
    except: atr['Y_UNIT'] = 'degrees'
    try:    atr['X_UNIT'] = atr_lut['X_UNIT']
    except: atr['X_UNIT'] = 'degrees'

    # Reference point from y/x to lat/lon
    if 'ref_y' in atr_rdr.keys() and 'ref_x' in atr_rdr.keys():
        ref_x_rdr = np.array(int(atr_rdr['ref_x']))
        ref_y_rdr = np.array(int(atr_rdr['ref_y']))
        trans_file = atr_lut['FILE_PATH']
        ref_lat, ref_lon = ut.radar2glob(ref_y_rdr, ref_x_rdr, trans_file, atr_rdr, print_msg=False)[0:2]
        if ~np.isnan(ref_lat) and ~np.isnan(ref_lon):
            ref_y = np.rint((ref_lat - float(atr['Y_FIRST'])) / float(atr['Y_STEP']))
            ref_x = np.rint((ref_lon - float(atr['X_FIRST'])) / float(atr['X_STEP']))
            atr['ref_lat'] = str(ref_lat)
            atr['ref_lon'] = str(ref_lon)
            atr['ref_y'] = str(int(ref_y))
            atr['ref_x'] = str(int(ref_x))
            if print_msg:
                print 'update ref_lat/lon/y/x'
        else:
            warnings.warn("original reference pixel is out of .trans file's coverage. Continue.")
            try: atr.pop('ref_y')
            except: pass
            try: atr.pop('ref_x')
            except: pass
            try: atr.pop('ref_lat')
            except: pass
            try: atr.pop('ref_lon')
            except: pass
    return atr


def geocode_file_geo_lut(fname, lookup_file, fname_out, inps):
    '''Geocode file using ROI_PAC/Gamma lookup table file.
    Related module: scipy.interpolate.RegularGridInterpolator

    Inputs:
        fname      : string, file to be geocoded
        lookup_file   : string, optional, lookup table file genereated by ROIPAC or Gamma
                     i.e. geomap_4rlks.trans           from ROI_PAC
                          sim_150911-150922.UTM_TO_RDC from Gamma
        interp_method     : string, optional, interpolation/resampling method, supporting nearest, linear
        fill_value : value used for points outside of the interpolation domain.
        fname_out  : string, optional, output geocoded filename
    Output:
        fname_out  : string, optional, output geocoded filename
    '''

    start = time.time()
    ## Default Inputs and outputs
    if not fname_out:
        fname_out = os.path.dirname(os.path.abspath(fname))+'/geo_'+os.path.basename(fname)

    ##### Interpolate value on irregular radar coordinates (from lookup table file value)
    ##### with known value on regular radar coordinates (from radar file attribute)
    ## Grid/regular coordinates from row/column number in radar file
    print '------------------------------------------------------'
    print 'geocoding file: '+fname
    atr_rdr = readfile.read_attribute(fname)
    len_rdr = int(atr_rdr['FILE_LENGTH'])
    wid_rdr = int(atr_rdr['WIDTH'])
    pts_old = (np.arange(len_rdr), np.arange(wid_rdr))

    ## Irregular coordinates from data value in lookup table
    print 'reading lookup table file: '+lookup_file
    atr_lut = readfile.read_attribute(lookup_file)
    rg = readfile.read(lookup_file, epoch='range')[0]
    az = readfile.read(lookup_file, epoch='azimuth')[0]
    len_geo = int(atr_lut['FILE_LENGTH'])
    wid_geo = int(atr_lut['WIDTH'])

    # adjustment if input radar file has been subseted.
    if 'subset_x0' in atr_rdr.keys():
        x0 = float(atr_rdr['subset_x0'])
        y0 = float(atr_rdr['subset_y0'])
        rg -= x0
        az -= y0
        print '\tinput radar coord file has been subsetted, adjust lookup table value'

    # extract pixels only available in radar file (get ride of invalid corners)
    idx = (az>0.0)*(az<=len_rdr)*(rg>0.0)*(rg<=wid_rdr)
    pts_new = np.hstack((az[idx].reshape(-1,1), rg[idx].reshape(-1,1)))
    del az, rg

    print 'geocoding using scipy.interpolate.RegularGridInterpolator ...'
    data_geo = np.empty((len_geo, wid_geo))
    data_geo.fill(inps.fill_value)
    k = atr_rdr['FILE_TYPE']
    ##### Multiple Dataset File
    if k in multi_group_hdf5_file+multi_dataset_hdf5_file:
        h5 = h5py.File(fname,'r')
        epoch_list = sorted(h5[k].keys())
        epoch_num = len(epoch_list)
        prog_bar = ptime.progress_bar(maxValue=epoch_num)

        h5out = h5py.File(fname_out,'w')
        group = h5out.create_group(k)
        print 'writing >>> '+fname_out

        if k in multi_dataset_hdf5_file:
            print 'number of datasets: '+str(epoch_num)
            for i in range(epoch_num):
                date = epoch_list[i]
                data = h5[k].get(date)[:]

                RGI_func = RGI(pts_old, data, method=inps.interp_method,\
                               bounds_error=False, fill_value=inps.fill_value)
                data_geo[idx] = RGI_func(pts_new)

                dset = group.create_dataset(date, data=data_geo, compression='gzip')
                prog_bar.update(i+1, suffix=date)
            prog_bar.close()

            print 'update attributes'
            atr = update_attribute_geo_lut(atr_rdr, atr_lut)
            for key,value in atr.iteritems():
                group.attrs[key] = value

        elif k in multi_group_hdf5_file:
            print 'number of interferograms: '+str(epoch_num)
            date12_list = ptime.list_ifgram2date12(epoch_list)
            for i in range(epoch_num):
                ifgram = epoch_list[i]
                data = h5[k][ifgram].get(ifgram)[:]

                RGI_func = RGI(pts_old, data, method=inps.interp_method,\
                               bounds_error=False, fill_value=inps.fill_value)
                data_geo[idx] = RGI_func(pts_new)

                gg = group.create_group(ifgram)
                dset = gg.create_dataset(ifgram, data=data_geo, compression='gzip')

                atr = update_attribute_geo_lut(h5[k][ifgram].attrs, atr_lut, print_msg=False)
                for key, value in atr.iteritems():
                    gg.attrs[key] = value
                prog_bar.update(i+1, suffix=date12_list[i])
        h5.close()
        h5out.close()

    ##### Single Dataset File
    else:
        print 'reading '+fname
        data = readfile.read(fname)[0]
        RGI_func = RGI(pts_old, data, method=inps.interp_method,\
                       bounds_error=False, fill_value=inps.fill_value)
        data_geo[idx] = RGI_func(pts_new)

        print 'update attributes'
        atr = update_attribute_geo_lut(atr_rdr, atr_lut)

        print 'writing >>> '+fname_out
        writefile.write(data_geo, atr, fname_out)

    del data_geo
    print 'finished writing file: %s' % (fname_out)
    s = time.time()-start;  m, s = divmod(s, 60);  h, m = divmod(m, 60)
    print 'Time used: %02d hours %02d mins %02d secs' % (h, m, s)
    return fname_out


############################ Geocoded with lut in radar coord #######################
## Reference: https://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-
## multiple-interpolations-between-two-irregular-grids
def interp_weights(xy, uv,d=2):
    '''calculate triangulation and coordinates transformation using qhull.Delaunay
    1) Triangulate the irregular grid coordinates xy;
    2) For each point in the new grid uv, search which simplex does it lay
    3) Calculate barycentric coordinates with respect to the vertices of enclosing simplex
    '''
    tri = qhull.Delaunay(xy)
    simplex = tri.find_simplex(uv)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uv - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))


def interpolate(values, vtx, wts, fill_value=np.nan):
    '''Interpolate values on new points'''
    ret = np.einsum('nj,nj->n', np.take(values, vtx), wts)
    ret[np.any(wts < 0, axis=1)] = fill_value
    return ret


def update_attribute_radar_lut(atr_rdr, inps, lat=None, lon=None, print_msg=True):
    '''Get attributes in geo coord from atr_rdr dict and geo_data matrix
    Inputs:
        atr_rdr - dict, attribute of file in radar coord
        inps    - Namespace, including items of the following:
                  lat0/lon0
                  lat_step/lon_step
                  lat_num/lon_num
        lat/lon - 2D np.array of lat/lon value
    Output:
        atr - dict, attributes of output file in geo coord.
    '''
    # copy atr_rdr
    atr = dict()
    for key, value in atr_rdr.iteritems():
        atr[key] = str(value)

    atr['FILE_LENGTH'] = str(inps.lat_num)
    atr['WIDTH'] = str(inps.lon_num)
    atr['Y_FIRST'] = str(inps.lat0)
    atr['X_FIRST'] = str(inps.lon0)
    atr['Y_STEP'] = str(inps.lat_step)
    atr['X_STEP'] = str(inps.lon_step)
    atr['Y_UNIT'] = 'degrees'
    atr['X_UNIT'] = 'degrees'

    ##Reference pixel
    if ('ref_y' in atr_rdr.keys() and lat is not None and\
        'ref_x' in atr_rdr.keys() and lon is not None):
        length_rdr = int(atr_rdr['FILE_LENGTH'])
        width_rdr = int(atr_rdr['WIDTH'])
        ref_y_rdr = int(atr_rdr['ref_y'])
        ref_x_rdr = int(atr_rdr['ref_x'])
        ref_lat = lat[ref_y_rdr, ref_x_rdr]
        ref_lon = lon[ref_y_rdr, ref_x_rdr]

        ref_y = int(np.rint((ref_lat - inps.lat0) / inps.lat_step))
        ref_x = int(np.rint((ref_lon - inps.lon0) / inps.lon_step))
        if 0 <= ref_y <= inps.lat_num and 0 <= ref_x <= inps.lon_num:
            ref_lat = inps.lat0 + ref_y * inps.lat_step
            ref_lon = inps.lon0 + ref_x * inps.lon_step
            atr['ref_lat'] = str(ref_lat)
            atr['ref_lon'] = str(ref_lon)
            atr['ref_y'] = str(ref_y)
            atr['ref_x'] = str(ref_x)
            if print_msg:
                print 'update ref_lat/lon/y/x'
        else:
            warnings.warn("original reference pixel is out of lookup file's coverage. Continue.")
            try: atr.pop('ref_y')
            except: pass
            try: atr.pop('ref_x')
            except: pass
            try: atr.pop('ref_lat')
            except: pass
            try: atr.pop('ref_lon')
            except: pass
    return atr


def geocode_file_radar_lut(fname, lookup_file, fname_out=None, inps=None):
    '''Geocode file using lookup table file in radar coordinates (isce).
    Two solutions:
    1) scipy.interpolate.griddata, with a speed up solution from Jaime and Jeff (Stack Overflow)
        https://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-multiple-interpo
        lations-between-two-irregular-grids
    2) matplotlib.tri, interpolation from triangular grid to quad grid, which is much slower than 1).

    Inputs:
        fname       : string, file to be geocoded
        lookup_file : string, lookup table file, geometryRadar.h5
        fname_out   : string, optional, output geocoded filename
        inps        : namespace, object with the following items:
                      interp_method : string, interpolation/resampling method, supporting linear
                      fill_value    : value used for points outside of the interpolation domain
    Output:
        fname_out  : string, optional, output geocoded filename
    '''
    start = time.time()
    ## Default Inputs and outputs
    if not inps:
        inps = cmdLineParse()

    if inps.interp_method != 'linear':
        print 'ERROR: Supported interpolation method: linear'
        print 'Input method is '+inps.interp_method
        sys.exit(-1)

    if not fname_out:
        fname_out = os.path.dirname(os.path.abspath(fname))+'/geo_'+os.path.basename(fname)

    ## Read lookup table file
    atr_rdr = readfile.read_attribute(fname)
    length = int(atr_rdr['FILE_LENGTH'])
    width = int(atr_rdr['WIDTH'])
    print 'reading lookup table file '+lookup_file
    lat = readfile.read(lookup_file, epoch='latitude')[0]
    lon = readfile.read(lookup_file, epoch='longitude')[0]

    #####Prepare output pixel grid: lat/lon range and step
    if os.path.isfile(inps.lalo_step):
        print 'use file %s as reference for output grid lat/lon range and step' % (inps.lalo_step)
        atr_ref = readfile.read_attribute(inps.lalo_step)
        inps.lat_step = float(atr_ref['Y_STEP'])
        inps.lon_step = float(atr_ref['X_STEP'])
        inps.lat_num = int(atr_ref['FILE_LENGTH'])
        inps.lon_num = int(atr_ref['WIDTH'])
        inps.lat0 = float(atr_ref['Y_FIRST'])
        inps.lon0 = float(atr_ref['X_FIRST'])
        inps.lat1 = inps.lat0 + inps.lat_step*inps.lat_num
        inps.lon1 = inps.lon0 + inps.lon_step*inps.lon_num
    else:
        try:
            inps.lat_step = -1*abs(float(inps.lalo_step))
            inps.lon_step = abs(float(inps.lalo_step))
            inps.lat0 = np.nanmax(lat)
            inps.lat1 = np.nanmin(lat)
            inps.lon0 = np.nanmin(lon)
            inps.lon1 = np.nanmax(lon)
            inps.lat_num = int((inps.lat1-inps.lat0)/inps.lat_step)
            inps.lon_num = int((inps.lon1-inps.lon0)/inps.lon_step)
            inps.lat_step = (inps.lat1 - inps.lat0)/inps.lat_num
            inps.lon_step = (inps.lon1 - inps.lon0)/inps.lon_num
        except ValueError:
            print 'Input lat/lon step is neither a float number nor a file in geo-coord, please try again.'

    print 'output lat range: %f - %f' % (inps.lat0, inps.lat1)
    print 'output lon range: %f - %f' % (inps.lon0, inps.lon1)
    print 'output lat_step : %f' % (inps.lat_step)
    print 'output lon_step : %f' % (inps.lon_step)
    print 'input  file size in   y/x  : %d/%d' % (length, width)
    print 'output file size in lat/lon: %d/%d' % (inps.lat_num, inps.lon_num)

    grid_lat, grid_lon = np.mgrid[inps.lat0:inps.lat1:inps.lat_num*1j,\
                                  inps.lon0:inps.lon1:inps.lon_num*1j]


    ##### Interpolate value on regular geo coordinates (from lookup table file attributes, 2D ndarray)
    ##### with known value on irregular geo coordinates (from lookup table file value, tuple of ndarray of float)

    ##Solution 1 - qhull
    print 'calculate triangulation and coordinates transformation using scipy.spatial.qhull.Delaunay ...'
    pts_old = np.hstack((lat.reshape(-1,1), lon.reshape(-1,1)))
    pts_new = np.hstack((grid_lat.reshape(-1,1), grid_lon.reshape(-1,1)))
    vtx, wts = interp_weights(pts_old, pts_new)
    del pts_old, pts_new, grid_lat, grid_lon

    ##Solution 2 - matplotlib.tri
    #triang = mtri.Triangulation(lat.flatten(),lon.flatten())

    data_geo = np.empty((inps.lat_num, inps.lon_num)).flatten()
    data_geo.fill(inps.fill_value)
    k = atr_rdr['FILE_TYPE']
    ##### Multiple Dataset File
    if k in multi_group_hdf5_file+multi_dataset_hdf5_file:
        h5 = h5py.File(fname,'r')
        epoch_list = sorted(h5[k].keys())
        epoch_num = len(epoch_list)
        prog_bar = ptime.progress_bar(maxValue=epoch_num)

        h5out = h5py.File(fname_out,'w')
        group = h5out.create_group(k)
        print 'writing >>> '+fname_out

        if k in multi_dataset_hdf5_file:
            print 'number of acquisitions: '+str(epoch_num)
            for i in range(epoch_num):
                date = epoch_list[i]
                data = h5[k].get(date)[:]

                data_geo = interpolate(data.flatten(), vtx, wts).reshape(inps.lat_num, inps.lon_num)

                dset = group.create_dataset(date, data=data_geo, compression='gzip')
                prog_bar.update(i+1, suffix=date)
            prog_bar.close()

            print 'update attributes'
            atr = update_attribute_radar_lut(atr_rdr, inps, lat, lon)
            for key,value in atr.iteritems():
                group.attrs[key] = value

        elif k in multi_group_hdf5_file:
            print 'number of interferograms: '+str(epoch_num)
            date12_list = ptime.list_ifgram2date12(epoch_list)
            for i in range(epoch_num):
                ifgram = epoch_list[i]
                data = h5[k][ifgram].get(ifgram)[:]

                data_geo = interpolate(data.flatten(), vtx, wts).reshape(inps.lat_num, inps.lon_num)

                gg = group.create_group(ifgram)
                dset = gg.create_dataset(ifgram, data=data_geo, compression='gzip')

                atr = update_attribute_radar_lut(h5[k][ifgram].attrs, inps, lat, lon, print_msg=False)
                for key, value in atr.iteritems():
                    gg.attrs[key] = value
                prog_bar.update(i+1, suffix=date12_list[i])
        h5.close()
        h5out.close()

    ##### Single Dataset File
    else:
        print 'reading '+fname
        data = readfile.read(fname)[0]

        ##Solution 1 - qhull
        data_geo = interpolate(data.flatten(), vtx, wts).reshape(inps.lat_num, inps.lon_num)

        ###Solution 2 - matplotlib.tri
        #interp_lin = mtri.LinearTriInterpolator(triang, data.flatten())
        #data_geo = interp_lin(grid_lat.flatten(), grid_lon.flatten())
        #interp_cubic = mtri.CubicTriInterpolator(triang, data, kind='geom')
        #data_geo = interp_cubic(grid_lat, grid_lon)

        print 'update attributes'
        atr = update_attribute_radar_lut(atr_rdr, inps, lat, lon)

        print 'writing >>> '+fname_out
        writefile.write(data_geo, atr, fname_out)

    del data_geo, vtx, wts
    print 'finished writing file: %s' % (fname_out)
    s = time.time()-start;  m, s = divmod(s, 60);  h, m = divmod(m, 60)
    print 'Time used: %02d hours %02d mins %02d secs' % (h, m, s)
    return fname_out


def geocode_file(fname, lookup_file, fname_out, inps):
    '''Geocode input file with lookup table file'''
    atr = readfile.read_attribute(lookup_file)
    if 'Y_FIRST' in atr.keys():
        if not inps.interp_method:
            inps.interp_method = 'nearest'
        print 'lookup table in geo coordinates: '+lookup_file
        print 'interpolation method: '+inps.interp_method
        fname_out = geocode_file_geo_lut(fname, lookup_file, fname_out, inps)
    else:
        if not inps.interp_method:
            inps.interp_method = 'linear'
        print 'lookup table in radar coordinates: '+lookup_file
        print 'interpolation method: '+inps.interp_method
        fname_out = geocode_file_radar_lut(fname, lookup_file, fname_out, inps)
    return fname_out


def read_template2inps(template_file, inps):
    '''Read input template options into Namespace inps'''
    print 'read input option from template file: '+template_file
    if not inps:
        inps = cmdLineParse()

    template = readfile.read_template(template_file)
    key_list = template.keys()

    # Coherence-based network modification
    prefix = 'pysar.geocode.'

    key = prefix+'resolution'
    if key in key_list:
        value = template[key]
        if value not in ['auto']:
            inps.lalo_step = value

    return inps


######################################################################################
TEMPLATE='''template:
## For data processed by ROI_PAC/Gamma, output resolution for geocoded file is the same as their lookup table file.
## For data processed by ISCE/Doris, output resolution is assign by user with resolution option:
## 1) float number - resolution in degree, 0.001 by default, around 100 m on equator
## 2) file name    - use the resolution from a file in geo coordinates, e.g. demGeo.h5
pysar.geocode            = auto  #[yes / no], auto for yes
pysar.geocode.resolution = auto  #[0.0-inf / filename], auto for 0.001 (~100 m), output resolution for ISCE processor
'''

EXAMPLE='''example:
  geocode.py  velocity.h5
  geocode.py  velocity.h5 temporalCoherence.h5 timeseries_ECMWF_demErr_refDate.h5
  geocode.py  velocity.h5 temporalCoherence.h5 timeseries_ECMWF_demErr_refDate.h5 -t pysarApp_template.txt

  geocode.py  101120-110220.unw   -l geomap_4rlks.trans
  geocode.py  velocity.h5         -l sim_150911-150922.UTM_TO_RDC
  geocode.py  coherence.h5        -l geometryRadar.h5   --lalo-step 0.0003333
  geocode.py  unwrapIfgram.h5     -l geometryRadar.h5   --lalo-step demGeo_tight.h5
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Geocode using lookup table',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=TEMPLATE+'\n'+EXAMPLE)

    parser.add_argument('file', nargs='+', help='File(s) to be geocoded')
    parser.add_argument('-l','--lookup', dest='lookup_file', help='Lookup table file generated by InSAR processors.')
    parser.add_argument('--template','-t', dest='template_file', help='Template file with geocoding options.')
    parser.add_argument('-i','--interpolate', dest='interp_method',\
                        choices={'nearest','linear'},\
                        help='interpolation/resampling method. Default:\n'+\
                             'nearest - lookup table in geo   coord (roipac, gamma)\n'+\
                             'linear  - lookup table in radar coord (isce)')
    parser.add_argument('--fill', dest='fill_value', type=float, default=np.nan,\
                        help='Value used for points outside of the interpolation domain.\n'+\
                             'Default: np.nan')
    parser.add_argument('--no-parallel',dest='parallel',action='store_false',default=True,\
                        help='Disable parallel processing. Diabled auto for 1 input file.')
    parser.add_argument('-o','--output', dest='outfile', help="output file name. Default: add prefix 'geo_'")
    parser.add_argument('--lalo-step', dest='lalo_step', type=str, default='0.001',\
                        help='output pixel size in degree, for lookup table in radar coord (isce)\n'+\
                             'Support 2 types of input:\n'+\
                             '1) float number - as the output pixel size in degree. Default: 0.001 (~100m)\n'+\
                             '2) DEM file in geo coord - as reference for a) geo coverage and b) pixel size')

    inps = parser.parse_args()
    return inps


######################################################################################
def main(argv):
    inps = cmdLineParse()
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    inps.file = ut.get_file_list(inps.file)
    print 'number of files to geocode: '+str(len(inps.file))
    print inps.file
    if len(inps.file) > 1:
        inps.outfile = None
    print 'fill_value: '+str(inps.fill_value)

    ##Check Lookup table
    inps.lookup_file = ut.get_lookup_file(inps.lookup_file)
    if not inps.lookup_file:
        sys.exit('No lookup table found! Can not geocode without it.')

    # check outfile and parallel option
    if inps.parallel:
        num_cores, inps.parallel, Parallel, delayed = ut.check_parallel(len(inps.file))

    #####
    if inps.parallel:
        Parallel(n_jobs=num_cores)(delayed(geocode_file)\
                                   (fname, inps.lookup_file, inps.outfile, inps) for fname in inps.file)
    else:
        for fname in inps.file:
            geocode_file(fname, inps.lookup_file, inps.outfile, inps)

    print 'Done.'
    return


######################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
