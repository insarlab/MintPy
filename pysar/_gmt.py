#! /usr/bin/env python2
############################################################
# _gmt.py is part of the GIAnT v1.0 
############################################################

import scipy.io.netcdf as netcdf
import numpy as np

'''This is based on the gdal2grd.py script found at http://http://www.vso.cape.com/~nhv/files/python/gdal/gdal2grd.py'''
def write_gmt_simple(lons, lats, z, fname, title='default', name='z', scale=1.0, offset=0, units='meters'):
    '''Writes a simple GMT grd file with one array.
    
    .. Args:
        
        * lons     -> 1D Array of lon values
        * lats     -> 1D Array of lat values
        * z        -> 2D slice to be saved
        * fname    -> Output file name
        
    .. Kwargs:
        
        * title    -> Title for the grd file
        * name     -> Name of the field in the grd file
        * scale    -> Scale value in the grd file
        * offset   -> Offset value in the grd file
        
    .. Returns:
        
        * None'''
    fid = netcdf.netcdf_file(fname,'w')

    ####Create a dimension variable
    fid.createDimension('side',2)
    fid.createDimension('xysize',np.prod(z.shape))

    ####Range variables
    fid.createVariable('x_range','d',('side',))
    fid.variables['x_range'].units = 'degrees'

    fid.createVariable('y_range','d',('side',))
    fid.variables['y_range'].units = 'degrees'

    fid.createVariable('z_range','d',('side',))
    fid.variables['z_range'].units = units

    #####Spacing
    fid.createVariable('spacing','d',('side',))
    fid.createVariable('dimension','i4',('side',))

    fid.createVariable('z','f',('xysize',))
    fid.variables['z'].long_name = name
    fid.variables['z'].scale_factor = scale
    fid.variables['z'].add_offset = offset
    fid.variables['z'].node_offset=0

    fid.title = title
    fid.source = 'PySAR v1.2'

    #####Filling in the actual data
    fid.variables['x_range'][0] = lons[0]
    fid.variables['x_range'][1] = lons[-1]
    fid.variables['spacing'][0] = lons[1]-lons[0]

    fid.variables['y_range'][0] = lats[0]
    fid.variables['y_range'][1] = lats[-1]
    fid.variables['spacing'][1] = lats[1]-lats[0]

    #####Range
    zmin = np.nanmin(z)
    zmax = np.nanmax(z)

    fid.variables['z_range'][0] = zmin
    fid.variables['z_range'][1] = zmax

    fid.variables['dimension'][:] = z.shape[::-1]
    fid.variables['z'][:] = np.flipud(z).flatten()
    fid.close()

############################################################
# Program is part of GIAnT v1.0                            #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
############################################################
