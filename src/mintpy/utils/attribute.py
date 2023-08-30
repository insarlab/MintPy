"""Update attributes for operations."""
#############################################################
# Program is part of MintPy                                 #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi          #
# Author: Zhang Yunjun, Dec 2020                            #
#############################################################
# Recommend import:
#   from mintpy.utils import attribute as attr


import warnings

import numpy as np

from mintpy.objects.coord import coordinate
from mintpy.utils import readfile


def update_attribute4resize(atr_in, resize2shape, print_msg=True):
    """update input dictionary of attributes due to resizing

    Parameters: atr_in       - dict, input dictionary of attributes
                resize2shape - tuple of 2 int, for the resized shape
    Returns:    atr          - dict, updated dictionary of attributes
    """
    vprint = print if print_msg else lambda *args, **kwargs: None
    # make a copy of original meta dict
    atr = {**atr_in}

    yscale = int(atr['LENGTH']) / resize2shape[0]
    xscale = int(atr['WIDTH']) / resize2shape[1]
    vprint(f'output data in size: {resize2shape[0]}, {resize2shape[1]}')

    atr['LENGTH'] = resize2shape[0]
    atr['WIDTH'] = resize2shape[1]
    atr['ALOOKS'] = np.rint(int(atr.get('ALOOKS', 1)) * yscale).astype(int)
    atr['RLOOKS'] = np.rint(int(atr.get('RLOOKS', 1)) * xscale).astype(int)
    vprint('update LENGTH, WIDTH, Y/XMIN/MAX, A/RLOOKS')

    if 'AZIMUTH_PIXEL_SIZE' in atr.keys():
        atr['AZIMUTH_PIXEL_SIZE'] = float(atr['AZIMUTH_PIXEL_SIZE']) * yscale
        vprint('update AZIMUTH_PIXEL_SIZE')

    if 'RANGE_PIXEL_SIZE' in atr.keys():
        atr['RANGE_PIXEL_SIZE'] = float(atr['RANGE_PIXEL_SIZE']) * xscale
        vprint('update RANGE_PIXEL_SIZE')

    if 'NCORRLOOKS' in atr.keys():
        atr['NCORRLOOKS'] = float(atr['NCORRLOOKS']) * yscale * xscale
        vprint('update NCORRLOOKS')

    return atr


def update_attribute4multilook(atr_in, lks_y, lks_x, box=None, print_msg=True):
    """update input dictionary of attributes due to multilooking

    Parameters: atr_in - dict, input dictionary of attributes
                lks_y  - int, number of looks in y/row/azimuth  direction
                lks_x  - int, number of looks in x/column/range direction
                box    - tuple of 4 int indicating (x0, y0, x1, y1)
                         if --margin option is used in multilook.py
    Returns:    atr    - dict, updated dictionary of attributes
    """
    vprint = print if print_msg else lambda *args, **kwargs: None

    # make a copy of original meta dict
    atr = {**atr_in}

    if box is None:
        box = (0, 0, int(atr['WIDTH']), int(atr['LENGTH']))
    length, width = box[3] - box[1], box[2] - box[0]

    length_mli = length // lks_y
    width_mli = width // lks_x
    vprint(f'output data in size: {length_mli}, {width_mli}')

    # Update attributes
    atr['LENGTH'] = str(length_mli)
    atr['WIDTH'] = str(width_mli)
    atr['XMIN'] = str(box[0])
    atr['YMIN'] = str(box[1])
    atr['XMAX'] = str(width_mli - 1 + box[0])
    atr['YMAX'] = str(length_mli - 1 + box[1])
    atr['RLOOKS'] = str(int(atr.get('RLOOKS', '1')) * lks_x)
    atr['ALOOKS'] = str(int(atr.get('ALOOKS', '1')) * lks_y)
    vprint('update LENGTH, WIDTH, Y/XMIN/MAX, A/RLOOKS')

    if 'Y_STEP' in atr.keys():
        atr['Y_STEP'] = str(lks_y * float(atr['Y_STEP']))
        atr['X_STEP'] = str(lks_x * float(atr['X_STEP']))
        vprint('update Y/X_STEP')

    if 'AZIMUTH_PIXEL_SIZE' in atr.keys():
        atr['AZIMUTH_PIXEL_SIZE'] = str(lks_y * float(atr['AZIMUTH_PIXEL_SIZE']))
        vprint('update AZIMUTH_PIXEL_SIZE')

    if 'RANGE_PIXEL_SIZE' in atr.keys():
        atr['RANGE_PIXEL_SIZE'] = str(lks_x * float(atr['RANGE_PIXEL_SIZE']))
        vprint('update RANGE_PIXEL_SIZE')

    if 'REF_Y' in atr.keys():
        atr['REF_Y'] = str( (int(atr['REF_Y']) - box[1]) // lks_y )
        atr['REF_X'] = str( (int(atr['REF_X']) - box[0]) // lks_x )
        vprint('update REF_Y/X')

    if 'SUBSET_XMIN' in atr.keys():
        atr['SUBSET_YMIN'] = str( (int(atr['SUBSET_YMIN']) - box[1]) // lks_y )
        atr['SUBSET_YMAX'] = str( (int(atr['SUBSET_YMAX']) - box[1]) // lks_y )
        atr['SUBSET_XMIN'] = str( (int(atr['SUBSET_XMIN']) - box[0]) // lks_x )
        atr['SUBSET_XMAX'] = str( (int(atr['SUBSET_XMAX']) - box[0]) // lks_x )
        vprint('update SUBSET_XMIN/XMAX/YMIN/YMAX')
    return atr


def update_attribute4geo2radar(atr_in, shape2d=None, res_obj=None, print_msg=True):
    """update input dictionary of attributes due to resampling from geo to radar coordinates

    Parameters: atr_in  - dict, input dictionary of attributes
                # combination 1
                shape2d - tuple of 2 int in (length, width)
                # combination 2
                res_obj   - mintpy.objects.resample.resample object
    Returns:    atr     - dict, updated dictionary of attributes
    """
    # make a copy of original meta dict
    atr = {**atr_in}

    # grab info from res_obj
    if res_obj is not None:
        shape2d = (res_obj.length, res_obj.width)

    # update shape
    atr['LENGTH'] = shape2d[0]
    atr['WIDTH'] = shape2d[1]

    # remove geo-coord related metadata
    for key in ['Y_FIRST', 'Y_STEP', 'Y_UNIT',
                'X_FIRST', 'X_STEP', 'X_UNIT',
                'REF_Y', 'REF_X', 'REF_LAT', 'REF_LON']:
        if key in atr.keys():
            atr.pop(key)
    return atr


def update_attribute4radar2geo(atr_in, shape2d=None, lalo_step=None, SNWE=None, lut_file=None,
                               res_obj=None, print_msg=True):
    """update input dictionary of attributes due to resampling from radar to geo coordinates

    Parameters: atr_in  - dict, input dictionary of attributes
                # combination 1
                shape2d - tuple of 2 int in (length, width)
                lalo_step - tuple of 2 float, step size in lat/lon direction
                SNWE      - tuple of 4 float
                lut_file  - str, path of lookup table file
                # combination 2
                res_obj   - mintpy.objects.resample.resample object
    Returns:    atr     - dict, updated dictionary of attributes
    """
    # make a copy of original meta dict
    atr = {**atr_in}

    # grab info from res_obj
    if res_obj is not None:
        shape2d = (res_obj.length, res_obj.width)
        lalo_step = res_obj.lalo_step
        SNWE = res_obj.SNWE
        lut_file = res_obj.lut_file

    atr['LENGTH'] = shape2d[0]
    atr['WIDTH'] = shape2d[1]
    atr['Y_STEP'] = lalo_step[0]
    atr['X_STEP'] = lalo_step[1]
    atr['Y_FIRST'] = SNWE[1]
    atr['X_FIRST'] = SNWE[2]

    lut_meta = readfile.read_attribute(lut_file)
    atr['Y_UNIT'] = lut_meta.get('Y_UNIT', 'degrees')
    atr['X_UNIT'] = lut_meta.get('X_UNIT', 'degrees')

    # Reference point from y/x to lat/lon
    if 'REF_Y' in atr.keys():
        coord = coordinate(atr_in, lookup_file=lut_file)
        ref_lat, ref_lon = coord.radar2geo(np.array(int(atr['REF_Y'])),
                                           np.array(int(atr['REF_X'])),
                                           print_msg=False)[0:2]
        if ~np.isnan(ref_lat) and ~np.isnan(ref_lon):
            ref_y = int(np.rint((ref_lat - float(atr['Y_FIRST'])) / float(atr['Y_STEP'])))
            ref_x = int(np.rint((ref_lon - float(atr['X_FIRST'])) / float(atr['X_STEP'])))
            atr['REF_LAT'] = str(ref_lat)
            atr['REF_LON'] = str(ref_lon)
            atr['REF_Y'] = str(ref_y)
            atr['REF_X'] = str(ref_x)
            if print_msg:
                print('update REF_LAT/LON/Y/X')
        else:
            warnings.warn("original reference pixel is out of .trans file's coverage. Continue.")
            try:
                atr.pop('REF_Y')
                atr.pop('REF_X')
            except:
                pass
            try:
                atr.pop('REF_LAT')
                atr.pop('REF_LON')
            except:
                pass
    return atr


def update_attribute4subset(atr_in, subset_box, print_msg=True):
    """update input dictionary of attributes due to subset

    Parameters: atr_in     - dict, data attributes to update
                subset_box - 4-tuple of int, subset box defined in (x0, y0, x1, y1)
    Returns:    atr        - dict, updated data attributes
    """
    vprint = print if print_msg else lambda *args, **kwargs: None

    if subset_box is None:
        return atr_in

    sub_x = [subset_box[0], subset_box[2]]
    sub_y = [subset_box[1], subset_box[3]]

    # Update attribute variable
    atr = {**atr_in}
    atr['LENGTH'] = str(sub_y[1]-sub_y[0])
    atr['WIDTH'] = str(sub_x[1]-sub_x[0])
    atr['YMAX'] = str(sub_y[1]-sub_y[0] - 1)
    atr['XMAX'] = str(sub_x[1]-sub_x[0] - 1)
    vprint('update LENGTH, WIDTH, Y/XMAX')

    # Subset attribute
    atr['SUBSET_YMAX'] = str(sub_y[1] + int(atr_in.get('SUBSET_YMIN', '0')))
    atr['SUBSET_YMIN'] = str(sub_y[0] + int(atr_in.get('SUBSET_YMIN', '0')))
    atr['SUBSET_XMAX'] = str(sub_x[1] + int(atr_in.get('SUBSET_XMIN', '0')))
    atr['SUBSET_XMIN'] = str(sub_x[0] + int(atr_in.get('SUBSET_XMIN', '0')))
    vprint(('update/add SUBSET_XMIN/YMIN/XMAX/YMAX: '
            '{x0}/{y0}/{x1}/{y1}').format(x0=atr['SUBSET_XMIN'],
                                          y0=atr['SUBSET_YMIN'],
                                          x1=atr['SUBSET_XMAX'],
                                          y1=atr['SUBSET_YMAX']))

    # Geo coord
    if 'Y_FIRST' in atr.keys():
        atr['Y_FIRST'] = str(float(atr['Y_FIRST']) + sub_y[0]*float(atr['Y_STEP']))
        atr['X_FIRST'] = str(float(atr['X_FIRST']) + sub_x[0]*float(atr['X_STEP']))
        vprint('update Y/X_FIRST')

    # Reference in space
    if 'REF_Y' in atr.keys():
        atr['REF_Y'] = str(int(atr['REF_Y']) - sub_y[0])
        atr['REF_X'] = str(int(atr['REF_X']) - sub_x[0])
        vprint('update REF_Y/X')

    # Starting Range for file in radar coord
    if 'Y_FIRST' not in atr_in.keys():
        try:
            atr['STARTING_RANGE'] = float(atr['STARTING_RANGE'])
            atr['STARTING_RANGE'] += float(atr['RANGE_PIXEL_SIZE'])*sub_x[0]
            vprint('update STARTING_RANGE')
        except:
            pass

    return atr
