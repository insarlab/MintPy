#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import sys
import argparse
import h5py
import numpy as np
import random
from mintpy.objects import timeseries
from mintpy.defaults.template import get_template_content
from mintpy.utils import readfile, writefile, utils as ut


#########################################  Usage  ##############################################
TEMPLATE = get_template_content('reference_point')

NOTE = """note: Reference value cannot be nan, thus, all selected reference point must be:
  a. non zero in mask, if mask is given
  b. non nan  in data (stack)
  
  Priority:
      input reference_lat/lon
      input reference_y/x
      input selection_method
      existing REF_Y/X attributes (can be ignored by --force option)
      default selection methods:
          maxCoherence
          random

  The recommended reference pixel should meets the following criteria:
  1) not in deforming areas
  2) not in areas affected by strong atmospheric turbulence, such as ionospheric streaks
  3) close but outside of deforming area of interest with similar elevation, to minimize
     the spatial correlation effect of atmosspheric delay, especially for shot-wavelength
     deformation (Chaussard et al., 2013; Morales-Rivera et al., 2016)
  4) in high coherent area to minimize the decorrelation effect
"""

EXAMPLE = """example:
  # for ifgramStack file, update metadata only
  # add --write-data to update data matrix value
  reference_point.py  inputs/ifgramStack.h5  -t smallbaselineApp.cfg  -c avgSpatialCoh.h5
  reference_point.py  inputs/ifgramStack.h5 --method manual
  reference_point.py  inputs/ifgramStack.h5 --method random

  # for all the other files, update both metadata and data matrix value
  reference_point.py  091120_100407.unw -y 257    -x 151      -m Mask.h5
  reference_point.py  geo_velocity.h5   -l 34.45  -L -116.23  -m Mask.h5
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Reference to the same pixel in space.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=NOTE+'\n'+TEMPLATE+'\n'+EXAMPLE)

    parser.add_argument('file', type=str, help='file to be referenced.')
    parser.add_argument('-t', '--template', dest='template_file',
                        help='template with reference info')
    parser.add_argument('-m', '--mask', dest='maskFile', help='mask file')

    parser.add_argument('-o', '--outfile', type=str, default=None,
                        help='output file name (default: %(default)s). This option is diabled for ifgramStack file.\n'
                             'None (default) for update data value directly without writing to a new file.\n')

    parser.add_argument('--write-data', dest='write_data', action='store_true',
                        help='(option for ifgramStack file only) update data value, in addition to update metadata.')

    parser.add_argument('--reset', action='store_true',
                        help='remove reference pixel information from attributes in the file')
    parser.add_argument('--force', action='store_true',
                        help='Enforce the re-selection of reference point.')

    coord = parser.add_argument_group('input coordinates')
    coord.add_argument('-y', '--row', dest='ref_y', type=int,
                       help='row/azimuth  number of reference pixel')
    coord.add_argument('-x', '--col', dest='ref_x', type=int,
                       help='column/range number of reference pixel')
    coord.add_argument('-l', '--lat', dest='ref_lat',
                       type=float, help='latitude  of reference pixel')
    coord.add_argument('-L', '--lon', dest='ref_lon',
                       type=float, help='longitude of reference pixel')

    coord.add_argument('-r', '--reference', dest='reference_file',
                       help='use reference/seed info of this file')
    coord.add_argument('--lookup', '--lookup-file', dest='lookup_file',
                       help='Lookup table file from SAR to DEM, i.e. geomap_4rlks.trans\n' +
                            'Needed for radar coord input file with --lat/lon seeding option.')

    parser.add_argument('-c', '--coherence', dest='coherenceFile', default='averageSpatialCoherence.h5',
                        help='use input coherence file to find the pixel with max coherence for reference pixel.')
    parser.add_argument('--min-coherence', dest='minCoherence', type=float, default=0.85,
                        help='minimum coherence of reference pixel for max-coherence method.')
    parser.add_argument('--method', type=str, choices=['maxCoherence', 'manual', 'random'],
                        help='methods to select reference pixel if not given in specific y/x or lat/lon:\n' +
                             'maxCoherence : select pixel with highest coherence value as reference point\n' +
                             '               enabled when there is --coherence option input\n' +
                             'manual       : display stack of input file and manually select reference point\n' +
                             'random       : random select pixel as reference point\n')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    
    atr = readfile.read_attribute(inps.file)
    if atr['FILE_TYPE'] != 'ifgramStack':
        # turn ON wirte_data for non-ifgramStack file by default
        inps.write_data = True

    else:
        # disable --output option for ifgramStack file
        if inps.outfile:
            raise SystemExit('--outfile is disabled for "ifgramStack" input file!')

    return inps


def read_template_file2inps(template_file, inps=None):
    """Read seed/reference info from template file and update input namespace"""
    if not inps:
        inps = cmd_line_parse([''])
    inps_dict = vars(inps)
    template = readfile.read_template(template_file)
    template = ut.check_template_auto_value(template)

    prefix = 'mintpy.reference.'
    key_list = [i for i in list(inps_dict)
                if prefix+i in template.keys()]
    for key in key_list:
        value = template[prefix+key]
        if value:
            if key in ['coherenceFile', 'maskFile']:
                inps_dict[key] = value
            elif key == 'minCoherence':
                inps_dict[key] = float(value)

    key = prefix+'yx'
    if key in template.keys():
        value = template[key]
        if value:
            value = value.replace('[','').replace(']','')
            inps.ref_y, inps.ref_x = [int(i) for i in value.split(',')]

    key = prefix+'lalo'
    if key in template.keys():
        value = template[key]
        if value:
            value = value.replace('[','').replace(']','')
            inps.ref_lat, inps.ref_lon = [float(i) for i in value.split(',')]

    return inps


########################################## Sub Functions #############################################
def nearest(x, tbase, xstep):
    """ find nearest neighbour """
    dist = np.sqrt((tbase - x)**2)
    if min(dist) <= np.abs(xstep):
        indx = dist == min(dist)
    else:
        indx = []
    return indx


def reference_file(inps):
    """Seed input file with option from input namespace
    Return output file name if succeed; otherwise, return None
    """
    if not inps:
        inps = cmd_line_parse([''])
    atr = readfile.read_attribute(inps.file)

    # update_mode
    if (not inps.force 
            and inps.ref_y is not None and inps.ref_y == int(atr.get('REF_Y', -999)) 
            and inps.ref_x is not None and inps.ref_x == int(atr.get('REF_X', -999))):
        print('SAME reference pixel is already selected/saved in file, skip updating.')
        return inps.file

    # Check 1 - stack and its non-nan mask pixel coverage
    # outFile=False --> no avgPhaseVelocity file is generated due to the lack of reference point info.
    # did not use maskConnComp.h5 because not all input dataset has connectComponent info
    stack = ut.temporal_average(inps.file, datasetName='unwrapPhase', updateMode=True, outFile=False)[0]
    mask = np.multiply(~np.isnan(stack), stack != 0.)
    if np.nansum(mask) == 0.0:
        raise ValueError('no pixel found with valid phase value in all datasets.')

    # Check 2 - input ref_y/x: location and validity
    if inps.ref_y is not None and inps.ref_x is not None:
        if mask[inps.ref_y, inps.ref_x] == 0.:
            raise ValueError('reference y/x have nan value in some dataset. Please re-select.')
    else:
        # Find reference y/x
        if inps.method == 'maxCoherence':
            inps.ref_y, inps.ref_x = select_max_coherence_yx(coh_file=inps.coherenceFile,
                                                             mask=mask,
                                                             min_coh=inps.minCoherence)
        elif inps.method == 'random':
            inps.ref_y, inps.ref_x = random_select_reference_yx(mask)
        elif inps.method == 'manual':
            inps = manual_select_reference_yx(stack, inps, mask)

        # Check ref_y/x from auto method
        if inps.ref_y is None or inps.ref_x is None:
            raise ValueError('ERROR: no reference y/x found.')

    # Seeding file with reference y/x
    atrNew = reference_point_attribute(atr, y=inps.ref_y, x=inps.ref_x)
    if not inps.write_data:
        print('Add/update ref_x/y attribute to file: '+inps.file)
        print(atrNew)
        inps.outfile = ut.add_attribute(inps.file, atrNew)

    else:
        if not inps.outfile:
            inps.outfile = inps.file

        k = atr['FILE_TYPE']
        fext = os.path.splitext(inps.file)[1]

        if fext == '.h5':
            if inps.outfile == inps.file:
                print('updating data value without re-writing to a new file')

                if k == 'ifgramStack':
                    with h5py.File(inps.file, 'r+') as f:
                        ds = f['unwrapPhase']
                        for i in range(ds.shape[0]):
                            ds[i, :, :] -= ds[i, inps.ref_y, inps.ref_x]

                        print('update metadata')
                        f.attrs.update(atrNew)

                else:
                    with h5py.File(inps.file, 'r+') as f:
                        ds = f[k]
                        if len(ds.shape) == 3:
                            # 3D matrix
                            for i in range(ds.shape[0]):
                                ds[i, :, :] -= ds[i, inps.ref_y, inps.ref_x]

                        else:
                            # 2D matrix
                            ds[:] -= ds[inps.ref_y, inps.ref_x]

                        print('update metadata')
                        f.attrs.update(atrNew)

            else:
                ## write to a new file
                print('writing the referenced data into file: {}'.format(inps.outfile))

                # 1. read and update data value
                data, atr = readfile.read(inps.file)
                if len(data.shape) == 3:
                    # 3D matrix
                    for i in range(data.shape[0]):
                        data[i, :, :] -= data[i, inps.ref_y, inps.ref_x]

                else:
                    # 2D matrix
                    data -= data[inps.ref_y, inps.ref_x]

                # 2. update metadata
                atr.update(atrNew)

                # 3. write to file
                writefile.write(data, inps.outfile, metadata=atr, ref_file=inps.file)

        else:
            # for binary file, over-write directly
            data = readfile.read(inps.file)[0]
            data -= data[inps.ref_y, inps.ref_x]
            atr.update(atrNew)
            writefile.write(data, out_file=inps.outfile, metadata=atr)

    ut.touch([inps.coherenceFile, inps.maskFile])
    return inps.outfile


###############################################################
def reference_point_attribute(atr, y, x):
    atrNew = dict()
    atrNew['REF_Y'] = str(y)
    atrNew['REF_X'] = str(x)
    coord = ut.coordinate(atr)
    if 'X_FIRST' in atr.keys():
        atrNew['REF_LAT'] = str(coord.yx2lalo(y, coord_type='y'))
        atrNew['REF_LON'] = str(coord.yx2lalo(x, coord_type='x'))
    return atrNew


###############################################################
def manual_select_reference_yx(data, inps, mask=None):
    """
    Input: 
        data4display : 2D np.array, stack of input file
        inps    : namespace, with key 'REF_X' and 'REF_Y', which will be updated
    """
    import matplotlib.pyplot as plt
    print('\nManual select reference point ...')
    print('Click on a pixel that you want to choose as the refernce ')
    print('    pixel in the time-series analysis;')
    print('Then close the displayed window to continue.\n')
    if mask is not None:
        data[mask == 0] = np.nan

    # Mutable object
    # ref_url: http://stackoverflow.com/questions/15032638/how-to-return
    #           -a-value-from-button-press-event-matplotlib
    # Display
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(data)

    # Selecting Point
    def onclick(event):
        if event.button == 1:
            print('click')
            x = int(event.xdata+0.5)
            y = int(event.ydata+0.5)

            if not np.isnan(data[y][x]):
                print('valid input reference y/x: '+str([y, x]))
                inps.ref_y = y
                inps.ref_x = x
                # plt.close(fig)
            else:
                print('\nWARNING:')
                print('The selectd pixel has NaN value in data.')
                print('Try a difference location please.')
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    print('y/x: {}'.format((inps.ref_y, inps.ref_x)))
    return inps


def select_max_coherence_yx(coh_file, mask=None, min_coh=0.85):
    """Select pixel with coherence > min_coh in random"""
    print('random select pixel with coherence > {}'.format(min_coh))
    print('\tbased on coherence file: '+coh_file)
    coh, coh_atr = readfile.read(coh_file)
    if mask is not None:
        coh[mask == 0] = 0.0
    coh_mask = coh >= min_coh
    if np.all(coh_mask == 0.):
        msg = ('No pixel with average spatial coherence > {} '
               'are found for automatic reference point selection!').format(min_coh)
        msg += '\nTry the following:'
        msg += '\n  1) manually specify the reference point using mintpy.reference.yx/lalo option.'
        msg += '\n  2) change mintpy.reference.minCoherence to a lower value.'
        raise RuntimeError(msg)

    y, x = random_select_reference_yx(coh_mask, print_msg=False)
    #y, x = np.unravel_index(np.argmax(coh), coh.shape)
    print('y/x: {}'.format((y, x)))
    return y, x


def random_select_reference_yx(data_mat, print_msg=True):
    nrow, ncol = np.shape(data_mat)
    y = random.choice(list(range(nrow)))
    x = random.choice(list(range(ncol)))
    while data_mat[y, x] == 0:
        y = random.choice(list(range(nrow)))
        x = random.choice(list(range(ncol)))
    if print_msg:
        print('random select pixel\ny/x: {}'.format((y, x)))
    return y, x


###############################################################
def read_reference_file2inps(reference_file, inps=None):
    """Read reference info from reference file and update input namespace"""
    if not inps:
        inps = cmd_line_parse([''])
    atrRef = readfile.read_attribute(inps.reference_file)
    if (inps.ref_y is None or inps.ref_x is None) and 'REF_X' in atrRef.keys():
        inps.ref_y = int(atrRef['REF_Y'])
        inps.ref_x = int(atrRef['REF_X'])
    if (inps.ref_lat is None or inps.ref_lon is None) and 'REF_LON' in atrRef.keys():
        inps.ref_lat = float(atrRef['REF_LAT'])
        inps.ref_lon = float(atrRef['REF_LON'])
    return inps


def read_reference_input(inps):
    atr = readfile.read_attribute(inps.file)
    length = int(atr['LENGTH'])
    width = int(atr['WIDTH'])
    inps.go_reference = True

    if inps.reset:
        remove_reference_pixel(inps.file)
        inps.outfile = inps.file
        inps.go_reference = False
        return inps

    print('-'*50)
    # Check Input Coordinates
    # Read ref_y/x/lat/lon from reference/template
    # priority: Direct Input > Reference File > Template File
    if inps.template_file:
        print('reading reference info from template: '+inps.template_file)
        inps = read_template_file2inps(inps.template_file, inps)
    if inps.reference_file:
        print('reading reference info from reference: '+inps.reference_file)
        inps = read_reference_file2inps(inps.reference_file, inps)

    # Convert ref_lat/lon to ref_y/x
    coord = ut.coordinate(atr, lookup_file=inps.lookup_file)
    if inps.ref_lat and inps.ref_lon:
        (inps.ref_y,
         inps.ref_x) = coord.geo2radar(np.array(inps.ref_lat),
                                       np.array(inps.ref_lon))[0:2]
        print('input reference point in lat/lon: {}'.format((inps.ref_lat, inps.ref_lon)))

    # Check input ref_y/x
    if inps.ref_y is not None and inps.ref_x is not None:
        print('input reference point in y/x: {}'.format((inps.ref_y, inps.ref_x)))
        # Do not use ref_y/x outside of data coverage
        if not (0 <= inps.ref_y < length and 0 <= inps.ref_x < width):
            inps.ref_y, inps.ref_x = None, None
            raise ValueError('input reference point is OUT of data coverage!')

        # Do not use ref_y/x in masked out area
        if inps.maskFile and os.path.isfile(inps.maskFile):
            print('mask: '+inps.maskFile)
            mask = readfile.read(inps.maskFile, datasetName='mask')[0]
            if mask[inps.ref_y, inps.ref_x] == 0:
                inps.ref_y, inps.ref_x = None, None
                msg = 'input reference point is in masked OUT area defined by {}!'.format(inps.maskFile)
                raise ValueError(msg)

    else:
        # Determine auto selection method
        print('no input reference y/x.')
        if not inps.method:
            # Use existing REF_Y/X if 1) no ref_y/x input and 2) no method input and 3) ref_yx is in coverage
            if (not inps.force 
                    and 'REF_X' in atr.keys()
                    and 0 <= float(atr['REF_Y']) <= length 
                    and 0 <= float(atr['REF_X']) <= width):
                print('REF_Y/X exists in input file, skip updating.')
                print('REF_Y: '+atr['REF_Y'])
                print('REF_X: '+atr['REF_X'])
                inps.outfile = inps.file
                inps.go_reference = False

            # method to select reference point
            elif inps.coherenceFile and os.path.isfile(inps.coherenceFile):
                inps.method = 'maxCoherence'
            else:
                inps.method = 'random'
        print('reference point selection method: '+str(inps.method))
    print('-'*50)
    return inps


def remove_reference_pixel(File):
    """Remove reference pixel info from input file"""
    print("remove REF_Y/X and/or REF_LAT/LON from file: "+File)
    atrDrop = {}
    for i in ['REF_X', 'REF_Y', 'REF_LAT', 'REF_LON']:
        atrDrop[i] = 'None'
    File = ut.add_attribute(File, atrDrop)
    return File


#######################################  Main Function  ########################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    inps.file = ut.get_file_list(inps.file)[0]
    inps = read_reference_input(inps)

    if inps.go_reference:
        reference_file(inps)
    print('Done.')
    return


################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
