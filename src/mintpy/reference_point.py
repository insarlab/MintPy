############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os
import random

import h5py
import numpy as np

from mintpy.utils import ptime, readfile, utils as ut, writefile


###############################################################
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
    if atr['FILE_TYPE'] == 'ifgramStack':
        ds_name = [i for i in readfile.get_dataset_list(inps.file)
                   if i in ['unwrapPhase', 'rangeOffset', 'azimuthOffset']][0]
    else:
        ds_name = None
    stack = ut.temporal_average(inps.file, datasetName=ds_name, updateMode=True, outFile=False)[0]
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
            inps.ref_y, inps.ref_x = select_max_coherence_yx(
                coh_file=inps.coherenceFile,
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
        inps.outfile = inps.outfile if inps.outfile else inps.file
        ftype = atr['FILE_TYPE']
        fext = os.path.splitext(inps.file)[1]

        if fext == '.h5':
            if inps.outfile == inps.file:
                print('updating dataset values without re-writing to a new file')

                if ftype == 'ifgramStack':
                    with h5py.File(inps.file, 'r+') as f:
                        ds = f['unwrapPhase']
                        num_date12 = ds.shape[0]
                        prog_bar = ptime.progressBar(maxValue=num_date12)
                        for i in range(num_date12):
                            prog_bar.update(i+1, suffix=f'{i+1} / {num_date12}')

                            # make a copy of ds[i] because h5py allows fancy indexing for 1D arrays only.
                            data_2d = ds[i, :, :]

                            # apply spatial referencing (skip pixels with no-data-value)
                            data_2d[data_2d != 0.] -= data_2d[inps.ref_y, inps.ref_x]
                            ds[i, :, :] = data_2d

                        prog_bar.close()

                        print('update metadata')
                        f.attrs.update(atrNew)

                else:
                    with h5py.File(inps.file, 'r+') as f:
                        ds = f[ftype]
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
                print(f'writing the referenced data into file: {inps.outfile}')

                # 1. read and update data value
                data, atr = readfile.read(inps.file, datasetName=ftype)
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
            dis_names = ['phase', 'displacement']
            ds_names = readfile.get_dataset_list(inps.file)
            ds_dict = {}
            for ds_name in ds_names:
                data = readfile.read(inps.file, datasetName=ds_name)[0]
                if ds_name in dis_names:
                    data -= data[inps.ref_y, inps.ref_x]
                else:
                    print(f"skip spatial referencing for {ds_name}, as it's not in {dis_names}")
                ds_dict[ds_name] = data
            atr.update(atrNew)
            writefile.write(ds_dict, out_file=inps.outfile, metadata=atr)

    ut.touch([inps.coherenceFile, inps.maskFile])
    return inps.outfile


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
    """Manually select reference point in row/column number.
    Parameters: data : 2D np.ndarray, stack of input file
                inps : namespace, with key 'REF_X' and 'REF_Y', which will be updated
                mask : 2D np.ndarray
    """
    from matplotlib import pyplot as plt
    print('\nManual select reference point ...')
    print('Click on a pixel that you want to choose as the reference ')
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
                print('The selected pixel has NaN value in data.')
                print('Try a difference location please.')

    fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    print(f'y/x: {(inps.ref_y, inps.ref_x)}')

    return inps


def select_max_coherence_yx(coh_file, mask=None, min_coh=0.85):
    """Select pixel with coherence > min_coh in random"""
    print(f'random select pixel with coherence > {min_coh}')
    print('\tbased on coherence file: '+coh_file)
    coh = readfile.read(coh_file)[0]
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
    print(f'y/x: {(y, x)}')
    return y, x


def random_select_reference_yx(data_mat, print_msg=True):
    nrow, ncol = np.shape(data_mat)
    y = random.choice(list(range(nrow)))
    x = random.choice(list(range(ncol)))
    while data_mat[y, x] == 0:
        y = random.choice(list(range(nrow)))
        x = random.choice(list(range(ncol)))
    if print_msg:
        print(f'random select pixel\ny/x: {(y, x)}')
    return y, x


###############################################################
def read_template2inps(template_file, inps):
    """Read seed/reference info from template file and update input namespace"""
    inps_dict = vars(inps)
    template = readfile.read_template(template_file, skip_chars=['[', ']'])
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
            inps.ref_y, inps.ref_x = (int(i) for i in value.split(','))

    key = prefix+'lalo'
    if key in template.keys():
        value = template[key]
        if value:
            inps.ref_lat, inps.ref_lon = (float(i) for i in value.split(','))

    return inps


def read_reference_file2inps(reference_file, inps=None):
    """Read reference info from reference file and update input namespace"""
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
        inps = read_template2inps(inps.template_file, inps)

    if inps.reference_file:
        print('reading reference info from reference: '+inps.reference_file)
        inps = read_reference_file2inps(inps.reference_file, inps)

    if inps.ref_lat and np.abs(inps.ref_lat) > 90 and 'UTM_ZONE' not in atr.keys():
        msg = f'input reference latitude ({inps.ref_lat}) > 90 deg in magnitude!'
        msg += ' This does not make sense, double check your inputs!'
        raise ValueError(msg)

    # Convert ref_lat/lon to ref_y/x
    coord = ut.coordinate(atr, lookup_file=inps.lookup_file)
    if inps.ref_lat and inps.ref_lon:
        (inps.ref_y,
         inps.ref_x) = coord.geo2radar(np.array(inps.ref_lat),
                                       np.array(inps.ref_lon))[0:2]
        print(f'input reference point in lat/lon: {(inps.ref_lat, inps.ref_lon)}')

    # Check input ref_y/x
    if inps.ref_y is not None and inps.ref_x is not None:
        print(f'input reference point in y/x: {(inps.ref_y, inps.ref_x)}')
        # Do not use ref_y/x outside of data coverage
        if not (0 <= inps.ref_y < length and 0 <= inps.ref_x < width):
            inps.ref_y, inps.ref_x = None, None
            raise ValueError('input reference point is OUT of data coverage!')

        # Do not use ref_y/x in masked out area
        if inps.maskFile and os.path.isfile(inps.maskFile):
            print('mask: '+inps.maskFile)
            ds_names = readfile.get_dataset_list(inps.maskFile)
            ds_name = [x for x in ds_names if x in ['mask', 'waterMask']][0]
            mask = readfile.read(inps.maskFile, datasetName=ds_name)[0]
            if mask[inps.ref_y, inps.ref_x] == 0:
                inps.ref_y, inps.ref_x = None, None
                msg = f'input reference point is in masked OUT area defined by {inps.maskFile}!'
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
