############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, Zhang Yunjun, 2013               #
############################################################
# Add --no-offset option by Robert Zinke, Nov 2020


import os

import matplotlib.pyplot as plt
import numpy as np
from skimage.transform import rescale

from mintpy.multilook import multilook_data
from mintpy.utils import plot as pp, readfile, writefile


#############################################################################################
def manual_offset_estimate(mat1, mat2):
    """Manually estimate offset between two data matrix.
    By manually selecting a line from each of them, and estimate the difference.
    It usually used when 2 input data matrix have no area in common.
    """
    def onclick(event):
        if event.button == 1:
            print('click')
            xc.append(int(event.xdata))
            yc.append(int(event.ydata))

    # Select line from data matrix 1
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(mat1)
    xc = []
    yc = []
    print('please click on start and end point of the desired profile/line')
    print('afterward close the figure to continue the process')
    fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    x0 = xc[0]
    x1 = xc[1]
    y0 = yc[0]
    y1 = yc[1]

    # Select line from data matrix 2
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(mat2)
    xc = []
    yc = []
    print('please click on start and end point of the desired profile')
    print('afterward close the figure to continue the process')
    fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    x00 = xc[0]
    x11 = xc[1]
    y00 = yc[0]
    y11 = yc[1]

    # Calculate the Offset - Difference
    # offset=V2[y00:y11,x00:x11]-V2[y0:y1,x0:x1]
    offset = (  np.nansum(mat2[y00:y11, x00:x11]) / np.sum(np.isfinite(mat2[y00:y11, x00:x11]))
              - np.nansum(mat1[y0:y1, x0:x1]) / np.sum(np.isfinite(mat1[y0:y1, x0:x1])))

    return offset


#############################################################################################
def rescale_data(data, meta, ref_meta):
    """rescale matrix into a different resolution"""

    # calc scaling factor
    scale = (float(meta['Y_STEP']) / float(ref_meta['Y_STEP']),
             float(meta['X_STEP']) / float(ref_meta['X_STEP']))
    # scale
    data_out = rescale(data, scale)
    # update metadata
    meta['Y_STEP'] = ref_meta['Y_STEP']
    meta['X_STEP'] = ref_meta['X_STEP']
    meta['LENGTH'], meta['WIDTH'] = data_out.shape

    return data_out, meta


def get_corners(atr):
    """Get corners coordinate."""
    length = int(atr['LENGTH'])
    width = int(atr['WIDTH'])
    W = float(atr['X_FIRST'])
    N = float(atr['Y_FIRST'])
    lon_step = float(atr['X_STEP'])
    lat_step = float(atr['Y_STEP'])
    S = N + lat_step * length
    E = W + lon_step * width

    return S, N, W, E, width, length


def stitch_two_matrices(mat1, atr1, mat2, atr2, apply_offset=True, print_msg=True):
    """Stitch two geocoded matrices
    and/or shift the 2nd matrix value to match with the 1st one.

    Parameters: mat1/2       - 2D np.ndarray
                atr1/2       - dict, attributes
                apply_offset - bool, estimate offset and adjust 2nd matrix value
    Returns:    mat          - 2D np.ndarray, stitched matrix
                atr          - dict, attributes of stitched matrix
    """
    vprint = print if print_msg else lambda *args, **kwargs: None

    # resize the 2nd matrix, if it has different spatial resolution
    ratio_x = abs((float(atr1['X_STEP']) - float(atr2['X_STEP'])) / float(atr1['X_STEP']))
    ratio_y = abs((float(atr1['Y_STEP']) - float(atr2['Y_STEP'])) / float(atr1['Y_STEP']))
    if any(i > 1e-3 for i in [ratio_x, ratio_y]):
        vprint('file 1: X_STEP - {}, Y_STEP - {}'.format(atr1['X_STEP'], atr1['Y_STEP']))
        vprint('file 2: X_STEP - {}, Y_STEP - {}'.format(atr2['X_STEP'], atr2['Y_STEP']))
        vprint('rescale the 2nd matrix into the same spatial resolution as the 1st one ...')
        mat2, atr2 = rescale_data(mat2, meta=atr2, ref_meta=atr1)

    # input spatial extents
    vprint('grab corners of input matrices')
    S1, N1, W1, E1, width1, length1 = get_corners(atr1)
    S2, N2, W2, E2, width2, length2 = get_corners(atr2)

    # output spatial extent
    vprint('calculate corners of output matrix')
    W, E = min(W1, W2), max(E1, E2)
    S, N = min(S1, S2), max(N1, N2)
    lon_step = float(atr1['X_STEP'])
    lat_step = float(atr1['Y_STEP'])
    width  = int(np.ceil((E - W) / lon_step))
    length = int(np.ceil((S - N) / lat_step))

    # index of input matrices in output matrix
    vprint('estimate difference in the overlapping area')
    lon_seq = np.arange(W, W + width  * lon_step, lon_step)
    lat_seq = np.arange(N, N + length * lat_step, lat_step)
    x1, y1 = np.argmin(np.square(lon_seq - W1)), np.argmin(np.square(lat_seq - N1))
    x2, y2 = np.argmin(np.square(lon_seq - W2)), np.argmin(np.square(lat_seq - N2))

    # estimate offset of the overlapping area
    mat11 = np.zeros([length, width]) * np.nan;
    mat22 = np.zeros([length, width]) * np.nan;
    mat11[y1:y1+length1, x1:x1+width1] = mat1
    mat22[y2:y2+length2, x2:x2+width2] = mat2
    mat_diff = mat22 - mat11

    # apply the offset
    if apply_offset:
        offset = np.nansum(mat_diff) / np.sum(np.isfinite(mat_diff))
        if ~np.isnan(offset):
            vprint(f'average offset between two matrices in the common area: {offset}')
            vprint(f'offset all pixel values in the 2nd matrix by {offset} ')
            mat2 -= offset
        else:
            print('*'*50)
            print('WARNING: NO common area found between two matrices!')
            print('    continue the stitching without applying offset')
            print('*'*50)

    # initiate output matrix
    # with the default value of NaN for float type and 0 for the other types
    vprint(f'create output metadata and matrix in shape of {(length, width)}')
    fill_value = np.nan if str(mat1.dtype).startswith('float') else 0
    mat = np.zeros([length, width], dtype=mat1.dtype) * fill_value

    # fill the output matrix
    flag2 = np.isfinite(mat2)
    mat[y1:y1+length1, x1:x1+width1] = mat1
    mat[y2:y2+length2, x2:x2+width2][flag2] = mat2[flag2]

    # output attributes
    atr = dict()
    for key, value in atr1.items():
        atr[key] = value
    atr['WIDTH'] = width
    atr['LENGTH'] = length
    atr['X_FIRST'] = W
    atr['Y_FIRST'] = N

    return mat, atr, mat11, mat22, mat_diff


def plot_stitch(mat11, mat22, mat, mat_diff, out_fig=None, disp_scale=1, disp_vlim=None, disp_cmap=None):
    """plot stitching result"""

    # plot settings
    titles = ['file 1', 'file 2', 'stitched', 'difference']
    if disp_scale != 1:
        print(f'scale the data by a factor of {disp_scale} for plotting')
    mat_mli = multilook_data(mat, 20, 20, method='mean')
    vmin = disp_vlim[0] if disp_vlim else np.nanmin(mat_mli) * disp_scale
    vmax = disp_vlim[1] if disp_vlim else np.nanmax(mat_mli) * disp_scale

    fig_size = pp.auto_figure_size(ds_shape=mat.shape, scale=1.4, disp_cbar=True, print_msg=True)

    # plot
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=fig_size, sharex=True, sharey=True)
    for ax, data, title in zip(axs.flatten(), [mat11, mat22, mat, mat_diff], titles):
        im = ax.imshow(data * disp_scale, vmin=vmin, vmax=vmax, cmap=disp_cmap, interpolation='nearest')
        ax.set_title(title, fontsize=12)
        ax.tick_params(which='both', direction='in', labelsize=12, left=True, right=True, top=True, bottom=True)
    fig.tight_layout()

    # colorbar
    fig.subplots_adjust(right=0.9)
    cax = fig.add_axes([0.901, 0.3, 0.01, 0.4])
    plt.colorbar(im, cax=cax)

    # output
    fig.savefig(out_fig, bbox_inches='tight', transparent=True, dpi=150)
    print(f'save figure to file: {out_fig}')

    return


def stitch_files(fnames, out_file, apply_offset=True, disp_fig=True, no_data_value=None,
                 disp_scale=1, disp_vlim=None, disp_cmap=None):
    """Stitch all input files into one
    """
    fext = os.path.splitext(fnames[0])[1]
    atr = readfile.read_attribute(fnames[0])

    # grab ds_names
    ds_names = set(readfile.get_dataset_list(fnames[0]))
    # get the common dataset list among all input files
    for fname in fnames[1:]:
        ds_names.intersection_update(readfile.get_dataset_list(fname))
    ds_names = sorted(list(ds_names))

    # special treatment for velocity/time_function files
    if atr['FILE_TYPE'] == 'velocity' and len(ds_names) > 1:
        ds_names = ['velocity']

    print(f'files to be stitched: {fnames}')
    print(f'datasets to be stitched: {ds_names}')

    # stitching
    dsDict = {}
    for ds_name in ds_names:
        # reading
        mat, atr = readfile.read(fnames[0], datasetName=ds_name)
        ds_name_out = ds_name if ds_name else atr['FILE_TYPE']
        print('#'*50)
        print(f'read {ds_name_out} from file: {fnames[0]}')

        # masking
        if no_data_value is not None:
            print(f'convert no_data_value from {no_data_value} to NaN')
            mat[mat==no_data_value] = np.nan

        # skip pixels with zero incidenceAngle for geometry files
        if atr['FILE_TYPE'] in ['geometry', 'los'] and 'incidenceAngle' in ds_names:
            print('ignore pixels with ZERO incidenceAngle')
            inc_angle = readfile.read(fnames[0], datasetName='incidenceAngle')[0]
            mat[inc_angle == 0] = np.nan

        for i, fname in enumerate(fnames[1:]):
            print('-'*30)
            print(f'read data from file: {fname}')
            # reading
            mat2, atr2 = readfile.read(fname, datasetName=ds_name)
            # masking
            if no_data_value is not None:
                mat2[mat2==no_data_value] = np.nan
            # skip pixels with zero incidenceAngle for geometry files
            if atr['FILE_TYPE'] in ['geometry', 'los'] and 'incidenceAngle' in ds_names:
                print('ignore pixels with ZERO incidenceAngle')
                inc_angle2 = readfile.read(fname, datasetName='incidenceAngle')[0]
                mat2[inc_angle2 == 0] = np.nan

            print('stitching ...')
            mat, atr, mat11, mat22, mat_diff = stitch_two_matrices(
                mat, atr,
                mat2, atr2,
                apply_offset=apply_offset)

            # plot
            if apply_offset:
                print('plot stitching & shifting result ...')
                out_fig = f'{os.path.splitext(out_file)[0]}_{i}{i+1}.png'
                plot_stitch(
                    mat11, mat22,
                    mat, mat_diff,
                    out_fig=out_fig,
                    disp_scale=disp_scale,
                    disp_vlim=disp_vlim,
                    disp_cmap=disp_cmap,
                )

        dsDict[ds_name_out] = mat

    # write output file
    print('#'*50)
    writefile.write(dsDict, out_file=out_file, metadata=atr)

    # plot
    if disp_fig:
        print('showing ...')
        plt.show()
    else:
        plt.close()

    return out_file
