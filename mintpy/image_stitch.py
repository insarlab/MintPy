#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, Zhang Yunjun, 2013               #
############################################################
# Add --no-offset option by Robert Zinke, Nov 2020


import os
import sys
import numpy as np
import matplotlib.pyplot as plt
try:
    from skimage.transform import rescale
except ImportError:
    raise ImportError('Could not import skimage!')

from mintpy.utils import readfile, writefile, plot as pp
from mintpy.utils.arg_utils import create_argument_parser
from mintpy.multilook import multilook_data


#############################################################################################
EXAMPLE = """example:
  image_stitch.py  vel_AlosAT422.h5  vel_AlosAT423.h5  vel_AlosAT424.h5  vel_AlosAT425.h5 -o  vel_AlosA.h5
  image_stitch.py geom_AlosAT422.h5 geom_AlosAT423.h5 geom_AlosAT424.h5 geom_AlosAT425.h5 -o geom_AlosA.h5 --no-offset
"""

NOTE = """
  The function automatically:
  1) finds the common area between adjacent input files
  2) calculates the average offset between them
  3) apply this average offset to the later file
"""

def create_parser(subparsers=None):
    synopsis = 'Stitch/mosaic multiple geocoded datasets into one.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis+NOTE, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file1', help='file to stitch')
    parser.add_argument('file2s', nargs='+', metavar='file2', help='file(s) to stitch')
    parser.add_argument('-o', '--output', dest='outfile', required=True, help='output file name')

    # stitch option
    parser.add_argument('--no-offset','--no-off', dest='apply_offset', action='store_false',
                        help='Do not apply offset if 1) data sets are merely to be stitched '
                             'AND 2) no adjustment of values needs to be made\n'
                             '(i.e., for two coherence maps), use this flag')

    # plot options
    parser.add_argument('--nodisplay', dest='disp_fig', action='store_false',
                        help='do not display the result plotting.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


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
    vprint('estimate difference in the overlaping area')
    lon_seq = np.arange(W, W + width  * lon_step, lon_step)
    lat_seq = np.arange(N, N + length * lat_step, lat_step)
    x1, y1 = np.argmin(np.square(lon_seq - W1)), np.argmin(np.square(lat_seq - N1))
    x2, y2 = np.argmin(np.square(lon_seq - W2)), np.argmin(np.square(lat_seq - N2))

    # estimate offset of the overlaping area
    mat11 = np.zeros([length, width]) * np.nan;
    mat22 = np.zeros([length, width]) * np.nan;
    mat11[y1:y1+length1, x1:x1+width1] = mat1
    mat22[y2:y2+length2, x2:x2+width2] = mat2
    mat_diff = mat22 - mat11

    # apply the offset
    if apply_offset:
        offset = np.nansum(mat_diff) / np.sum(np.isfinite(mat_diff))
        if ~np.isnan(offset):
            vprint('average offset between two matrices in the common area: {}'.format(offset))
            vprint('offset all pixel values in the 2nd matrix by {} '.format(offset))
            mat2 -= offset
        else:
            print('*'*50)
            print('WARNING: NO common area found between two matrices!')
            print('    continue the stitching without applying offset')
            print('*'*50)

    # initiate output matrix
    # with the default value of NaN for float type and 0 for the other types
    vprint('create output metadata and matrix in shape of {}'.format((length, width)))
    fill_value = np.nan if str(mat1.dtype).startswith('float') else 0
    mat = np.zeros([length, width]) * fill_value

    # fill the output matrix
    flag2 = np.isfinite(mat2)
    mat[y1:y1+length1, x1:x1+width1] = mat1
    mat[y2:y2+length2, x2:x2+width2][flag2] = mat2[flag2]
    mat = np.array(mat, dtype=mat1.dtype)

    # output attributes
    atr = dict()
    for key, value in atr1.items():
        atr[key] = value
    atr['WIDTH'] = width
    atr['LENGTH'] = length
    atr['X_FIRST'] = W
    atr['Y_FIRST'] = N

    return mat, atr, mat11, mat22, mat_diff


def plot_stitch(mat11, mat22, mat, mat_diff, out_fig=None):
    """plot stitching result"""

    # plot settings
    titles = ['file 1', 'file 2', 'stitched', 'difference']
    mat_mli = multilook_data(mat, 20, 20, method='mean')
    vmin = np.nanmin(mat_mli)
    vmax = np.nanmax(mat_mli)
    fig_size = pp.auto_figure_size(ds_shape=mat.shape, scale=1.4, disp_cbar=True, print_msg=True)

    # plot
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=fig_size, sharex=True, sharey=True)
    for ax, data, title in zip(axs.flatten(), [mat11, mat22, mat, mat_diff], titles):
        im = ax.imshow(data, vmin=vmin, vmax=vmax, interpolation='nearest')
        ax.set_title(title, fontsize=12)
        ax.tick_params(which='both', direction='in', labelsize=12, left=True, right=True, top=True, bottom=True)
    fig.tight_layout()

    # colorbar
    fig.subplots_adjust(right=0.9)
    cax = fig.add_axes([0.901, 0.3, 0.01, 0.4])
    plt.colorbar(im, cax=cax)

    # output
    fig.savefig(out_fig, bbox_inches='tight', transparent=True, dpi=150)
    print('save figure to file: {}'.format(out_fig))

    #if disp_fig:
    #    print('showing ...')
    #    plt.show()
    #else:
    #    plt.close()

    return


def stitch_files(fnames, out_file, apply_offset=True, disp_fig=True, no_data_value=None):
    """Stitch all input files into one
    """
    fext = os.path.splitext(fnames[0])[1]
    atr = readfile.read_attribute(fnames[0])

    # grab ds_names
    ds_names = [None]
    if fext in ['.h5', '.he5']:
        # get the common dataset list among all input files
        ds_names = set(readfile.get_dataset_list(fnames[0]))
        for fname in fnames[1:]:
            ds_names.intersection_update(readfile.get_dataset_list(fname))
        ds_names = sorted(list(ds_names))

        # special treatment for velocity/time_function files
        if atr['FILE_TYPE'] == 'velocity' and len(ds_names) > 1:
            ds_names = ['velocity']

    print('files to be stitched: {}'.format(fnames))
    print('datasets to be stitched: {}'.format(ds_names))

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
            print('convert no_data_value from {} to NaN'.format(no_data_value))
            mat[mat==no_data_value] = np.nan

        # skip pixels with zero incidenceAngle for geometry files
        if atr['FILE_TYPE'] == 'geometry' and 'incidenceAngle' in ds_names:
            print('ignore pixels with ZERO incidenceAngle')
            inc_angle = readfile.read(fnames[0], datasetName='incidenceAngle')[0]
            mat[inc_angle == 0] = np.nan

        for i, fname in enumerate(fnames[1:]):
            print('-'*30)
            print('read data from file: {}'.format(fname))
            # reading
            mat2, atr2 = readfile.read(fname, datasetName=ds_name)
            # masking
            if no_data_value is not None:
                mat2[mat2==no_data_value] = np.nan
            # skip pixels with zero incidenceAngle for geometry files
            if atr['FILE_TYPE'] == 'geometry' and 'incidenceAngle' in ds_names:
                print('ignore pixels with ZERO incidenceAngle')
                inc_angle2 = readfile.read(fname, datasetName='incidenceAngle')[0]
                mat2[inc_angle2 == 0] = np.nan

            print('stitching ...')
            (mat, atr,
             mat11,
             mat22,
             mat_diff) = stitch_two_matrices(mat, atr,
                                             mat2, atr2,
                                             apply_offset=apply_offset)

            # plot
            if apply_offset:
                print('plot stitching & shifting result ...')
                suffix = '{}{}'.format(i, i+1)
                out_fig = '{}_{}.png'.format(os.path.splitext(out_file)[0], suffix)
                plot_stitch(mat11, mat22, mat, mat_diff, out_fig=out_fig)

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


#############################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    stitch_files(fnames=[inps.file1] + inps.file2s,
                 out_file=inps.outfile,
                 apply_offset=inps.apply_offset,
                 disp_fig=inps.disp_fig)

    return


#############################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
