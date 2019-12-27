#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, Zhang Yunjun, 2013               #
############################################################


import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
try:
    from skimage.transform import rescale
except ImportError:
    raise ImportError('Could not import skimage!')
from mintpy.utils import readfile, writefile


#############################################################################################
EXAMPLE = """example:
  match.py  vel_AlosAT42*.h5
  match.py  vel_AlosAT42*.h5  -o vel_AlosA.h5
  match.py  vel_AlosAT422.h5  vel_AlosAT423.h5  vel_AlosAT424.h5  vel_AlosAT425.h5
  match.py  vel_AlosAT422.h5  vel_AlosAT423.h5
  match.py  vel_AlosAT422.h5  vel_AlosAT423.h5  --manual
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Merge 2 or more geocoded datasets sharing common area into one.\n'
                                                 '\tFunction automatically finds the common area and calculates\n'
                                                 '\tthe average offset between the two velocity.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)
    parser.add_argument('file1')
    parser.add_argument('file2s', nargs='+', help='file(s) to match')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output file name')

    parser.add_argument('--manual', dest='manual_match', action='store_true',
                        help='manually select lines to estimate offset for matching.')
    parser.add_argument('--nodisplay', dest='disp_fig', action='store_false',
                        help='do not display the result ploting.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


#############################################################################################
def corners(atr):
    """Get corners coordinate."""
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    W = float(atr['X_FIRST'])
    N = float(atr['Y_FIRST'])
    lon_step = float(atr['X_STEP'])
    lat_step = float(atr['Y_STEP'])
    S = N + lat_step * (length - 1)
    E = W + lon_step * (width - 1)
    #lon_seq = np.arange(West, West +width *lon_step,lon_step)
    #lat_seq = np.arange(North,North+length*lat_step,lat_step)
    return S, N, W, E, width, length


#############################################################################################
def nearest(x, x_seq):
    """ find nearest neighbour """
    dist = np.sqrt((x_seq - x)**2)
    indx = np.where(dist == min(dist))
    return indx[0]


#############################################################################################
def manual_offset_estimate(matrix1, matrix2):
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
    ax.imshow(matrix1)
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
    ax.imshow(matrix2)
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
    offset = np.nansum(V2[y00:y11, x00:x11]) / np.sum(np.isfinite(V2[y00:y11, x00:x11])) - np.nansum(V1[y0:y1, x0:x1]) / np.sum(np.isfinite(V1[y0:y1, x0:x1]))

    return offset


def rescale_data_resolution(data, metadata, ref_metadata):
    scale = (float(metadata['Y_STEP']) / float(ref_metadata['Y_STEP']),
             float(metadata['X_STEP']) / float(ref_metadata['X_STEP']))
    data_out = rescale(data, scale)
    metadata['Y_STEP'] = ref_metadata['Y_STEP']
    metadata['X_STEP'] = ref_metadata['X_STEP']
    metadata['LENGTH'], metadata['WIDTH'] = data_out.shape
    return data_out, metadata


#############################################################################################
def match_two_files(file1, file2, out_file=None, manual_match=False, disp_fig=False):
    """Match two geocoded files by estimating their offset.
    Better for two files with common area overlaping.
    """
    # Read Input Files
    V1, atr1 = readfile.read(file1)
    V2, atr2 = readfile.read(file2)

    k = atr1['FILE_TYPE']
    print('---------------------------')
    print('matching 2 {} files: {} and {}'.format(k, file1, file2))

    # if two files are in different spatial resolution
    if any(atr1[key] != atr2[key] for key in ['X_STEP', 'Y_STEP']):
        print('resizing file2 into the same spatial resolution as file1 ...')
        V2, atr2 = rescale_data_resolution(V2, metadata=atr2, ref_metadata=atr1)

    # Get Coverage Info
    # Boundary Info - 2 Input Files
    S1, N1, W1, E1, width1, length1 = corners(atr1)
    S2, N2, W2, E2, width2, length2 = corners(atr2)

    # Boundary Info - Output File
    print('finding the corners of the whole area')
    W, E = min(W1, W2), max(E1, E2)
    S, N = min(S1, S2), max(N1, N2)
    lon_step = float(atr1['X_STEP'])
    lat_step = float(atr1['Y_STEP'])
    width = int(round((E - W) / lon_step + 1.0))
    length = int(round((S - N) / lat_step + 1.0))

    # Get Index of Input Files in Output Files
    lon_seq = np.arange(W, W + width * lon_step, lon_step)
    lat_seq = np.arange(N, N + length * lat_step, lat_step)
    indx1, indy1 = np.argmin(np.square(lon_seq - W1)), np.argmin(np.square(lat_seq - N1))
    indx2, indy2 = np.argmin(np.square(lon_seq - W2)), np.argmin(np.square(lat_seq - N2))

    # Estimate Offset of overlaping area
    VV1 = np.zeros([length, width])
    VV2 = np.zeros([length, width])
    VV1[:, :] = np.nan
    VV2[:, :] = np.nan
    VV1[indy1:indy1+length1, indx1:indx1+width1] = V1
    VV2[indy2:indy2+length2, indx2:indx2+width2] = V2

    if not manual_match:
        VV_diff = VV2 - VV1
        offset = np.nansum(VV_diff) / np.sum(np.isfinite(VV_diff))

    if np.isnan(offset):
        print('**************************************************')
        print('WARNING:')
        print('')
        print('No common area found between two velocity maps')
        print('At least one common pixel is required.')
        print('No matching applied. ')
        print('Continue with manual matching ...')
        print('    by selecting two line from each dataset to calculate the offset')
        print('**************************************************')
        manual_matching = True

    if manual_match:
        offset = manual_offset_estimate(V1, V2)

    # Adjust file2 value using offset
    if np.isnan(offset):
        print('**************************************************')
        print('WARNING:')
        print('')
        print('No offset is estimated and no matching applied.')
        print('continue to merge two input files without any adjustment.')
        print('**************************************************')
    else:
        print('average offset between two velocity in the common area is: ' + str(offset))
        V2 = V2 - offset

    # Get merged data matrix value
    indv2 = np.isfinite(V2)
    VV = np.zeros([length, width])
    VV[:, :] = np.nan
    VV[indy1:indy1+length1, indx1:indx1+width1] = V1
    VV[indy2:indy2+length2, indx2:indx2+width2][indv2] = V2[indv2]

    # Write Output File
    if not out_file:
        out_file = '{}_{}{}'.format(os.path.splitext(os.path.basename(file1))[0],
                                    os.path.splitext(os.path.basename(file2))[0],
                                    os.path.splitext(file1)[1])
    atr = atr1.copy()
    atr['WIDTH'] = width
    atr['LENGTH'] = length
    atr['X_FIRST'] = W
    atr['Y_FIRST'] = N
    writefile.write(VV, out_file=out_file, metadata=atr)

    # Display
    fig = plt.figure(figsize=[15.0, 8.0])
    print('plotting result ...')
    fig=plt.subplot(2,2,1);  plt.imshow(VV1);      plt.title('File1');   plt.colorbar()
    fig=plt.subplot(2,2,2);  plt.imshow(VV2);      plt.title('File2');   plt.colorbar()
    fig=plt.subplot(2,2,3);  plt.imshow(VV);       plt.title('Merged');  plt.colorbar()
    fig=plt.subplot(2,2,4);  plt.imshow(VV_diff);  plt.title('Offset');  plt.colorbar()

    plt.tight_layout()
    plt.savefig(out_file+'.png', bbox_inches='tight', transparent=True, dpi=150)
    print('save figure to '+out_file+'.png')

    if disp_fig:
        print('showing ...')
        plt.show()

    return out_file


#############################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    print('\n**************** Match Files *********************')
    print('Files to be matched:\n{}, {}'.format(inps.file1, inps.file2s))

    # Loop to match two files at a time
    file1 = inps.file1
    for file2 in inps.file2s:
        # use inps.outfile for the last match
        if file2 == inps.file2s[-1]:
            out_file = inps.outfile
        else:
            out_file = None

        # match two files
        file1 = match_two_files(file1, file2, out_file, inps.manual_match, inps.disp_fig)

    print('Done.')
    return file1


#############################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
