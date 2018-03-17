#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Jan 2016: put manual matching code to manual_offset_estimate()
#                   put two files matching code into match_two_files()
#                   add cmdLineParse(), merge matching_all.py to it.


import os
import sys
import argparse

import h5py
import numpy as np
import matplotlib.pyplot as plt

import pysar.utils.readfile as readfile
import pysar.utils.writefile as writefile


#############################################################################################
def corners(atr):
    '''Get corners coordinate.'''
    width  = int(atr['WIDTH'])
    length = int(atr['FILE_LENGTH'])
    West  = float(atr['X_FIRST'])
    North = float(atr['Y_FIRST'])
    lon_step = float(atr['X_STEP'])
    lat_step = float(atr['Y_STEP'])
    South = North + lat_step*(length-1)
    East  = West + lon_step*(width-1)
    #lon_seq = np.arange(West, West +width *lon_step,lon_step)
    #lat_seq = np.arange(North,North+length*lat_step,lat_step)

    return West,East,North,South,width,length


#############################################################################################
def nearest(x, X):
    """ find nearest neighbour """
    dist = np.sqrt((X - x)**2)
    indx=np.where(dist==min(dist))
  
    return indx[0]


#############################################################################################
def manual_offset_estimate(matrix1, matrix2):
    '''Manually estimate offset between two data matrix.
    By manually selecting a line from each of them, and estimate the difference.
    It usually used when 2 input data matrix have no area in common.
    '''
    # Select line from data matrix 1
    fig = plt.figure()
    ax=fig.add_subplot(111)
    ax.imshow(matrix1)
    xc=[] 
    yc=[] 
    print('please click on start and end point of the desired profile/line')
    print('afterward close the figure to continue the process')
    def onclick(event):
        if event.button==1:
            print('click')
            xc.append(int(event.xdata))
            yc.append(int(event.ydata))
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    x0=xc[0];x1=xc[1]
    y0=yc[0];y1=yc[1]

    # Select line from data matrix 2
    fig = plt.figure()
    ax=fig.add_subplot(111)
    ax.imshow(matrix2)
    xc=[]
    yc=[]
    print('please click on start and end point of the desired profile')
    print('afterward close the figure to continue the process')
    def onclick(event):
        if event.button==1:
            print('click')
            xc.append(int(event.xdata))
            yc.append(int(event.ydata))
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    x00=xc[0];x11=xc[1]
    y00=yc[0];y11=yc[1]

    # Calculate the Offset - Difference
    #offset=V2[y00:y11,x00:x11]-V2[y0:y1,x0:x1]
    offset = np.nansum(V2[y00:y11,x00:x11]) / np.sum(np.isfinite(V2[y00:y11,x00:x11]))\
             - np.nansum(V1[y0:y1,x0:x1]) / np.sum(np.isfinite(V1[y0:y1,x0:x1]))

    return offset


#############################################################################################
def match_two_files(File1, File2, outName=None, manual_match=False, disp_fig=False):
    '''Match two geocoded files by estimating their offset.
    Better for two files with common area overlaping.
    '''
    
    # Read Input Files
    V1, atr1 = readfile.read(File1)
    V2, atr2 = readfile.read(File2)
    k = atr1['FILE_TYPE']
    print('---------------------------')
    print('matching 2 '+k+' files:\n'+File1+'\n'+File2)
    
    # Get Coverage Info 
    # Boundary Info - 2 Input Files
    West1,East1,North1,South1,width1,length1 = corners(atr1)
    West2,East2,North2,South2,width2,length2 = corners(atr2)
    # Boundary Info - Output File
    print('finding the corners of the whole area')
    West  = min(West1, West2)
    East  = max(East1, East2)
    North = max(North1,North2)
    South = min(South1,South2)
    lon_step = float(atr1['X_STEP'])
    lat_step = float(atr1['Y_STEP'])
    width  = int(round((East - West )/lon_step + 1.0))
    length = int(round((South - North)/lat_step + 1.0))

    # Get Index of Input Files in Output Files
    lon_seq = np.arange(West, West +width *lon_step, lon_step) 
    lat_seq = np.arange(North, North+length*lat_step, lat_step)
    indx1 = nearest(West1,  lon_seq)[0]
    indy1 = nearest(North1, lat_seq)[0]
    indx2 = nearest(West2,  lon_seq)[0]
    indy2 = nearest(North2, lat_seq)[0]

    # Estimate Offset of overlaping area
    VV1 = np.zeros([length,width])
    VV2 = np.zeros([length,width])
    VV1[:,:] = np.nan
    VV2[:,:] = np.nan
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

    # Adjust File2 value using offset
    if np.isnan(offset):
        print('**************************************************')
        print('WARNING:')
        print('')
        print('No offset is estimated and no matching applied.')
        print('Continue to merge two input files without any adjustment.')
        print('**************************************************')   
    else:
        print('Average offset between two velocity in the common area is: ' + str(offset))
        V2 = V2 - offset

    # Get merged data matrix value
    indv2 = np.isfinite(V2)
    VV = np.zeros([length,width])
    VV[:,:] = np.nan
    VV[indy1:indy1+length1, indx1:indx1+width1] = V1
    VV[indy2:indy2+length2, indx2:indx2+width2][indv2] = V2[indv2]
    
    # Write Output File
    if not outName:
        ext = os.path.splitext(File1)[1]
        outName = os.path.splitext(os.path.basename(File1))[0]+'_'+\
                  os.path.splitext(os.path.basename(File2))[0]+ext
    print('writing >>> '+outName)
    atr = atr1.copy()
    atr['WIDTH'] = width
    atr['FILE_LENGTH'] = length
    atr['X_FIRST'] = West
    atr['Y_FIRST'] = North
    writefile.write(VV, atr, outName)

    # Display
    fig_size = [16.0,16.0]
    fig = plt.figure(figsize=fig_size)
    print('plotting result ...')
    fig=plt.subplot(2,2,1);  plt.imshow(VV1);      plt.title(File1);     plt.colorbar()
    fig=plt.subplot(2,2,2);  plt.imshow(VV2);      plt.title(File2);     plt.colorbar()
    fig=plt.subplot(2,2,3);  plt.imshow(VV);       plt.title(outName);   plt.colorbar()
    fig=plt.subplot(2,2,4);  plt.imshow(VV_diff);  plt.title('Offset');  plt.colorbar()
    plt.tight_layout()
    plt.savefig(outName+'.png', bbox_inches='tight', transparent=True, dpi=150)
    print('save figure to '+outName+'.png')

    if disp_fig:
        print('showing ...')
        plt.show()

    return outName


#############################################################################################
EXAMPLE='''example:
  match.py  vel_AlosAT42*.h5
  match.py  vel_AlosAT42*.h5  -o vel_AlosA.h5
  match.py  vel_AlosAT422.h5  vel_AlosAT423.h5  vel_AlosAT424.h5  vel_AlosAT425.h5
  match.py  vel_AlosAT422.h5  vel_AlosAT423.h5
  match.py  vel_AlosAT422.h5  vel_AlosAT423.h5  --manual
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Match 2 or more geocoded datasets sharing common area.\n'
                                                 'Function automatically finds the common area and calculates\n'
                                                 'the average offset between the two velocity.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)
    parser.add_argument('first_file')
    parser.add_argument('other_file', nargs='+', help='file(s) to match')
    parser.add_argument('-o','--output', dest='outfile', help='output file name')
    parser.add_argument('--manual', dest='manual_match', action='store_true',\
                        help='manually select lines to estimate offset for matching.')
    parser.add_argument('--display', dest='disp_fig', action='store_true',\
                        help='do not display the result ploting.')

    inps = parser.parse_args()
    return inps


#############################################################################################
def main(argv):
    # Inputs
    inps = cmdLineParse()
    print('\n**************** Match Files *********************')
    print('Files to be matched:\n'+inps.first_file+', '+str(inps.other_file))
    
    # Matching file by file
    file1 = inps.first_file
    for file2 in inps.other_file:
        if file2 == inps.other_file[-1]:
            outName = inps.outfile
        else:
            outName = None
        file1 = match_two_files(file1, file2, outName, inps.manual_match, inps.disp_fig)
    print('Done.')
    return file1

#############################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])



