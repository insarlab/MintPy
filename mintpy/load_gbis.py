#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2019                               #
############################################################


import os
import argparse
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
# suppress UserWarning from matplotlib
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")

from mintpy.utils import writefile


EXAMPLE = """example:
  load_gbis.py invert_1_2_C.mat
  load_gbis.py invert_1_2_C.mat --nodisplay
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Load GBIS inversion result to HDF5 format.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='GBIS inversion mat file.')
    parser.add_argument('-o', '--output', dest='outfile', help='output file name.')
    parser.add_argument('--nodisplay', dest='disp_fig', action='store_false', help='do not display the figure')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    inps.file = os.path.abspath(inps.file)

    # Backend setting
    if not inps.disp_fig:
            plt.switch_backend('Agg')

    return inps


def gbis_mat2hdf5(inv_mat_file, display=True):
    """Convert InSAR related GBIS inversion result .mat file into HDF5 file."""
    out_dir = os.path.dirname(inv_mat_file)
    out_files = []

    print('read mat file: {}'.format(inv_mat_file))
    mat = sio.loadmat(inv_mat_file, struct_as_record=False, squeeze_me=True)

    # when num_file == 1
    if isinstance(mat['insar'], sio.matlab.mio5_params.mat_struct):
        mat['insar'] = [mat['insar']]
        mat['insarPlot'] = [mat['insarPlot']]
    num_file = len(mat['insar'])
    print('number of output HDF5 file: {}'.format(num_file))

    # grab optimal model parameter
    parValue = mat['invResults'].optimalmodel
    parName = mat['invResults'].model.parName
    modelName = parName[0].split()[0]

    mDict = {}
    for i in range(len(parValue)):
        key = parName[i].replace(' ', '_')
        mDict[key] = parValue[i]

    # convert model x/y to lat/lon
    #ref_lon, ref_lat = mat['geo'].referencePoint
    #mDict['{}_REF_LAT'.format(modelName)] = ref_lat
    #mDict['{}_REF_LON'.format(modelName)] = ref_lon

    #import pyproj
    #geod = pyproj.Geod(ellps='WGS84')
    #mX, mY = parValue[0:2]
    #mLon, mLat = geod.fwd(ref_lon, ref_lat,
    #                      az=np.arctan(mX/mY) * 180 / np.pi,
    #                      dist=(mX**2 + mY**2)**0.5,
    #                      radians=False)[0:2]
    #mDict['{}_LAT'.format(modelName)] = mLat
    #mDict['{}_LON'.format(modelName)] = mLon
    mDict['DEFORMATION_MODEL'] = modelName

    if display:
        fig_size = [12, 3*num_file]
        fig, axs = plt.subplots(nrows=num_file, ncols=3, figsize=fig_size)
        axs = axs.reshape(-1,3)   #convert to 2D array when num_file is 1.
        print('creating figure in size of {}'.format(fig_size))

    for i in range(num_file):
        insar_mat_file = mat['insar'][i].dataPath
        print('-'*30)
        print('read mask from file: {}'.format(insar_mat_file))

        # read mask
        mask = sio.loadmat(insar_mat_file, struct_as_record=False, squeeze_me=True)['Mask']
        length, width = mask.shape

        # convert to 2D matrix
        insarPlot = mat['insarPlot'][i]
        out_file = os.path.join(out_dir, '{}.h5'.format(insarPlot.name))
        out_files.append(out_file)

        data = np.zeros((length, width), dtype=np.float32) * np.nan
        model = np.zeros((length, width), dtype=np.float32) * np.nan
        residual = np.zeros((length, width), dtype=np.float32) * np.nan
        data[mask!=0] = insarPlot.data
        model[mask!=0] = insarPlot.model
        residual[mask!=0] = insarPlot.residual

        # prepare metadata
        meta = vars(sio.loadmat(insar_mat_file, struct_as_record=False, squeeze_me=True)['Metadata'])
        temp = meta.pop('_fieldnames') # remote _fieldnames added by Matlab
        meta['UNIT'] = 'm'
        meta['FILE_TYPE'] = 'displacement'
        meta['PROCESSOR'] = 'GBIS'
        for key, value in mDict.items():
            meta[key] = value

        # write to HDF5 file
        dsDict = {}
        dsDict['data'] = data
        dsDict['model'] = model
        dsDict['residual'] = residual
        writefile.write(dsDict, out_file=out_file, metadata=meta)

        # plot
        if display:
            dlim = np.nanmax(np.abs(data))
            for ax, d in zip(axs[i, :], [data, model, residual]):
                im = ax.imshow(d, vmin=-dlim, vmax=dlim, cmap='jet')
                fig.colorbar(im, ax=ax)
            axs[i,0].set_ylabel('\n'.join(insarPlot.name.split('_', 1)))

    if display:
        print('showing ...')
        axs[0,0].set_title('data')
        axs[0,1].set_title('model')
        axs[0,2].set_title('residual')
        plt.show()
    return out_files


##############################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    out_files = gbis_mat2hdf5(inps.file, display=inps.disp_fig)

    return out_files


##########################################################################
if __name__ == '__main__':
    main()
