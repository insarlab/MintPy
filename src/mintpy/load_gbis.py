############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2019                               #
############################################################


import os
import warnings  # suppress UserWarning from matplotlib

warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

from mintpy.utils import writefile


##############################################################################
def grab_data_paths_from_inp_file(inp_file):
    """Grab data paths from inp file."""
    data_paths = []
    with open(inp_file) as f:
        lines = f.readlines()
        lines = [i.strip() for i in lines if not i.startswith('%')]
        lines = [i for i in lines if i]
        for line in lines:
            c = [i.strip() for i in line.strip().split('=', 1)]
            if (len(c) >= 2 and c[0] == 'insar{insarID}.dataPath'):
                mat_file = c[1].replace('\n','').split(";")[0].strip()
                mat_file = mat_file.replace("'","").replace('"','')
                data_paths.append(mat_file)
    return data_paths


def gbis_mat2hdf5(inv_mat_file, display=True):
    """Convert InSAR related GBIS inversion result .mat file into HDF5 file."""
    out_dir = os.path.dirname(inv_mat_file)
    out_files = []

    print(f'read mat file: {inv_mat_file}')
    mat = sio.loadmat(inv_mat_file, struct_as_record=False, squeeze_me=True)

    # when num_file == 1
    if isinstance(mat['insar'], sio.matlab.mio5_params.mat_struct):
        mat['insar'] = [mat['insar']]
        mat['insarPlot'] = [mat['insarPlot']]
    num_file = len(mat['insar'])
    print(f'number of output HDF5 file: {num_file}')

    # grab optimal model parameter
    parValue = mat['invResults'].optimalmodel
    parName = mat['invResults'].model.parName
    modelName = parName[0].split()[0]

    mDict = {}
    for i in range(len(parValue)):
        key = parName[i].replace(' ', '_')
        mDict[key] = parValue[i]
    mDict['DEFORMATION_MODEL'] = modelName

    if display:
        fig_size = [12, 3*num_file]
        fig, axs = plt.subplots(nrows=num_file, ncols=3, figsize=fig_size)
        axs = axs.reshape(-1,3)   #convert to 2D array when num_file is 1.
        print(f'creating figure in size of {fig_size}')

    for i in range(num_file):
        insar_mat_file = mat['insar'][i].dataPath
        if not os.path.isfile(insar_mat_file):
            inp_file = os.path.dirname(os.path.dirname(inv_mat_file))+'.inp'
            insar_mat_file = grab_data_paths_from_inp_file(inp_file)[i]
        print('-'*30)
        print(f'read mask from file: {insar_mat_file}')

        # read 1D height and 2D mask
        hgt1 = sio.loadmat(insar_mat_file, struct_as_record=False, squeeze_me=True)['Height']
        mask = sio.loadmat(insar_mat_file, struct_as_record=False, squeeze_me=True)['Mask']
        length, width = mask.shape

        # convert to 2D matrix
        insarPlot = mat['insarPlot'][i]
        out_file = os.path.join(out_dir, f'{insarPlot.name}.h5')
        out_files.append(out_file)

        hgt2 = np.zeros((length, width), dtype=np.float32) * np.nan
        data = np.zeros((length, width), dtype=np.float32) * np.nan
        model = np.zeros((length, width), dtype=np.float32) * np.nan
        residual = np.zeros((length, width), dtype=np.float32) * np.nan
        hgt2[mask!=0] = hgt1
        data[mask!=0] = insarPlot.data
        model[mask!=0] = insarPlot.model
        residual[mask!=0] = insarPlot.residual

        # prepare metadata
        meta = vars(sio.loadmat(insar_mat_file, struct_as_record=False, squeeze_me=True)['Metadata'])
        if '_fieldnames' in meta.keys():
            meta.pop('_fieldnames')
        meta['UNIT'] = 'm'
        meta['FILE_TYPE'] = 'displacement'
        meta['PROCESSOR'] = 'GBIS'
        for key, value in mDict.items():
            meta[key] = value

        # write to HDF5 file
        dsDict = {}
        dsDict['hgt'] = hgt2
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
