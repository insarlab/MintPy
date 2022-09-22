############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2016                               #
############################################################


import os
import time

import h5py
import numpy as np

from mintpy.objects import ifgramStack
from mintpy.objects.conncomp import connectComponent
from mintpy.utils import ptime, readfile, utils as ut, writefile

# key configuration parameter name
key_prefix = 'mintpy.unwrapError.'
config_keys = [
    'waterMaskFile',
    'connCompMinArea',
    'ramp',
    'bridgePtsRadius',
]


####################################################################################################
def run_or_skip(inps):
    print('-'*50)
    print('update mode: ON')
    flag = 'skip'

    # check output dataset
    with h5py.File(inps.ifgram_file, 'r') as f:
        if inps.datasetNameOut not in f.keys():
            flag = 'run'
            print(f'1) output dataset: {inps.datasetNameOut} NOT found.')
        else:
            print(f'1) output dataset: {inps.datasetNameOut} exists')
            ti = float(f[inps.datasetNameIn].attrs.get('MODIFICATION_TIME', os.path.getmtime(inps.ifgram_file)))
            to = float(f[inps.datasetNameOut].attrs.get('MODIFICATION_TIME', os.path.getmtime(inps.ifgram_file)))
            if ti > to:
                flag = 'run'
                print(f'2) output dataset is NOT newer than input dataset: {inps.datasetNameIn}.')
            else:
                print(f'2) output dataset is newer than input dataset: {inps.datasetNameIn}')

    # check configuration
    if flag == 'skip':
        # convert inps value to common str format
        inps_dict = dict(vars(inps))
        specialValues = {True : 'yes',
                         False: 'no',
                         None : 'no'}
        for key, value in inps_dict.items():
            if value in specialValues.keys():
                inps_dict[key] = specialValues[value]
        atr = readfile.read_attribute(inps.ifgram_file)

        # check all keys
        changed_keys = [key for key in config_keys
                        if str(inps_dict[key]) != atr.get(key_prefix+key, 'no')]
        if len(changed_keys) > 0:
            flag = 'run'
            print(f'3) NOT all key configuration parameters are the same: {config_keys}.')
            for key in changed_keys:
                print('\t{}\t: {} --> {}'.format(key, atr.get(key_prefix+key, 'no'), str(inps_dict[key])))
        else:
            print(f'3) all key configuration parameters are the same: {config_keys}.')

    # result
    print(f'run or skip: {flag}.')
    return flag


##########################################################################################
def run_unwrap_error_bridging(ifgram_file, water_mask_file, ramp_type=None, radius=50, cc_min_area=2.5e3,
                              ccName='connectComponent', dsNameIn='unwrapPhase',
                              dsNameOut='unwrapPhase_bridging', inps=None):
    """Run unwrapping error correction with bridging
    Parameters: ifgram_file     : str, path of ifgram stack file
                water_mask_file : str, path of water mask file
                ramp_type       : str, name of phase ramp to be removed during the phase jump estimation
                cc_min_area     : float, minimum region/area size
                ccName          : str, dataset name of connected components
                dsNameIn        : str, dataset name of unwrap phase to be corrected
                dsNameOut       : str, dataset name of unwrap phase to be saved after correction
                inps            : Namespace object, optional
    Returns:    ifgram_file     : str, path of ifgram stack file
    """
    start_time = time.time()
    print('-'*50)
    print(f'correct unwrapping error in {ifgram_file} with bridging ...')
    if ramp_type is not None:
        print(f'estimate and remove a {ramp_type} ramp while calculating phase offset')

    # read water mask
    if water_mask_file and os.path.isfile(water_mask_file):
        print('read water mask from file:', water_mask_file)
        water_mask = readfile.read(water_mask_file)[0]
    else:
        water_mask = None

    # file info
    fbase, fext = os.path.splitext(inps.ifgram_file)
    atr = readfile.read_attribute(ifgram_file)
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    k = atr['FILE_TYPE']

    # correct unwrap error ifgram by ifgram
    if k == 'ifgramStack':
        date12_list = ifgramStack(ifgram_file).get_date12_list(dropIfgram=False)
        date12_list_kept = ifgramStack(ifgram_file).get_date12_list(dropIfgram=True)
        num_ifgram = len(date12_list)
        shape_out = (num_ifgram, length, width)

        # prepare output data writing
        print(f'open {ifgram_file} with r+ mode')
        with h5py.File(ifgram_file, 'r+') as f:
            print('input  dataset:', dsNameIn)
            print('output dataset:', dsNameOut)
            if dsNameOut in f.keys():
                ds = f[dsNameOut]
                print(f'access /{dsNameOut} of np.float32 in size of {shape_out}')
            else:
                ds = f.create_dataset(
                    dsNameOut,
                    shape_out,
                    maxshape=(None, None, None),
                    chunks=True,
                    compression=None,
                )
                print(f'create /{dsNameOut} of np.float32 in size of {shape_out}')

            # correct unwrap error ifgram by ifgram
            prog_bar = ptime.progressBar(maxValue=num_ifgram)
            for i in range(num_ifgram):
                # read unwrapPhase
                date12 = date12_list[i]
                unw = np.squeeze(f[dsNameIn][i, :, :])

                # skip dropped interferograms
                if date12 not in date12_list_kept:
                    ds[i, :, :] = unw

                else:
                    # read connectComponent
                    cc = np.squeeze(f[ccName][i, :, :])
                    if water_mask is not None:
                        cc[water_mask == 0] = 0

                    # bridging
                    cc_obj = connectComponent(conncomp=cc, metadata=atr)
                    cc_obj.label(min_area=cc_min_area)
                    cc_obj.find_mst_bridge()
                    unw_cor = cc_obj.unwrap_conn_comp(unw, radius=radius, ramp_type=ramp_type)

                    # write to hdf5 file
                    ds[i, :, :] = unw_cor

                prog_bar.update(i+1, suffix=date12)
            prog_bar.close()
            ds.attrs['MODIFICATION_TIME'] = str(time.time())
        print(f'close {ifgram_file} file.')

    if k == '.unw':
        # read unwrap phase
        unw = readfile.read(ifgram_file)[0]

        # read connected components
        cc_files0 = [ifgram_file+'.conncomp', f'{fbase}_snap_connect.byt']
        cc_files = [i for i in cc_files0 if os.path.isfile(i)]
        if len(cc_files) == 0:
            raise FileNotFoundError(cc_files0)
        cc = readfile.read(cc_files[0])[0]
        if water_mask is not None:
            cc[water_mask == 0] = 0

        # bridging
        cc_obj = connectComponent(conncomp=cc, metadata=atr)
        cc_obj.label(min_area=cc_min_area)
        cc_obj.find_mst_bridge()
        unw_cor = cc_obj.unwrap_conn_comp(unw, ramp_type=ramp_type)

        # write to hdf5 file
        out_file = f'{fbase}_unwCor{fext}'
        print(f'writing >>> {out_file}')
        writefile.write(unw_cor, out_file=out_file, ref_file=ifgram_file)

    # update config parameter for HDF5 file
    if fext in ['.h5', '.he5'] and inps is not None:
        print('add/update the following configuration metadata to file:')
        config_meta = dict()
        for key in config_keys:
            config_meta[key_prefix+key] = str(vars(inps)[key])
        ut.add_attribute(inps.ifgram_file, config_meta, print_msg=True)

    # used time
    m, s = divmod(time.time() - start_time, 60)
    print(f'\ntime used: {m:02.0f} mins {s:02.1f} secs\nDone.')

    return ifgram_file
