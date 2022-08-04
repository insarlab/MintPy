#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2016                               #
############################################################


import os
import sys
import time
import h5py
import numpy as np

from mintpy.objects import ifgramStack
from mintpy.objects.conncomp import connectComponent
from mintpy.defaults.template import get_template_content
from mintpy.utils import ptime, readfile, writefile, utils as ut
from mintpy.utils.arg_utils import create_argument_parser


# key configuration parameter name
key_prefix = 'mintpy.unwrapError.'
configKeys = [
    'waterMaskFile',
    'connCompMinArea',
    'ramp',
    'bridgePtsRadius',
]


####################################################################################################
TEMPLATE = get_template_content('correct_unwrap_error')

REFERENCE = """reference:
  Yunjun, Z., H. Fattahi, and F. Amelung (2019), Small baseline InSAR time series analysis:
  Unwrapping error correction and noise reduction, Computers & Geosciences, 133, 104331,
  doi:10.1016/j.cageo.2019.104331.
"""

EXAMPLE = """Example:
  unwrap_error_bridging.py  ./inputs/ifgramStack.h5  -t GalapagosSenDT128.template --update
  unwrap_error_bridging.py  ./inputs/ifgramStack.h5  --water-mask waterMask.h5
  unwrap_error_bridging.py  20180502_20180619.unw    --water-mask waterMask.h5
"""

NOTE = """
  by connecting reliable regions with MST bridges. This method assumes the phase differences
  between neighboring regions are less than pi rad in magnitude.
"""

def create_parser(subparsers=None):
    synopsis = 'Unwrapping Error Correction with Bridging'
    epilog = REFERENCE + '\n' + TEMPLATE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis+NOTE, epilog=epilog, subparsers=subparsers)

    parser.add_argument('ifgram_file', type=str, help='interferograms file to be corrected')
    parser.add_argument('-r','--radius', dest='bridgePtsRadius', type=int, default=50,
                        help='radius of the end point of bridge to search area to get median representative value\n'+
                             'default: 50.')
    parser.add_argument('--ramp', dest='ramp', choices=['linear', 'quadratic'],
                          help='type of phase ramp to be removed before correction.')
    parser.add_argument('--water-mask','--wm', dest='waterMaskFile', type=str, help='path of water mask file.')
    parser.add_argument('-m', '--min-area', dest='connCompMinArea', type=float, default=2.5e3,
                        help='minimum region/area size of a single connComponent.')

    parser.add_argument('-t', '--template', dest='template_file', type=str,
                          help='template file with bonding point info, e.g.\n' +
                               'mintpy.unwrapError.yx = 283,1177,305,1247;350,2100,390,2200')

    parser.add_argument('-i','--in-dataset', dest='datasetNameIn', default='unwrapPhase',
                        help='name of dataset to be corrected, default: unwrapPhase')
    parser.add_argument('-o','--out-dataset', dest='datasetNameOut',
                        help='name of dataset to be written after correction, default: {}_bridging')
    parser.add_argument('--update', dest='update_mode', action='store_true',
                        help='Enable update mode: if unwrapPhase_unwCor dataset exists, skip the correction.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    # check input file type
    k = readfile.read_attribute(inps.ifgram_file)['FILE_TYPE']
    if k not in ['ifgramStack', '.unw']:
        raise ValueError('input file is not ifgramStack: {}'.format(k))

    # default output dataset name
    if not inps.datasetNameOut:
        inps.datasetNameOut = '{}_bridging'.format(inps.datasetNameIn)

    # discard water mask file is not found
    if inps.waterMaskFile and not os.path.isfile(inps.waterMaskFile):
        inps.waterMaskFile = None

    return inps


def read_template2inps(template_file, inps=None):
    """Read input template options into Namespace inps"""
    if not inps:
        inps = cmd_line_parse()
    inpsDict = vars(inps)
    print('read options from template file: '+os.path.basename(template_file))
    template = readfile.read_template(inps.template_file)
    template = ut.check_template_auto_value(template)

    key_list = [i for i in list(inpsDict.keys()) if key_prefix+i in template.keys()]
    for key in key_list:
        value = template[key_prefix+key]
        if key in ['update']:
            inpsDict[key] = value
        elif value:
            if key in ['waterMaskFile', 'ramp']:
                inpsDict[key] = value
            elif key in ['bridgePtsRadius']:
                inpsDict[key] = int(value)
            elif key in ['connCompMinArea']:
                inpsDict[key] = float(value)
    return inps


def run_or_skip(inps):
    print('-'*50)
    print('update mode: ON')
    flag = 'skip'

    # check output dataset
    with h5py.File(inps.ifgram_file, 'r') as f:
        if inps.datasetNameOut not in f.keys():
            flag = 'run'
            print('1) output dataset: {} NOT found.'.format(inps.datasetNameOut))
        else:
            print('1) output dataset: {} exists'.format(inps.datasetNameOut))
            ti = float(f[inps.datasetNameIn].attrs.get('MODIFICATION_TIME', os.path.getmtime(inps.ifgram_file)))
            to = float(f[inps.datasetNameOut].attrs.get('MODIFICATION_TIME', os.path.getmtime(inps.ifgram_file)))
            if ti > to:
                flag = 'run'
                print('2) output dataset is NOT newer than input dataset: {}.'.format(inps.datasetNameIn))
            else:
                print('2) output dataset is newer than input dataset: {}'.format(inps.datasetNameIn))

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
        changed_keys = [key for key in configKeys 
                        if str(inps_dict[key]) != atr.get(key_prefix+key, 'no')]
        if len(changed_keys) > 0:
            flag = 'run'
            print('3) NOT all key configuration parameters are the same: {}.'.format(configKeys))
            for key in changed_keys:
                print('\t{}\t: {} --> {}'.format(key, atr.get(key_prefix+key, 'no'), str(inps_dict[key])))
        else:
            print('3) all key configuration parameters are the same: {}.'.format(configKeys))

    # result
    print('run or skip: {}.'.format(flag))
    return flag


##########################################################################################
def run_unwrap_error_bridge(ifgram_file, water_mask_file, ramp_type=None, radius=50, cc_min_area=2.5e3,
                            ccName='connectComponent', dsNameIn='unwrapPhase',
                            dsNameOut='unwrapPhase_bridging'):
    """Run unwrapping error correction with bridging
    Parameters: ifgram_file     : str, path of ifgram stack file
                water_mask_file : str, path of water mask file
                ramp_type       : str, name of phase ramp to be removed during the phase jump estimation
                cc_min_area     : float, minimum region/area size
                ccName          : str, dataset name of connected components
                dsNameIn        : str, dataset name of unwrap phase to be corrected
                dsNameOut       : str, dataset name of unwrap phase to be saved after correction
    Returns:    ifgram_file     : str, path of ifgram stack file
    """
    print('-'*50)
    print('correct unwrapping error in {} with bridging ...'.format(ifgram_file))
    if ramp_type is not None:
        print('estimate and remove a {} ramp while calculating phase offset'.format(ramp_type))

    # read water mask
    if water_mask_file and os.path.isfile(water_mask_file):
        print('read water mask from file:', water_mask_file)
        water_mask = readfile.read(water_mask_file)[0]
    else:
        water_mask = None

    # file info
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
        print('open {} with r+ mode'.format(ifgram_file))
        with h5py.File(ifgram_file, 'r+') as f:
            print('input  dataset:', dsNameIn)
            print('output dataset:', dsNameOut)
            if dsNameOut in f.keys():
                ds = f[dsNameOut]
                print('access /{d} of np.float32 in size of {s}'.format(d=dsNameOut, s=shape_out))
            else:
                ds = f.create_dataset(dsNameOut,
                                      shape_out,
                                      maxshape=(None, None, None),
                                      chunks=True,
                                      compression=None)
                print('create /{d} of np.float32 in size of {s}'.format(d=dsNameOut, s=shape_out))

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
        print('close {} file.'.format(ifgram_file))

    if k == '.unw':
        # read unwrap phase
        unw = readfile.read(ifgram_file)[0]

        # read connected components
        cc_files0 = [ifgram_file+'.conncomp', os.path.splitext(ifgram_file)[0]+'_snap_connect.byt']
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
        out_file = '{}_unwCor{}'.format(os.path.splitext(ifgram_file)[0],
                                        os.path.splitext(ifgram_file)[1])
        print('writing >>> {}'.format(out_file))
        writefile.write(unw_cor, out_file=out_file, ref_file=ifgram_file)

    return ifgram_file


####################################################################################################
def main(iargs=None):
    # check inputs
    inps = cmd_line_parse(iargs)

    # update mode
    if inps.update_mode and run_or_skip(inps) == 'skip':
        return inps.ifgram_file

    start_time = time.time()
    # run bridging
    run_unwrap_error_bridge(
        ifgram_file=inps.ifgram_file,
        water_mask_file=inps.waterMaskFile,
        ramp_type=inps.ramp,
        radius=inps.bridgePtsRadius,
        cc_min_area=inps.connCompMinArea,
        dsNameIn=inps.datasetNameIn,
        dsNameOut=inps.datasetNameOut)

    # config parameter
    if os.path.splitext(inps.ifgram_file)[1] in ['.h5', '.he5']:
        print('add/update the following configuration metadata to file:')
        config_metadata = dict()
        for key in configKeys:
            config_metadata[key_prefix+key] = str(vars(inps)[key])
        ut.add_attribute(inps.ifgram_file, config_metadata, print_msg=True)

    m, s = divmod(time.time()-start_time, 60)
    print('\ntime used: {:02.0f} mins {:02.1f} secs\nDone.'.format(m, s))
    return


####################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
