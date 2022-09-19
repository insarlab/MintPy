############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import os
import sys
import time
from mintpy.defaults.template import get_template_content
from mintpy.utils.arg_utils import create_argument_parser


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
    from mintpy.utils import readfile
    from mintpy.unwrap_error_bridging import read_template2inps

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


####################################################################################################
def main(iargs=None):
    from mintpy.utils import utils as ut
    from mintpy.unwrap_error_bridging import run_unwrap_error_bridge, run_or_skip, config_keys, key_prefix

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
        for key in config_keys:
            config_metadata[key_prefix+key] = str(vars(inps)[key])
        ut.add_attribute(inps.ifgram_file, config_metadata, print_msg=True)

    m, s = divmod(time.time()-start_time, 60)
    print('\ntime used: {:02.0f} mins {:02.1f} secs\nDone.'.format(m, s))


####################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
