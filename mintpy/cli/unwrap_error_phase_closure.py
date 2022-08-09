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


##########################################################################################
TEMPLATE1 = get_template_content('quick_overview')
TEMPLATE2 = get_template_content('correct_unwrap_error')

REFERENCE = """reference:
  Yunjun, Z., H. Fattahi, and F. Amelung (2019), Small baseline InSAR time series analysis:
  Unwrapping error correction and noise reduction, Computers & Geosciences, 133, 104331,
  doi:10.1016/j.cageo.2019.104331.
"""

EXAMPLE = """example:
  # correct phase unwrapping error with phase closure
  unwrap_error_phase_closure.py  ./inputs/ifgramStack.h5  --cc-mask maskConnComp.h5  -t smallbaselineApp.cfg   --update
  unwrap_error_phase_closure.py  ./inputs/ifgramStack.h5  --cc-mask maskConnComp.h5  --water-mask waterMask.h5 --update

  # calculate the number of non-zero closure phase
  unwrap_error_phase_closure.py  ./inputs/ifgramStack.h5  --action calculate
  unwrap_error_phase_closure.py  ./inputs/ifgramStack.h5  --action calculate  --water-mask waterMask.h5
"""

NOTE = """
  by exploiting the conservertiveness of the integer ambiguity of interferograms triplets.
  This method assumes:
  a. abundance of network: for interferogram with unwrapping error, there is
     at least of one triangular connection to form a closed circle; with more
     closed circles comes better constrain.
  b. majority rightness: most of interferograms have to be right (no unwrapping
     error) to correct the wrong minority. And if most of interferograms have 
     unwrapping errors, then the minor right interferograms will turn into wrong.
"""

def create_parser(subparsers=None):
    synopsis = 'Unwrapping Error Correction based on Phase Closure'
    epilog = REFERENCE + '\n' + TEMPLATE1 + '\n' + TEMPLATE2 + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis+NOTE, epilog=epilog, subparsers=subparsers)

    parser.add_argument('ifgram_file', help='interferograms file to be corrected')
    parser.add_argument('-c','--cc-mask', dest='cc_mask_file', default='maskConnComp.h5',
                        help='common connected components file, required for --action correct')
    parser.add_argument('-n','--num-sample', dest='numSample', type=int, default=100,
                        help='Number of randomly samples/pixels for each common connected component.')
    parser.add_argument('-m', '--min-area', dest='connCompMinArea', type=float, default=2.5e3,
                        help='minimum region/area size of a single connComponent.')

    parser.add_argument('-a','--action', dest='action', type=str, default='correct',
                        choices={'correct', 'calculate'},
                        help='action to take (default: %(default)s):\n'+
                             'correct   - correct phase unwrapping error\n'+
                             'calculate - calculate the number of non-zero closure phase')

    # IO
    parser.add_argument('-i','--in-dataset', dest='datasetNameIn', default='unwrapPhase',
                        help="name of dataset to be corrected, default: unwrapPhase")
    parser.add_argument('-o','--out-dataset', dest='datasetNameOut',
                        help='name of dataset to be written after correction, default: {}_phaseClosure')

    # mask
    mask = parser.add_argument_group('mask')
    mask.add_argument('--water-mask','--wm', dest='waterMaskFile', type=str,
                      help='path of water mask file.')
    mask.add_argument('-t', '--template', dest='template_file',
                      help='template file with options for setting.')

    parser.add_argument('--update', dest='update_mode', action='store_true',
                        help='Enable update mode: if unwrapPhase_phaseClosure dataset exists, skip the correction.')
    return parser


def cmd_line_parse(iargs=None):
    from ..utils import readfile
    from ..unwrap_error_phase_closure import read_template2inps
    from matplotlib import pyplot as plt

    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # backend setting for matplotlib
    plt.switch_backend('Agg')

    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    # check 1 input file type
    k = readfile.read_attribute(inps.ifgram_file)['FILE_TYPE']
    if k not in ['ifgramStack']:
        raise ValueError('input file is not ifgramStack: {}'.format(k))

    # check 2 cc_mask_file
    if inps.action == 'correct' and not os.path.isfile(inps.cc_mask_file):
        raise FileNotFoundError(inps.cc_mask_file)

    if not inps.datasetNameOut:
        inps.datasetNameOut = '{}_phaseClosure'.format(inps.datasetNameIn)

    # discard water mask file is not found
    if inps.waterMaskFile and not os.path.isfile(inps.waterMaskFile):
        inps.waterMaskFile = None

    return inps


####################################################################################################
def main(iargs=None):
    from mintpy.utils import plot as pp
    from ..unwrap_error_phase_closure import (
        run_or_skip,
        get_common_region_int_ambiguity,
        run_unwrap_error_phase_closure,
        calc_num_triplet_with_nonzero_integer_ambiguity,
    )

    inps = cmd_line_parse(iargs)
    start_time = time.time()

    if inps.action == 'correct':
        # update mode
        if inps.update_mode and run_or_skip(inps) == 'skip':
            return inps.ifgram_file

        # solve integer ambiguity for common connected components
        common_regions = get_common_region_int_ambiguity(ifgram_file=inps.ifgram_file,
                                                         cc_mask_file=inps.cc_mask_file,
                                                         water_mask_file=inps.waterMaskFile,
                                                         num_sample=inps.numSample,
                                                         dsNameIn=inps.datasetNameIn,
                                                         cc_min_area=inps.connCompMinArea)

        # apply the integer ambiguity from common conn comp to the whole ifgram
        if len(common_regions) == 0:
            print('skip phase closure correction ...')
        else:
            run_unwrap_error_phase_closure(inps.ifgram_file, common_regions,
                                           water_mask_file=inps.waterMaskFile,
                                           dsNameIn=inps.datasetNameIn,
                                           dsNameOut=inps.datasetNameOut)

    else:
        # calculate the number of triplets with non-zero integer ambiguity
        out_file = calc_num_triplet_with_nonzero_integer_ambiguity(inps.ifgram_file,
                                                                   mask_file=inps.waterMaskFile,
                                                                   dsName=inps.datasetNameIn,
                                                                   update_mode=inps.update_mode)
        # for debug
        debug_mode = False
        if debug_mode:
            pp.plot_num_triplet_with_nonzero_integer_ambiguity(out_file)

    m, s = divmod(time.time()-start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs\nDone.'.format(m, s))


####################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
