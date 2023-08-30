#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import os
import sys

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
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    from mintpy.utils import readfile

    # check: -t / --template option (read template content)
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    # check: input file type (ifgramStack is required)
    ftype = readfile.read_attribute(inps.ifgram_file)['FILE_TYPE']
    if ftype not in ['ifgramStack']:
        raise ValueError(f'input file is not ifgramStack: {ftype}')

    # check: --cc-mask option (required for "--action correct")
    if inps.action == 'correct' and not os.path.isfile(inps.cc_mask_file):
        raise FileNotFoundError(inps.cc_mask_file)

    # check: --water-mask option (file existence)
    if inps.waterMaskFile and not os.path.isfile(inps.waterMaskFile):
        inps.waterMaskFile = None

    # default: --out-dataset option
    if not inps.datasetNameOut:
        inps.datasetNameOut = f'{inps.datasetNameIn}_phaseClosure'

    return inps


def read_template2inps(template_file, inps):
    """Read input template options into Namespace inps"""
    print('read options from template file: '+os.path.basename(inps.template_file))

    from mintpy.unwrap_error_phase_closure import key_prefix
    from mintpy.utils import readfile, utils1 as ut

    inpsDict = vars(inps)
    template = readfile.read_template(template_file)
    template = ut.check_template_auto_value(template)

    key_list = [i for i in list(inpsDict.keys()) if key_prefix+i in template.keys()]
    for key in key_list:
        value = template[key_prefix+key]
        if value:
            if key in ['waterMaskFile']:
                inpsDict[key] = value
            elif key in ['numSample']:
                inpsDict[key] = int(value)
            elif key in ['connCompMinArea']:
                inpsDict[key] = float(value)

    return inps


####################################################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.unwrap_error_phase_closure import (
        run_unwrap_error_phase_closure,
    )

    # run
    run_unwrap_error_phase_closure(inps)


####################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
