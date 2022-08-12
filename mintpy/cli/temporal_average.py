############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import sys
import time
from mintpy.utils.arg_utils import create_argument_parser


#################################  Usage  ####################################
EXAMPLE = """example:
  temporal_average.py ./inputs/ifgramStack.h5 -d unwrapPhase -o avgPhaseVelocity.h5
  temporal_average.py ./inputs/ifgramStack.h5 -d coherence   -o avgSpatialCoh.h5
"""

def create_parser(subparsers=None):
    synopsis = 'Calculate temporal average (stacking) of multi-temporal datasets'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', type=str, help='input file with multi-temporal datasets')
    parser.add_argument('-d', '--ds', '--dataset', dest='datasetName', default='coherence',
                        help='dataset name to be averaged, for file with multiple dataset family,\n'+
                        'e.g. ifgramStack.h5\n' +
                        'default: coherence')
    parser.add_argument('-o', '--outfile', help='output file name')
    parser.add_argument('--update', dest='update_mode', action='store_true',
                        help='Enable update checking for --nonzero option.')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


#############################  Main Function  ################################
def main(iargs=None):
    from ..utils import utils as ut
    from ..temporal_average import check_output_filename, run_or_skip

    start_time = time.time()
    inps = cmd_line_parse(iargs)

    inps.outfile = check_output_filename(inps)

    if inps.update_mode and run_or_skip(inps) == 'skip':
        return inps.outfile

    ut.temporal_average(inps.file, datasetName=inps.datasetName, outFile=inps.outfile)

    m, s = divmod(time.time()-start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs\n'.format(m, s))
    return


##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
