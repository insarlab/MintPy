############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import sys
from mintpy.utils.arg_utils import create_argument_parser


#################################  Usage  ####################################
EXAMPLE = """example:
  spatial_average.py inputs/ifgramStack.h5  -d coherence -m maskConnComp.h5
  spatial_average.py timeseries_ERA5_demErr.h5 -m maskTempCoh.h5
"""

def create_parser(subparsers=None):
    synopsis = 'Calculate average in space'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', help='File to calculate spatial average')
    parser.add_argument('-d', '--dset', '--dataset', dest='datasetName',
                        help='dataset used to calculate, for ifgramStack file only.')
    parser.add_argument('-m', '--mask', dest='mask_file',
                        help='Mask file for the calculation')
    parser.add_argument('--nodisplay', dest='disp_fig',
                        action='store_false', help='save and do not display the figure')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


#############################  Main Function  ################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.spatial_average import run_spatial_average

    # run
    run_spatial_average(inps)


##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
