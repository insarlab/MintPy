#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import os
import sys

from mintpy.utils.arg_utils import create_argument_parser

###########################  Sub Function  ##############################
BL_LIST = """
070106     0.0   0.03  0.0000000  0.00000000000 2155.2 /scratch/SLC/070106/
070709  2631.9   0.07  0.0000000  0.00000000000 2155.2 /scratch/SLC/070709/
070824  2787.3   0.07  0.0000000  0.00000000000 2155.2 /scratch/SLC/070824/
"""

DATE12_LIST = """
20070709_20100901
20070709_20101017
20070824_20071009
"""

TEMPLATE = """
mintpy.network.maskFile  = auto  #[file name, no], auto for waterMask.h5 or no for all pixels
mintpy.network.aoiYX     = auto  #[y0:y1,x0:x1 / no], auto for no, area of interest for coherence calculation
mintpy.network.aoiLALO   = auto  #[lat0:lat1,lon0:lon1 / no], auto for no - use the whole area
"""

EXAMPLE = """example:
  plot_network.py inputs/ifgramStack.h5
  plot_network.py inputs/ifgramStack.h5 -t smallbaselineApp.cfg --nodisplay   #Save figures to files without display
  plot_network.py inputs/ifgramStack.h5 -t smallbaselineApp.cfg --show-kept   #Do not plot dropped ifgrams
  plot_network.py inputs/ifgramStack.h5 -d tbase -v 0 365.25 -c RdYlBu_r      #Color-code lines by temporal      baseline
  plot_network.py inputs/ifgramStack.h5 -d pbase -v 0 180    -c RdYlBu_r      #Color-code lines by perpendicular baseline
  plot_network.py coherenceSpatialAvg.txt

  # offsetSNR
  plot_network.py inputs/ifgramStack.h5 -d offsetSNR -v 0 20 --cmap-vlist 0 0.2 1
"""


def create_parser(subparsers=None):
    synopsis = 'Display Network of Interferograms'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', help='file with network information, ifgramStack.h5 or coherenceSpatialAvg.txt')
    parser.add_argument('--show-kept', dest='disp_drop', action='store_false',
                        help='display kept interferograms only, without dropped interferograms')

    # Color-code coherence/baseline in the network/matrix plot
    color = parser.add_argument_group('Color-code network/matrix plot',
                                      'color-code phase/offset pairs with coherence/baseline in network/matrix plot')
    color.add_argument('-d', '--dset', type=str, dest='dsetName', default='coherence',
                       choices={'coherence','offsetSNR','pbase','tbase'},
                       help='dataset used to calculate the mean. (default: %(default)s)')
    color.add_argument('-v', '--vlim', nargs=2, type=float, default=(0.2, 1.0), help='display range')

    color.add_argument('-t', '--template', dest='template_file',
                       help='template file with options below:\n'+TEMPLATE)
    color.add_argument('--mask', dest='maskFile', default='waterMask.h5',
                       help='mask file used to calculate the coherence. Default: waterMask.h5 or None.')

    color.add_argument('-c', '--colormap', dest='cmap_name', default='RdBu_truncate',
                       help='colormap name for the network display. Default: RdBu_truncate')
    color.add_argument('--cmap-vlist', dest='cmap_vlist', type=float, nargs=3, default=[0.2, 0.4, 1.0],
                       help='normalized start/jump/end value for truncated colormap (default: %(default)s).')

    # Figure  Setting
    fig = parser.add_argument_group('Figure', 'Figure settings for display')
    fig.add_argument('--fs', '--fontsize', type=int,
                     default=12, help='font size in points')
    fig.add_argument('--lw', '--linewidth', dest='linewidth',
                     type=int, default=2, help='line width in points')
    fig.add_argument('--mc', '--markercolor', dest='markercolor',
                     default='orange', help='marker color')
    fig.add_argument('--ms', '--markersize', dest='markersize',
                     type=int, default=8, help='marker size in points (default: %(default)s).')
    fig.add_argument('--every-year', dest='every_year', type=int,
                     default=1, help='number of years per major tick on x-axis')

    fig.add_argument('--dpi', dest='fig_dpi', type=int, default=1200,
                     help='DPI - dot per inch - for display/write')
    fig.add_argument('--figsize', dest='fig_size', type=float, nargs=2,
                     help='figure size in inches - width and length.')
    fig.add_argument('--notitle', dest='disp_title', action='store_false',
                     help='Do not display figure title.')
    fig.add_argument('--number', dest='number', type=str,
                     help='number mark to be plot at the corner of figure.')
    fig.add_argument('--nosplit-cmap', dest='split_cmap', action='store_false',
                     help='do not split colormap for coherence color')

    fig.add_argument('--save', dest='save_fig',
                     action='store_true', help='save the figure')
    fig.add_argument('--nodisplay', dest='disp_fig',
                     action='store_false', help='save and do not display the figure')
    return parser


def cmd_line_parse(iargs=None):
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    from mintpy.utils import readfile

    # save argv (to check the manually specified arguments)
    # use iargs        for python call
    # use sys.argv[1:] for command line call
    inps.argv = iargs if iargs else sys.argv[1:]

    # check: input file type
    if inps.file.endswith(('.h5','.he5')):
        k = readfile.read_attribute(inps.file)['FILE_TYPE']
        if k != 'ifgramStack':
            raise ValueError('input HDF5 file is NOT ifgramStack.')

    # check: --nodisplay option
    if not inps.disp_fig:
        inps.save_fig = True

    # check: -t / --template option
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)
    else:
        inps.template = {}

    # check: --mask option (file existence)
    if not os.path.isfile(inps.maskFile):
        inps.maskFile = None

    return inps


def read_template2inps(template_file, inps):
    """Read input template options into Namespace inps"""
    print('read options from template file: '+os.path.basename(template_file))

    from mintpy.utils import readfile, utils1 as ut

    inps.template = readfile.read_template(inps.template_file)
    inps.template = ut.check_template_auto_value(inps.template)

    # coherence-based network modification
    prefix = 'mintpy.network.'
    key = prefix+'maskFile'
    if key in inps.template.keys():
        if inps.template[key]:
            inps.maskFile = inps.template[key]

    return inps


##########################  Main Function  ##############################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.plot_network import plot_network

    # run
    plot_network(inps)


#########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
