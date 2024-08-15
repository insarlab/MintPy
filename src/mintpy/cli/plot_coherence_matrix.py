#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import os
import sys

from mintpy.utils.arg_utils import create_argument_parser

###########################  Sub Function  #############################
EXAMPLE = """example:
  plot_coherence_matrix.py inputs/ifgramStack.h5
  plot_coherence_matrix.py inputs/ifgramStack.h5 --yx 277 1069
  plot_coherence_matrix.py inputs/ifgramStack.h5 --lalo -0.8493 -91.1510 -c RdBu

  # left: map view
  plot_coherence_matrix.py inputs/ifgramStack.h5 --view-cmd "view.py {} --dem inputs/gsi10m.dem.wgs84"
  plot_coherence_matrix.py inputs/ifgramStack.h5 --view-cmd 'view.py {} --wrap --wrap-range -3 3"
  plot_coherence_matrix.py inputs/ifgramStack.h5 --view-cmd 'view.py {} --sub-x 900 1400 --sub-y 0 500'

  # right: matrix view
  # show color jump same as the coherence threshold in network inversion with pixel-wised masking
  plot_coherence_matrix.py inputs/ifgramStack.h5 --cmap-vlist 0 0.4 1
"""


def create_parser(subparsers=None):
    synopsis = 'Plot the coherence matrix of one pixel (interactive)'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('ifgram_file', help='interferogram stack file')
    parser.add_argument('--yx', type=int, metavar=('Y', 'X'), nargs=2,
                        help='Point of interest in y(row)/x(col)')
    parser.add_argument('--lalo', type=float, metavar=('LAT','LON'), nargs=2,
                        help='Point of interest in lat/lon')
    parser.add_argument('--lookup','--lut', dest='lookup_file',
                        help='Lookup file to convert lat/lon into y/x')
    parser.add_argument('-c','--cmap', dest='cmap_name', default='RdBu_truncate',
                        help='Colormap for coherence matrix.\nDefault: RdBu_truncate')
    parser.add_argument('--cmap-vlist', dest='cmap_vlist', type=float, nargs=3, default=[0.0, 0.7, 1.0],
                        help='start/jump/end fraction for truncated colormap. Default: 0.0 0.7 1.0')
    parser.add_argument('--figsize','--fs', dest='fig_size', metavar=('WID', 'LEN'), type=float, nargs=2,
                        help='figure size in inches. Default: [8, 4]')

    parser.add_argument('--img-file', dest='img_file',
                        help='dataset to show in map to facilitate point selection. Default: velocity.h5')
    parser.add_argument('--view-cmd', dest='view_cmd', default='view.py {} --wrap --noverbose ',
                        help='view.py command to plot the input map file\n'+
                             'Default: view.py img_file --wrap --noverbose')

    # aux files
    parser.add_argument('--tcoh', dest='tcoh_file', default='temporalCoherence.h5',
                        help='temporal coherence file.')
    parser.add_argument('-t','--template', dest='template_file',
                        help='temporal file.')

    parser.add_argument('--save', dest='save_fig',
                        action='store_true', help='save the figure')
    parser.add_argument('--nodisplay', dest='disp_fig',
                        action='store_false', help='save and do not display the figure')
    parser.add_argument('--noverbose', dest='print_msg', action='store_false',
                        help='Disable the verbose message printing.')
    return parser


def cmd_line_parse(iargs=None):
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # save argv (to check the manually specified arguments)
    # use iargs        for python call
    # use sys.argv[1:] for command line call
    inps.argv = iargs if iargs else sys.argv[1:]

    # default: auxiliary file paths (velocity and template)
    mintpy_dir = os.path.dirname(os.path.dirname(inps.ifgram_file))
    if not inps.img_file:
        inps.img_file = os.path.join(mintpy_dir, 'velocity.h5')
    if not inps.template_file:
        inps.template_file = os.path.join(mintpy_dir, 'smallbaselineApp.cfg')

    # check: existence of auxliary files
    if not os.path.isfile(inps.img_file):
        raise SystemExit(f'ERROR: input image file not found: {inps.img_file}')
    if not os.path.isfile(inps.tcoh_file):
        inps.tcoh_file = None
    if not os.path.isfile(inps.template_file):
        inps.template_file = None

    # check: --nodisplay option
    if not inps.disp_fig:
        inps.save_fig = True

    return inps


##########################  Main Function  ##############################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.plot_coherence_matrix import coherenceMatrixViewer

    # run
    obj = coherenceMatrixViewer(inps)
    obj.open()
    obj.plot()
    obj.fig.canvas.mpl_disconnect(obj.cid)


############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
