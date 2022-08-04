#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import sys
import argparse
import numpy as np
from matplotlib import pyplot as plt, ticker

from mintpy.utils import arg_utils, readfile, utils as ut, plot as pp
from mintpy import view


#####################################################################
# Only one line is supported right now.
GMT_FILE = """GMT xy file, i.e. transect_lonlat.xy:
>
131.1663    33.1157
131.2621    33.0860
"""

EXAMPLE = """example:
  plot_transection.py velocity.h5 --start-yx 5290 5579 --end-yx 12177 482
  plot_transection.py velocity.h5 --start-lalo 30.125 129.988 --end-lalo 30.250 130.116
  plot_transection.py velocity.h5 --line-file  transect_lonlat.xy --dem gsi10m.dem

  # multiple files
  plot_transection.py AlosA*/velocity.h5 AlosD*/velocity.h5 --off 2
  plot_transection.py Kirishima2017*.h5 Kirishima2008*.h5 --off 0 0 10 10
  plot_transection.py Kirishima2017*.h5 Kirishima2008*.h5 --off 0 0 10 10 --lalo0 31.947 130.843 --lalo1 31.947 130.860

  # interactive plot: click two points to draw a profile
"""


def create_parser(subparsers=None):
    synopsis = 'Generate transect/profile along a line'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = arg_utils.create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', nargs='+',
                        help='input file to show transection')
    parser.add_argument('--dset', dest='dset', help='Dataset name to read')
    parser.add_argument('-v','--vlim', dest='vlim', nargs=2, metavar=('VMIN', 'VMAX'), type=float,
                        help='Display limits for matrix plotting.')
    parser.add_argument('--offset','--off', dest='offset', type=float, nargs='+', default=[0.05],
                        help='offset between transects [for multiple files only; default: %(default)s m].\n'
                             'number of input offsets should be:\n'
                             '    1 - same (sequential) offset between adjacent transects OR\n'
                             '    num_file - different (cumulative) offset for each file, starting from 0.')
    parser.add_argument('--noverbose', dest='print_msg', action='store_false',
                        help='Disable the verbose message printing.')

    lines = parser.add_argument_group('Profile location', 'Start/end points of profile')
    lines.add_argument('--start-yx','--yx0', dest='start_yx', metavar=('Y0', 'X0'), type=int, nargs=2,
                       help='start point of the profile in pixel number [y, x]')
    lines.add_argument('--end-yx','--yx1', dest='end_yx', metavar=('Y1', 'X1'), type=int, nargs=2,
                       help='end   point of the profile in pixel number [y, x]')
    lines.add_argument('--start-lalo','--lalo0', dest='start_lalo', metavar=('LAT0', 'LON0'), type=float, nargs=2,
                       help='start point of the profile in [lat, lon]')
    lines.add_argument('--end-lalo','--lalo1', dest='end_lalo', metavar=('LAT1', 'LON1'), type=float, nargs=2,
                       help='end   point of the profile in [lat, lon]')
    lines.add_argument('--line-file', dest='lola_file',
                       help='file with start and end point info in lon lat, same as GMT format.\n'+GMT_FILE)

    lines.add_argument('--interpolation', default='nearest', choices=['nearest', 'bilinear', 'cubic'],
                       help='interpolation method while extacting profile along the line. Default: nearest.')
    lines.add_argument('--ms', '--markersize', dest='marker_size', type=float, default=2.0,
                       help='Point marker size. Default: 2.0')

    parser = arg_utils.add_figure_argument(parser)
    parser = arg_utils.add_save_argument(parser)
    return parser


def cmd_line_parse(iargs=None):
    inps = create_parser().parse_args(args=iargs)

    # save argv (to check the manually specified arguments)
    # use iargs        for python call
    # use sys.argv[1:] for command line call
    inps.argv = iargs if iargs else sys.argv[1:]

    if inps.outfile or not inps.disp_fig:
        inps.save_fig = True

    # input file info
    inps.file = ut.get_file_list(inps.file)
    inps.atr = readfile.read_attribute(inps.file[0])
    inps.coord = ut.coordinate(inps.atr)
    inps.num_file = len(inps.file)

    if inps.num_file > 1:
        # a) one input: it's interval between adjacent files
        if len(inps.offset) == 1:
            inps.offset = np.ones(inps.num_file, dtype=np.float32) * inps.offset
            inps.offset[0] = 0.
            inps.offset = np.cumsum(inps.offset)

        # b) multiple input: it's exact offset of all files
        elif len(inps.offset) == inps.num_file:
            inps.offset = np.array(inps.offset, dtype=np.float32)

        # c) do not support any other numbers of inputs
        else:
            msg = 'input number of offsets: {}.'.format(len(inps.offset))
            msg += '\nIt should be 1 or number of files: {}'.format(inps.num_file)
            raise ValueError(msg)
    else:
        # disable offset for single input file
        inps.offset = np.array([0], dtype=np.float32)

    if not inps.dset:
        inps.dset = readfile.get_slice_list(inps.file[0])[0]

    # lola_file --> start/end_lalo
    if inps.lola_file:
        inps.start_lalo, inps.end_lalo = read_lonlat_file(inps.lola_file)

    # start/end_lalo --> start/end_yx
    if inps.start_lalo and inps.end_lalo:
        [y0, y1] = inps.coord.lalo2yx([inps.start_lalo[0], inps.end_lalo[0]], coord_type='lat')
        [x0, x1] = inps.coord.lalo2yx([inps.start_lalo[1], inps.end_lalo[1]], coord_type='lon')
        inps.start_yx = [y0, x0]
        inps.end_yx = [y1, x1]

    # verbose print using --noverbose option
    global vprint
    vprint = print if inps.print_msg else lambda *args, **kwargs: None

    if not inps.disp_fig:
        plt.switch_backend('Agg')
    return inps


#####################################################################
def read_lonlat_file(lonlat_file):
    """Read Start/End lat/lon from lonlat text file in gmt format.
    Inputs:
        lonlat_file : text file in gmt lonlat point file
    Outputs:
        start/end_lalo : list of 2 float
    """
    fll = open(lonlat_file, 'r')
    lines = fll.read().splitlines()
    [lon0, lat0] = [float(i) for i in lines[1].split()]
    [lon1, lat1] = [float(i) for i in lines[2].split()]
    fll.close()

    start_lalo = [lat0, lon0]
    end_lalo = [lat1, lon1]
    return start_lalo, end_lalo


def get_view_cmd(iargs):
    """Assemble view.py command line from input arguments"""
    # define ALL parsing options from create_parser() that are common to view.py
    parser = argparse.ArgumentParser(description='view.py parser')
    parser.add_argument('-v','--vlim', dest='vlim', nargs=2, metavar=('VMIN', 'VMAX'), type=float,
                        help='Display limits for matrix plotting.')
    parser.add_argument('--noverbose', dest='print_msg', action='store_false',
                        help='Disable the verbose message printing.')
    parser = arg_utils.add_figure_argument(parser)
    parser = arg_utils.add_save_argument(parser)

    # get args that are applicable to view.py
    unique_args = parser.parse_known_args(iargs)[1]
    view_args = [i for i in iargs if i not in unique_args]

    # assemble view.py command line
    inps = cmd_line_parse(iargs)
    view_cmd = 'view.py {} '.format(inps.file[0])
    if inps.dset:
        view_cmd += ' {} '.format(inps.dset)
    view_cmd += ' '.join(view_args)
    return view_cmd


#####################################################################
class transectionViewer():
    """class for plot_transection
    Example:
        from mintpy.plot_transection import transectionViewer
        cmd = 'plot_transection.py velocity.h5 --noverbose --start-yx 10 10 --end-yx 200 300'
        obj = transectionViewer(cmd)
        obj.configure()
        obj.plot()
    """

    def __init__(self, cmd=None, iargs=None):
        if cmd:
            iargs = cmd.split()[1:]
        self.cmd = cmd
        self.iargs = iargs

        # figure variables
        self.figname = 'Transection'
        self.fig = None
        self.ax_img = None
        self.ax_txn = None

        self.img = None
        self.line_ann = None
        self.pts_idx = 0
        return

    def configure(self):
        # copy inps to self object
        inps = cmd_line_parse(self.iargs)
        for key, value in inps.__dict__.items():
            setattr(self, key, value)

        # copy inps from view.py to self object
        view_cmd = get_view_cmd(self.argv)
        self.data_img, atr, inps_view = view.prep_slice(view_cmd)
        for key, value in inps_view.__dict__.items():
            setattr(self, key, value)

        self.offset *= self.disp_scale  # due to unit change from view.prep_slice()

        # do not update the following setting from view.py
        self.file = inps.file
        self.dset = inps.dset
        self.fig_size = inps.fig_size

        # auto figure size
        if not self.fig_size:
            length, width = int(self.atr['LENGTH']), int(self.atr['WIDTH'])
            fig_size = pp.auto_figure_size((length, width), disp_cbar=True)
            self.fig_size = [fig_size[0] + fig_size[1], fig_size[1]]
        return

    def plot(self):
        # Read data for transection
        self.data_list = []
        self.atr_list = []
        for fname in self.file:
            data, atr = readfile.read(fname, datasetName=self.dset)
            data = pp.scale_data2disp_unit(data, metadata=atr, disp_unit=self.disp_unit)[0]
            self.data_list.append(data)
            self.atr_list.append(atr)

        # Figure
        self.fig, (self.ax_img, self.ax_txn) = plt.subplots(1, 2, num=self.figname, figsize=self.fig_size)

        # Axes 1 - map with view.prep/plot_slice()
        self.ax_img = view.plot_slice(self.ax_img, self.data_img, self.atr, self)[0]

        # Axes 2 - transection
        self.ax_txn.yaxis.tick_right()
        self.ax_txn.yaxis.set_label_position("right")

        # plot initial input transect
        if self.start_yx and self.end_yx:
            self.draw_line(self.start_yx, self.end_yx)
            self.draw_transection(self.start_yx, self.end_yx, self.start_lalo, self.end_lalo)

        self.fig.subplots_adjust(left=0.08, wspace=0.25, bottom=0.15)

        # save
        if self.save_fig:
            outfile = '{}.pdf'.format(self.outfile_base)
            self.fig.savefig(outfile, bbox_inches='tight', transparent=True, dpi=self.fig_dpi)
            vprint('saved transect to', outfile)

        self.cid = self.fig.canvas.mpl_connect('button_release_event', self.select_point)

        if self.disp_fig:
            vprint('showing ...')
            plt.show()
        return


    ##---------- event function
    def select_point(self, event):
        """Event handling function for points selection"""
        if event.inaxes == self.ax_img:
            # get row/col number
            if 'Y_FIRST' in self.atr.keys():
                lalo = [event.ydata, event.xdata]
                yx = self.coord.geo2radar(event.ydata, event.xdata, print_msg=False)[0:2]
            else:
                lalo = None
                yx = [int(event.ydata+0.5), int(event.xdata+0.5)]

            # insert selected points into self.start/end_yx member
            # print('pts_idx: {}'.format(self.pts_idx)) #for debug
            if self.pts_idx == 0:
                self.start_lalo = lalo
                self.start_yx = yx
            else:
                self.end_lalo = lalo
                self.end_yx = yx

            # update transection for every two clicks
            self.pts_idx += 1
            if self.pts_idx >= 2:
                self.pts_idx = 0
                self.draw_line(self.start_yx, self.end_yx)
                self.draw_transection(self.start_yx, self.end_yx, self.start_lalo, self.end_lalo)
        return


    ##---------- plot functions
    def draw_line(self, start_yx, end_yx):
        """Draw the transect line in the map axes"""
        # erase existing line
        if self.line_ann is not None:
            self.line_ann.remove()

        # convert coordinates accordingly
        if 'Y_FIRST' in self.atr.keys():
            ys = self.coord.yx2lalo([start_yx[0], end_yx[0]], coord_type='y')
            xs = self.coord.yx2lalo([start_yx[1], end_yx[1]], coord_type='x')
        else:
            ys = [start_yx[0], end_yx[0]]
            xs = [start_yx[1], end_yx[1]]

        # plot
        line = self.ax_img.plot(xs, ys, 'k--', alpha=0)[0]
        self.line_ann = pp.add_arrow(line, position=xs[1])

        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()
        return


    def draw_transection(self, start_yx, end_yx, start_lalo=None, end_lalo=None):
        """Plot the transect as dots"""
        self.ax_txn.cla()

        # loop - extract transection data
        txn_list = []
        min_dist = 0
        for i in range(self.num_file):
            # get transection data
            if start_lalo is not None:
                # use lat/lon whenever it's possible to support files with different resolutions
                txn = ut.transect_lalo(self.data_list[i],
                                       self.atr_list[i],
                                       start_lalo, end_lalo,
                                       interpolation=self.interpolation)
            else:
                txn = ut.transect_yx(self.data_list[i],
                                     self.atr_list[i],
                                     start_yx, end_yx,
                                     interpolation=self.interpolation)

            # save txn
            txn_list.append(txn)
            min_dist = max(min_dist, txn['distance'][0])

        # loop - plot transection
        for i, txn in enumerate(txn_list):
            # distance unit and scaling
            if txn.get('distance_unit', 'm') == 'pixel':
                dist_scale = 1.0
                dist_unit = 'pixel'
            else:
                dist_scale = 0.001
                dist_unit = 'km'

            # plot
            # update distance values by excluding the commonly masked out pixels in the begining
            self.ax_txn.scatter(x=(txn['distance'] - min_dist) * dist_scale,
                                y=txn['value'] - self.offset[i],
                                c=pp.mplColors[i],
                                s=self.marker_size**2)

        y0, x0, y1, x1 = start_yx + end_yx
        self.outfile_base = f'transect_Y{y0}X{x0}_Y{y1}X{x1}'

        # title
        msg = f'y/x: ({y0}, {x0}) --> ({y1}, {x1})'
        if 'Y_FIRST' in self.atr.keys():
            lat0, lon0 = self.coord.radar2geo(start_yx[0], start_yx[1])[0:2]
            lat1, lon1 = self.coord.radar2geo(end_yx[0], end_yx[1])[0:2]
            msg += f'\nlat/lon: ({lat0:.4f}, {lon0:.4f}) --> ({lat1:.4f}, {lon1:.4f})'
        self.ax_txn.set_title(msg, fontsize=self.font_size)

        # axis format
        self.ax_txn.tick_params(which='both', direction='in', labelsize=self.font_size,
                                bottom=True, top=True, left=True, right=True)
        self.ax_txn.yaxis.set_minor_locator(ticker.AutoMinorLocator(10))
        self.ax_txn.set_ylabel(self.disp_unit, fontsize=self.font_size)
        self.ax_txn.set_xlabel(f'Distance [{dist_unit}]', fontsize=self.font_size)
        self.ax_txn.set_xlim(0, (txn['distance'][-1] - min_dist) * dist_scale)

        # update figure
        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()
        return



############################ Main ###################################
def main(iargs=None):
    obj = transectionViewer(iargs=iargs)
    obj.configure()
    obj.plot()
    obj.fig.canvas.mpl_disconnect(obj.cid)
    return


#####################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
