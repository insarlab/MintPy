#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2018, Zhang Yunjun                          #
# Author:  Zhang Yunjun, Nov 2018                          #
############################################################


import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pysar.objects import ifgramStack
from pysar.utils import readfile, plot as pp, utils as ut
from pysar import view


###########################  Sub Function  #############################
EXAMPLE = """example:
  plot_coherence_matrix.py  INPUTS/ifgramStack.h5  --yx 277 1069
  plot_coherence_matrix.py  INPUTS/ifgramStack.h5  --yx 277 1069  --map-file velocity.h5 --map-plot-cmd 'view.py {} --wrap --wrap-range -3 3 --sub-x 900 1400 --sub-y 0 500'
  plot_coherence_matrix.py  INPUTS/ifgramStack.h5  --lalo -0.8493 -91.1510 -c RdBu
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Plot the coherence matrix of one pixel (interactive)',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('ifgram_file', help='interferogram stack file')
    parser.add_argument('--yx', type=int, metavar=('Y', 'X'), nargs=2, 
                        help='Point of interest in y(row)/x(col)')
    parser.add_argument('--lalo', type=float, metavar=('LAT','LON'), nargs=2,
                        help='Point of interest in lat/lon')
    parser.add_argument('--lookup','--lut', dest='lookup_file',
                        help='Lookup file to convert lat/lon into y/x')
    parser.add_argument('-c','--cmap', dest='colormap', default='truncate_RdBu',
                        help='Colormap for coherence matrix. Default: truncate_RdBu')
    parser.add_argument('--figsize','--fs', dest='fig_size', metavar=('WID', 'LEN'), type=float, nargs=2,
                        help='figure size in inches. Default: [8, 4]')

    parser.add_argument('--img-file', dest='img_file', default='velocity.h5',
                        help='dataset to show in map to facilitate point selection')
    parser.add_argument('--view-cmd', dest='view_cmd', default='view.py {} --wrap --noverbose ',
                        help='view.py command to plot the input map file')

    parser.add_argument('--save', dest='save_fig',
                        action='store_true', help='save the figure')
    parser.add_argument('--nodisplay', dest='disp_fig',
                        action='store_false', help='save and do not display the figure')
    parser.add_argument('--noverbose', dest='print_msg', action='store_false',
                        help='Disable the verbose message printing.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    if not os.path.isfile(inps.img_file):
        raise SystemExit('ERROR: input image file not found: {}'.format(inps.img_file))

    # verbose print using --noverbose option
    global vprint
    vprint = print if inps.print_msg else lambda *args, **kwargs: None

    if not inps.disp_fig:
        inps.save_fig = True
        plt.switch_backend('Agg')
    return inps


def read_network_info(inps):
    k = readfile.read_attribute(inps.ifgram_file)['FILE_TYPE']
    if k != 'ifgramStack':
        raise ValueError('input file {} is not ifgramStack: {}'.format(inps.ifgram_file, k))

    obj = ifgramStack(inps.ifgram_file)
    obj.open(print_msg=inps.print_msg)
    inps.date12_list    = obj.get_date12_list(dropIfgram=False)
    date12_kept = obj.get_date12_list(dropIfgram=True)
    inps.ex_date12_list = sorted(list(set(inps.date12_list) - set(date12_kept)))
    inps.date_list = obj.get_date_list(dropIfgram=False)
    vprint('number of all     interferograms: {}'.format(len(inps.date12_list)))
    vprint('number of dropped interferograms: {}'.format(len(inps.ex_date12_list)))
    vprint('number of kept    interferograms: {}'.format(len(inps.date12_list) - len(inps.ex_date12_list)))
    vprint('number of acquisitions: {}'.format(len(inps.date_list)))

    if inps.lalo:
        if not inps.lookup_file:            
            lookup_file = os.path.join(os.path.dirname(inps.ifgram_file), 'geometry*.h5')
            inps.lookup_file = ut.get_lookup_file(filePattern=lookup_file)
        coord = ut.coordinate(obj.metadata, lookup_file=inps.lookup_file)
        inps.yx = coord.geo2radar(inps.lalo[0], inps.lalo[1])[0:2]

    if not inps.yx:
        inps.yx = (obj.refY, obj.refX)
        vprint('plot initial coherence matrix at reference pixel: {}'.format(inps.yx))
    return inps


class networkViewer():
    """class for plot_coherence_matrix (and plot_network)
    Example:
        cmd = 'plot_coherence_matrix.py ./INPUTS/ifgramStack.h5 --noverbose --figsize 9 3 --yx 216 310'
        obj = networkViewer(cmd)
        obj.configure()
        obj.plot()
    """
    def __init__(self, cmd=None, iargs=None):
        if cmd:
            iargs = cmd.split()[1:]
        self.cmd = cmd
        self.iargs = iargs

        # figure variables
        self.figname = 'Coherence matrix'
        self.fig_size = None
        self.fig = None
        self.ax_img = None
        self.ax_mat = None
        return

    def configure(self):
        inps = cmd_line_parse(self.iargs)
        # read network info
        inps = read_network_info(inps)
        # copy inps to self object
        for key, value in inps.__dict__.items():
            setattr(self, key, value)

        # auto figure size
        if not self.fig_size:
            ds_shape = readfile.read(self.img_file)[0].shape
            fig_size = pp.auto_figure_size(ds_shape, disp_cbar=True, ratio=0.7)
            self.fig_size = [fig_size[0]+fig_size[1], fig_size[1]]
            vprint('create figure in size of {} inches'.format(self.fig_size))
        return

    def plot(self):
        # Figure 1
        self.fig = plt.figure(self.figname, figsize=self.fig_size)
        # Axes 1 - Image
        self.ax_img = self.fig.add_axes([0.05, 0.1, 0.4, 0.8])
        view_cmd = self.view_cmd.format(self.img_file)
        d_img, atr, inps_img = view.prep_slice(view_cmd)
        if self.yx:
            inps_img.pts_yx = np.array(self.yx).reshape(-1, 2)
            inps_img.pts_marker = 'r^'
        inps_img.print_msg = self.print_msg
        self.ax_img = view.plot_slice(self.ax_img, d_img, atr, inps_img)[0]

        # coordinate info
        self.coord = ut.coordinate(atr)
        self.fig_coord = inps_img.fig_coord

        # Axes 2 - coherence matrix
        self.ax_mat = self.fig.add_axes([0.55, 0.125, 0.40, 0.75])
        if self.yx:
            self.plot_coherence_matrix4pixel(self.yx)

        # Link the canvas to the plots.
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.update_coherence_matrix)
        if self.disp_fig:
            plt.show()
        return

    def plot_coherence_matrix4pixel(self, yx):
        """Plot coherence matrix for one pixel
        Parameters: yx : list of 2 int
        """
        # read coherence
        box = (yx[1], yx[0], yx[1]+1, yx[0]+1)
        coh = readfile.read(self.ifgram_file, datasetName='coherence', box=box)[0]
        # prep metadata
        plotDict = {}
        plotDict['fig_title'] = 'Y = {}, X = {}'.format(yx[0], yx[1])
        plotDict['colormap'] = self.colormap
        plotDict['disp_legend'] = False
        # plot
        coh_mat = pp.plot_coherence_matrix(self.ax_mat,
                                           date12List=self.date12_list,
                                           cohList=coh.tolist(),
                                           date12List_drop=self.ex_date12_list,
                                           plot_dict=plotDict)[1]
        self.fig.canvas.draw()

        # status bar
        def format_coord(x, y):
            row, col = int(y+0.5), int(x+0.5)
            date12 = sorted([self.date_list[row], self.date_list[col]])
            date12 = ['{}-{}-{}'.format(i[0:4], i[4:6], i[6:8]) for i in date12]
            return 'x={}, y={}, v={:.3f}'.format(date12[0], date12[1], coh_mat[row, col])
        self.ax_mat.format_coord = format_coord
        # info
        vprint('-'*30)
        vprint('pixel: yx = {}'.format(yx))
        vprint('min/max coherence: {:.2f} / {:.2f}'.format(np.min(coh), np.max(coh)))
        return

    def update_coherence_matrix(self, event):
        if event.inaxes == self.ax_img:
            if self.fig_coord == 'geo':
                yx = [self.coord.lalo2yx(event.ydata, coord_type='lat'),
                      self.coord.lalo2yx(event.xdata, coord_type='lon')]
            else:
                yx = [int(event.ydata+0.5),
                      int(event.xdata+0.5)]
            self.plot_coherence_matrix4pixel(yx)
        return


##########################  Main Function  ##############################
def main(iargs=None):
    obj = networkViewer(cmd=iargs)
    obj.configure()
    obj.plot()
    obj.fig.canvas.mpl_disconnect(obj.cid)
    return


############################################################
if __name__ == '__main__':
    main()
