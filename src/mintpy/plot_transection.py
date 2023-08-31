############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Heresh Fattahi, 2013               #
############################################################


import os

from matplotlib import pyplot as plt, ticker

from mintpy import view
from mintpy.utils import plot as pp, readfile, utils as ut


#####################################################################
class transectionViewer():
    """class for plot_transection.

    Example:
        from mintpy.cli.plot_transection import cmd_line_parse, get_view_cmd
        from mintpy.plot_transection import transectionViewer

        cmd = 'velocity.h5 --noverbose --start-yx 10 10 --end-yx 200 300'
        inps = cmd_line_parse(cmd.split())
        view_cmd = get_view_cmd(cmd.split())
        obj = transectionViewer(inps, view_cmd)
        obj.open()
        obj.plot()
        obj.fig.canvas.mpl_disconnect(obj.cid)
    """

    def __init__(self, inps, view_cmd):

        # figure variables
        self.figname = 'Transection'
        self.fig = None
        self.ax_img = None
        self.ax_txn = None

        self.img = None
        self.line_ann = None
        self.pts_idx = 0

        self.view_cmd = view_cmd

        # copy inps to self object
        for key, value in inps.__dict__.items():
            setattr(self, key, value)


    def open(self):
        global vprint
        vprint = print if self.print_msg else lambda *args, **kwargs: None

        # print command line
        if self.argv is not None:
            print(f'{os.path.basename(__file__)} ' + ' '.join(self.argv))

        # matplotlib backend setting
        if not self.disp_fig:
            plt.switch_backend('Agg')

        # copy inps from view.py to self object
        self.data_img, self.atr, view_inps = view.prep_slice(self.view_cmd)
        for key, value in view_inps.__dict__.items():
            # do not update the following setting from view.py
            if key not in ['file', 'dset', 'fig_size']:
                setattr(self, key, value)

        self.offset *= self.disp_scale  # due to unit change from view.prep_slice()

        # auto figure size
        if not self.fig_size:
            length, width = int(self.atr['LENGTH']), int(self.atr['WIDTH'])
            fig_size = pp.auto_figure_size((length, width), disp_cbar=True)
            self.fig_size = [fig_size[0] + fig_size[1], fig_size[1]]


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
            outfile = f'{self.outfile_base}.pdf'
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

                self.draw_line(
                    start_yx=self.start_yx,
                    end_yx=self.end_yx,
                )

                self.draw_transection(
                    start_yx=self.start_yx,
                    end_yx=self.end_yx,
                    start_lalo=self.start_lalo,
                    end_lalo=self.end_lalo,
                )

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
                txn = ut.transect_lalo(
                    self.data_list[i],
                    self.atr_list[i],
                    start_lalo, end_lalo,
                    interpolation=self.interpolation)
            else:
                txn = ut.transect_yx(
                    self.data_list[i],
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
            # update distance values by excluding the commonly masked out pixels in the beginning
            self.ax_txn.scatter(
                x=(txn['distance'] - min_dist) * dist_scale,
                y=txn['value'] - self.offset[i],
                c=f'C{i}',
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
