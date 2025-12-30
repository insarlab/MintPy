############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Nov 2018                           #
############################################################


import os

import matplotlib.pyplot as plt
import numpy as np

from mintpy import view
from mintpy.objects import ifgramStack
from mintpy.utils import plot as pp, readfile, utils as ut, ptime


###########################  Sub Function  #############################
def read_network_info(inps):
    """Read the network information"""
    ftype = readfile.read_attribute(inps.ifgram_file)['FILE_TYPE']
    if ftype != 'ifgramStack':
        raise ValueError(f'input file {inps.ifgram_file} is not ifgramStack: {ftype}')

    obj = ifgramStack(inps.ifgram_file)
    obj.open(print_msg=inps.print_msg)
    inps.date12_list = obj.get_date12_list(dropIfgram=False)
    date12_kept = obj.get_date12_list(dropIfgram=True)
    inps.ex_date12_list = sorted(list(set(inps.date12_list) - set(date12_kept)))
    inps.date_list = obj.get_date_list(dropIfgram=False)
    vprint(f'number of all     interferograms: {len(inps.date12_list)}')
    vprint(f'number of dropped interferograms: {len(inps.ex_date12_list)}')
    vprint(f'number of kept    interferograms: {len(inps.date12_list) - len(inps.ex_date12_list)}')
    vprint(f'number of acquisitions: {len(inps.date_list)}')

    if inps.lalo:
        if not inps.lookup_file:
            lookup_file = os.path.join(os.path.dirname(inps.ifgram_file), 'geometry*.h5')
            inps.lookup_file = ut.get_lookup_file(filePattern=lookup_file)
        coord = ut.coordinate(obj.metadata, lookup_file=inps.lookup_file)
        inps.yx = coord.geo2radar(inps.lalo[0], inps.lalo[1])[0:2]

    if not inps.yx:
        inps.yx = (obj.refY, obj.refX)
        vprint(f'plot initial coherence matrix at reference pixel: {inps.yx}')
    return inps


class coherenceMatrixViewer():
    """class for plot_coherence_matrix.

    Example:
        from mintpy.cli.plot_coherence_matrix import cmd_line_parse
        from mintpy.plot_coherence_matrix import coherenceMatrixViewer

        cmd = './inputs/ifgramStack.h5 --noverbose --figsize 9 3 --yx 216 310'
        inps = cmd_line_parse(cmd.split())
        obj = coherenceMatrixViewer(inps)
        obj.open()
        obj.plot()
    """

    def __init__(self, inps):
        # figure variables
        self.figname_img = 'Image'
        self.figsize_img = None
        self.fig_img = None
        self.ax_img = None
        self.cbar_img = None
        self.img = None

        self.figname_mat = 'Coherence Matrix'
        self.figsize_mat = None
        self.fig_mat = None
        self.ax_mat = None

        self.time_axis = getattr(inps, 'time_axis', False)

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

        # read network info
        self = read_network_info(self)

        # auto figure size
        if not self.figsize_img:
            ds_shape = readfile.read(self.img_file)[0].shape
            self.figsize_img = pp.auto_figure_size(ds_shape, disp_cbar=True, scale=0.7)
            vprint(f'create image figure in size of {self.figsize_img} inches')
        
        if not self.figsize_mat:
            num_ifg = len(self.date12_list)
            if num_ifg <= 50:
                self.figsize_mat = [6, 5]
            elif num_ifg <= 100:
                self.figsize_mat = [8, 6]
            else:
                self.figsize_mat = [10, 8]
            vprint(f'create matrix figure in size of {self.figsize_mat} inches')

        if not hasattr(self, 'cmap_name'):
            # Default colormap: use 'RdBu_truncate' for both timeaxis and normal mode (from CLI default)
            # This matches the CLI default value
            if self.time_axis:
                self.cmap_name = 'RdBu_truncate'
            else:
                self.cmap_name = 'viridis'
        if not hasattr(self, 'cmap_vlist'):
            self.cmap_vlist = [0.0, 1.0]
        self.colormap = pp.ColormapExt(self.cmap_name, vlist=self.cmap_vlist).colormap

        # read aux data
        # 1. temporal coherence value
        self.tcoh = None
        if self.tcoh_file:
            self.tcoh = readfile.read(self.tcoh_file)[0]
        # 2. minimum used coherence from template file
        self.min_coh_used = 0.0
        if self.template_file:
            template = readfile.read_template(self.template_file)
            template = ut.check_template_auto_value(template)
            if template['mintpy.networkInversion.maskDataset'] == 'coherence':
                self.min_coh_used = float(template['mintpy.networkInversion.maskThreshold'])
                vprint('Pixel-wised masking is applied in invert_network step')


    def plot(self):
        
        # Figure 1 - Image
        self.fig_img, self.ax_img = plt.subplots(num=self.figname_img, figsize=self.figsize_img)
        self.plot_init_image()
        
        # Figure 2 - Coherence Matrix
        self.fig_mat, self.ax_mat = plt.subplots(num=self.figname_mat, figsize=self.figsize_mat)
        if all(i is not None for i in self.yx):
            self.plot_coherence_matrix4pixel(self.yx)

        # Link the canvas to the plots.
        self.cid_img = self.fig_img.canvas.mpl_connect('button_press_event', self.update_coherence_matrix)
        self.cid_mat = self.fig_mat.canvas.mpl_connect('button_press_event', self.update_coherence_matrix)
        
        if self.disp_fig:
            plt.show()
        return

    def plot_init_image(self):
        """Plot the initial image"""

        view_cmd = self.view_cmd.format(self.img_file)
        d_img, atr, view_inps = view.prep_slice(view_cmd)
        self.coord = ut.coordinate(atr)

        if all(i is not None for i in self.yx):
            view_inps.pts_marker = 'r^'
            view_inps.pts_yx = np.array(self.yx).reshape(-1, 2)

            # point yx --> lalo for geocoded product
            if 'Y_FIRST' in atr.keys():
                view_inps.pts_lalo = np.array(
                    self.coord.radar2geo(
                        self.yx[0],
                        self.yx[1],
                    )[0:2],
                ).reshape(-1,2)

        view_inps.print_msg = self.print_msg
        self.ax_img = view.plot_slice(self.ax_img, d_img, atr, view_inps)[0]
        self.fig_coord = view_inps.fig_coord
        

        self.fig_img.canvas.manager.set_window_title(self.figname_img)
        self.fig_img.tight_layout()

    def plot_coherence_matrix4pixel_time_axis(self, yx):
        """Plot coherence matrix with continuous time axis for one pixel
        Parameters: yx : list of 2 int
        """
        self.ax_mat.cla()

        # read coherence
        box = (yx[1], yx[0], yx[1]+1, yx[0]+1)
        coh = readfile.read(self.ifgram_file, datasetName='coherence', box=box)[0]

        # ex_date for pixel-wise masking during network inversion
        ex_date12_list = self.ex_date12_list[:]   #local copy
        if self.min_coh_used > 0.:
            ex_date12_list += np.array(self.date12_list)[coh < self.min_coh_used].tolist()
            ex_date12_list = sorted(list(set(ex_date12_list)))

        # prep metadata
        plotDict = {}
        plotDict['fig_title'] = f'Y = {yx[0]}, X = {yx[1]}'
        # display temporal coherence value of the pixel
        if self.tcoh_file:
            tcoh = self.tcoh[yx[0], yx[1]]
            plotDict['fig_title'] += f', tcoh = {tcoh:.2f}'
        plotDict['colormap'] = self.colormap
        # cmap_vlist is [start, jump, end] for truncated colormap, but vlim needs [vmin, vmax]
        if len(self.cmap_vlist) >= 2:
            plotDict['vlim'] = [self.cmap_vlist[0], self.cmap_vlist[-1]]
        else:
            plotDict['vlim'] = [0.0, 1.0]
        plotDict['cbar_label'] = 'Coherence'
        plotDict['disp_legend'] = False

        # plot using the utility function
        Z, mesh = pp.plot_coherence_matrix_time_axis(
            self.ax_mat,
            date12List=self.date12_list,
            cohList=coh.tolist(),
            date12List_drop=ex_date12_list,
            p_dict=plotDict,
        )[1:3]

        # Info
        msg = f'pixel in yx = {tuple(yx)}, '
        msg += f'min/max spatial coherence: {np.nanmin(coh):.2f} / {np.nanmax(coh):.2f}, '
        if self.tcoh_file:
            tcoh = self.tcoh[yx[0], yx[1]]
            msg += f'temporal coherence: {tcoh:.2f}'
        vprint(msg)

        self.ax_mat.annotate('ifgrams\navailable', xy=(0.05, 0.05), xycoords='axes fraction', fontsize=12)
        self.ax_mat.annotate('ifgrams\nused', ha='right', xy=(0.95, 0.85), xycoords='axes fraction', fontsize=12)

        self.fig_mat.canvas.manager.set_window_title(self.figname_mat)
        self.fig_mat.tight_layout()

        # Update figure
        self.fig_mat.canvas.draw_idle()
        self.fig_mat.canvas.flush_events()
        return

    def plot_coherence_matrix4pixel(self, yx):
        """Plot coherence matrix for one pixel
        Parameters: yx : list of 2 int
        """
        # Use time axis mode if enabled
        if self.time_axis:
            return self.plot_coherence_matrix4pixel_time_axis(yx)
        
        self.ax_mat.cla()

        # read coherence
        box = (yx[1], yx[0], yx[1]+1, yx[0]+1)
        coh = readfile.read(self.ifgram_file, datasetName='coherence', box=box)[0]

        # ex_date for pixel-wise masking during network inversion
        ex_date12_list = self.ex_date12_list[:]   #local copy
        if self.min_coh_used > 0.:
            ex_date12_list += np.array(self.date12_list)[coh < self.min_coh_used].tolist()
            ex_date12_list = sorted(list(set(ex_date12_list)))

        # prep metadata
        plotDict = {}
        plotDict['fig_title'] = f'Y = {yx[0]}, X = {yx[1]}'
        # display temporal coherence value of the pixel
        if self.tcoh_file:
            tcoh = self.tcoh[yx[0], yx[1]]
            plotDict['fig_title'] += f', tcoh = {tcoh:.2f}'
        plotDict['colormap'] = self.colormap
        plotDict['cmap_vlist'] = self.cmap_vlist
        plotDict['disp_legend'] = False

        # plot
        coh_mat = pp.plot_coherence_matrix(
            self.ax_mat,
            date12List=self.date12_list,
            cohList=coh.tolist(),
            date12List_drop=ex_date12_list,
            p_dict=plotDict,
        )[1]

        self.ax_mat.annotate('ifgrams\navailable', xy=(0.05, 0.05), xycoords='axes fraction', fontsize=12)
        self.ax_mat.annotate('ifgrams\nused', ha='right', xy=(0.95, 0.85), xycoords='axes fraction', fontsize=12)

        # status bar
        def format_coord(x, y):
            row, col = int(y+0.5), int(x+0.5)
            date12 = sorted([self.date_list[row], self.date_list[col]])
            date12 = [f'{i[0:4]}-{i[4:6]}-{i[6:8]}' for i in date12]
            return f'x={date12[0]}, y={date12[1]}, v={coh_mat[row, col]:.3f}'
        self.ax_mat.format_coord = format_coord

        # info
        msg = f'pixel in yx = {tuple(yx)}, '
        msg += f'min/max spatial coherence: {np.min(coh):.2f} / {np.max(coh):.2f}, '
        if self.tcoh_file:
            msg += f'temporal coherence: {tcoh:.2f}'
        vprint(msg)

        self.fig_mat.canvas.manager.set_window_title(self.figname_mat)
        self.fig_mat.tight_layout()

        # update figure
        self.fig_mat.canvas.draw_idle()
        self.fig_mat.canvas.flush_events()
        return

    def update_coherence_matrix(self, event):
        """Update coherence matrix when clicking on either window"""
        if event.inaxes == self.ax_img:
            if self.fig_coord == 'geo':
                yx = self.coord.lalo2yx(event.ydata, event.xdata)
            else:
                yx = [int(event.ydata+0.5),
                      int(event.xdata+0.5)]
            self.plot_coherence_matrix4pixel(yx)
            
            self.update_image_marker(yx)
        elif event.inaxes == self.ax_mat:
            pass

    def update_image_marker(self, yx):
        """Update the marker point in the image window"""
        if hasattr(self, 'pts_yx'):
            for artist in self.ax_img.get_children():
                if hasattr(artist, 'get_marker') and artist.get_marker() == '^':
                    artist.remove()
            
            self.ax_img.plot(yx[1], yx[0], 'r^', markersize=10, markeredgecolor='black')
            self.fig_img.canvas.draw_idle()
