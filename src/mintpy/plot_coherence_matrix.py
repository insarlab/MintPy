############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Nov 2018                           #
############################################################


import os
from datetime import datetime, timedelta

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
            # Default colormap: use 'RdBu' for timeaxis mode, 'viridis' for normal mode
            if self.time_axis:
                self.cmap_name = 'RdBu'
            else:
                self.cmap_name = 'viridis'
        elif self.time_axis and self.cmap_name == 'RdBu_truncate':
            # If using timeaxis mode and user didn't specify colormap (using default 'RdBu_truncate'),
            # change to 'RdBu' for timeaxis mode
            self.cmap_name = 'RdBu'
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

        # Convert date strings to datetime objects using MintPy utilities
        # Normalize date format to YYYYMMDD
        date_list_normalized = ptime.yyyymmdd(self.date_list)
        date_objs = {}
        for date_str in date_list_normalized:
            try:
                date_objs[date_str] = datetime.strptime(date_str, '%Y%m%d')
            except ValueError:
                # Fallback: try with YYMMDD format
                if len(date_str) == 6:
                    date_objs[date_str] = datetime.strptime('20' + date_str, '%Y%m%d')
                else:
                    raise

        # Normalize date12 format
        date12_list_normalized = ptime.yyyymmdd_date12(self.date12_list)
        ex_date12_list_normalized = ptime.yyyymmdd_date12(ex_date12_list) if ex_date12_list else []
        
        # Create coherence dictionary
        coh_dict = {}
        excluded_pairs = set()
        for date12, coh_val in zip(date12_list_normalized, coh):
            date1_str, date2_str = date12.split('_')
            # Ensure we have datetime objects
            if date1_str not in date_objs:
                date1_str = ptime.yyyymmdd([date1_str])[0]
            if date2_str not in date_objs:
                date2_str = ptime.yyyymmdd([date2_str])[0]
            
            date1 = date_objs.get(date1_str)
            date2 = date_objs.get(date2_str)
            
            if date1 is None:
                date1 = datetime.strptime(date1_str, '%Y%m%d')
                date_objs[date1_str] = date1
            if date2 is None:
                date2 = datetime.strptime(date2_str, '%Y%m%d')
                date_objs[date2_str] = date2
            
            # Store as tuple (date1, date2) where date1 <= date2 for consistency
            pair = (min(date1, date2), max(date1, date2))
            coh_dict[pair] = float(coh_val)
            
            # Mark excluded pairs
            if date12 in ex_date12_list_normalized:
                excluded_pairs.add(pair)

        # Get all unique dates
        all_dates = set()
        for d1, d2 in coh_dict.keys():
            all_dates.add(d1)
            all_dates.add(d2)
        date_list = sorted(list(all_dates))

        # Create continuous time grid (based on actual data points)
        grid_points = [date_list[0]]  # starting point
        for i in range(len(date_list)-1):
            mid_point = date_list[i] + (date_list[i+1] - date_list[i])/2
            grid_points.append(mid_point)
        grid_points.append(date_list[-1])  # ending point

        # Convert to days for plotting
        base_date = min(date_list)
        days_grid = [(d - base_date).days for d in grid_points]

        # Create meshgrid
        X, Y = np.meshgrid(days_grid, days_grid)

        # Create value matrix, initialized with NaN (will display as white)
        Z = np.full((len(grid_points)-1, len(grid_points)-1), np.nan)

        # Create a mapping from date to grid index for fast lookup
        # Each date in date_list corresponds to a grid point
        date_to_grid_idx = {}
        for date in date_list:
            # Find the grid index where this date falls
            # Grid cells are defined by [grid_points[i], grid_points[i+1]]
            for grid_idx in range(len(grid_points)-1):
                # Check if date falls within this grid cell
                # Note: grid_points includes midpoints, so we need to check carefully
                if grid_idx == 0:
                    # First cell: from grid_points[0] to grid_points[1]
                    if grid_points[0] <= date <= grid_points[1]:
                        date_to_grid_idx[date] = grid_idx
                        break
                elif grid_idx == len(grid_points) - 2:
                    # Last cell: from grid_points[-2] to grid_points[-1]
                    if grid_points[grid_idx] < date <= grid_points[grid_idx+1]:
                        date_to_grid_idx[date] = grid_idx
                        break
                else:
                    # Middle cells: from grid_points[i] to grid_points[i+1]
                    if grid_points[grid_idx] < date <= grid_points[grid_idx+1]:
                        date_to_grid_idx[date] = grid_idx
                        break

        # Fill value matrix directly from coherence dictionary
        # This is much faster: O(n_ifgrams) instead of O(n_grid^2 * n_ifgrams)
        for (d1, d2), cor in coh_dict.items():
            # Find grid indices for both dates
            idx1 = date_to_grid_idx.get(d1)
            idx2 = date_to_grid_idx.get(d2)
            
            # Only fill if both dates are in valid grid cells
            if idx1 is not None and idx2 is not None:
                # Fill both upper and lower triangle (symmetric)
                Z[idx1, idx2] = cor
                if idx1 != idx2:
                    Z[idx2, idx1] = cor
                # Otherwise, Z[i, j] remains NaN (blank)

        # Create diagonal matrix for black diagonal cells (like in original implementation)
        # The diagonal should be black to distinguish from un-selected interferograms
        # In timeaxis mode, diagonal corresponds to same date pairs (idx1 == idx2)
        diag_Z = np.full_like(Z, np.nan)
        num_cells = len(grid_points) - 1
        for i in range(min(num_cells, Z.shape[0], Z.shape[1])):
            diag_Z[i, i] = 1.0
        
        # Plot diagonal as black first (using gray_r colormap, where 1.0 = black)
        if np.any(~np.isnan(diag_Z)):
            self.ax_mat.pcolormesh(X, Y, diag_Z,
                                 cmap='gray_r',
                                 vmin=0.0,
                                 vmax=1.0,
                                 shading='auto',
                                 zorder=1)
        
        # Plot using pcolormesh for coherence values
        cmap = self.colormap.copy()
        cmap.set_bad('white')  # NaN values will be white
        
        mesh = self.ax_mat.pcolormesh(X, Y, Z,
                                     cmap=cmap,
                                     vmin=self.cmap_vlist[0],
                                     vmax=self.cmap_vlist[-1],
                                     shading='auto',
                                     zorder=0)

        # Generate month ticks
        min_date = min(date_list)
        max_date = max(date_list)

        # If min_date is not the first day of month, start from next month 1st
        if min_date.day > 1:
            if min_date.month == 12:
                current_date = min_date.replace(year=min_date.year+1, month=1, day=1)
            else:
                current_date = min_date.replace(month=min_date.month+1, day=1)
        else:
            current_date = min_date.replace(day=1)

        tick_dates = []
        while current_date <= max_date:
            tick_dates.append(current_date)
            # Get next month
            if current_date.month == 12:
                current_date = current_date.replace(year=current_date.year+1, month=1)
            else:
                current_date = current_date.replace(month=current_date.month+1)

        # Calculate tick positions (at month start)
        tick_positions = [(d - base_date).days for d in tick_dates]

        # Calculate label positions (at middle of adjacent ticks)
        label_positions = []
        month_labels = []
        is_january = []

        for i in range(len(tick_dates)-1):
            # Calculate middle point of adjacent ticks
            mid_point = (tick_positions[i] + tick_positions[i+1]) / 2

            # Only add label for odd months
            if tick_dates[i].month % 2 == 1:
                label_positions.append(mid_point)
                month_labels.append(tick_dates[i].strftime('%-m'))
                is_january.append(tick_dates[i].month == 1)

        # Separate January and other month tick positions
        major_ticks = [pos for pos, date in zip(tick_positions, tick_dates) if date.month == 1]
        minor_ticks = [pos for pos, date in zip(tick_positions, tick_dates) if date.month != 1]

        # Separate January and other odd month label positions
        major_label_pos = [pos for pos, is_jan in zip(label_positions, is_january) if is_jan]
        minor_label_pos = [pos for pos, is_jan in zip(label_positions, is_january) if not is_jan]

        # Separate labels
        major_labels = [label for label, is_jan in zip(month_labels, is_january) if is_jan]
        minor_labels = [label for label, is_jan in zip(month_labels, is_january) if not is_jan]

        # Set tick positions (no labels)
        self.ax_mat.set_xticks(major_ticks)  # major ticks (January)
        self.ax_mat.set_xticks(minor_ticks, minor=True)  # minor ticks (other months)
        self.ax_mat.set_yticks(major_ticks)
        self.ax_mat.set_yticks(minor_ticks, minor=True)

        # Set empty labels (we'll add labels separately with text)
        self.ax_mat.set_xticklabels([''] * len(major_ticks))
        self.ax_mat.set_xticklabels([''] * len(minor_ticks), minor=True)
        self.ax_mat.set_yticklabels([''] * len(major_ticks))
        self.ax_mat.set_yticklabels([''] * len(minor_ticks), minor=True)

        # Normal (non-rotated) label positions
        offset = (self.ax_mat.get_ylim()[1] - self.ax_mat.get_ylim()[0]) * 0.03
        year_offset = offset * 2.2  # Year labels below month labels

        # Add month labels
        for pos, label in zip(major_label_pos, major_labels):
            self.ax_mat.text(pos, self.ax_mat.get_ylim()[1] + offset, label,
                            horizontalalignment='center', verticalalignment='top', fontsize=10)
            self.ax_mat.text(self.ax_mat.get_xlim()[0] - offset, pos, label,
                            horizontalalignment='right', verticalalignment='center', fontsize=10)

        for pos, label in zip(minor_label_pos, minor_labels):
            self.ax_mat.text(pos, self.ax_mat.get_ylim()[1] + offset, label,
                            horizontalalignment='center', verticalalignment='top', fontsize=10)
            self.ax_mat.text(self.ax_mat.get_xlim()[0] - offset, pos, label,
                            horizontalalignment='right', verticalalignment='center', fontsize=10)

        # Set tick line style
        self.ax_mat.tick_params(which='major', direction='out', length=6, width=1.1,
                               bottom=True, top=True, left=True, right=True)
        self.ax_mat.tick_params(which='minor', direction='out', length=3, width=1,
                               bottom=True, top=True, left=True, right=True)

        # Add year labels (below month labels)
        years = []
        year_positions = []
        prev_year = None
        for i, d in enumerate(tick_dates):
            if prev_year != d.year:
                years.append(str(d.year))
                year_positions.append(tick_positions[i])
                prev_year = d.year

        # Display year labels
        # Normal (non-rotated) year label positions
        for pos, year in zip(year_positions, years):
            # X-axis: year labels at bottom (below month labels)
            self.ax_mat.text(pos, self.ax_mat.get_ylim()[1] + year_offset, year,
                            horizontalalignment='center', verticalalignment='top', fontsize=10)
            # Y-axis: year labels at left (below month labels)
            self.ax_mat.text(self.ax_mat.get_xlim()[0] - year_offset, pos, year,
                            horizontalalignment='right', verticalalignment='center', fontsize=10)
        # Invert Y axis
        self.ax_mat.invert_yaxis()

        # Add colorbar
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(self.ax_mat)
        cax = divider.append_axes("right", "3%", pad="3%")
        cbar = plt.colorbar(mesh, cax=cax)
        cbar.set_label('Coherence', fontsize=12)

        # Set title
        title = f'Y = {yx[0]}, X = {yx[1]}'
        if self.tcoh_file:
            tcoh = self.tcoh[yx[0], yx[1]]
            title += f', tcoh = {tcoh:.2f}'
        self.ax_mat.set_title(title, fontsize=12)

        # Status bar - format coordinate display
        def format_coord(x, y):
            x_idx = np.argmin(np.abs(np.array(days_grid) - x))
            y_idx = np.argmin(np.abs(np.array(days_grid) - y))
            
            # Clamp indices to valid range
            x_idx = min(max(0, x_idx), len(grid_points) - 1)
            y_idx = min(max(0, y_idx), len(grid_points) - 1)
            
            if x_idx < len(grid_points) and y_idx < len(grid_points):
                date1 = grid_points[x_idx]
                date2 = grid_points[y_idx]
                date1_str = date1.strftime('%Y-%m-%d')
                date2_str = date2.strftime('%Y-%m-%d')
                
                # Find coherence value (Z has shape (len(grid_points)-1, len(grid_points)-1))
                coh_val = np.nan
                if x_idx < len(grid_points) - 1 and y_idx < len(grid_points) - 1:
                    coh_val = Z[y_idx, x_idx]
                
                if not np.isnan(coh_val):
                    return f'x={date1_str}, y={date2_str}, v={coh_val:.3f}'
                else:
                    return f'x={date1_str}, y={date2_str}, v=NaN'
            return ''
        
        self.ax_mat.format_coord = format_coord

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
