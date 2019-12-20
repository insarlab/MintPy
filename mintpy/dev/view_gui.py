#! /usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Joshua Zahner, 2017                              #
############################################################


from tkinter import *

import h5py
import numpy
import matplotlib
import argparse

matplotlib.use('TkAgg')
from tkinter import filedialog
from mintpy.utils import readfile, utils as ut
from mintpy import info, view

canvas, frame, h5_file, h5_file_short, pick_h5_file_button, mask_file, mask_short, \
pick_mask_file_button, starting_upper_lim, starting_lower_lim, y_lim_upper, y_lim_upper_slider, y_lim_lower, y_lim_lower_slider, unit, \
colormap, projection, lr_flip, ud_flip, wrap, opposite, transparency, show_info, dem_file, dem_short, \
pick_dem_file_button, shading, countours, countour_smoothing, countour_step, subset_x_from, subset_x_to, subset_y_from, \
subset_y_to, subset_lat_from, subset_lat_to, subset_lon_from, subset_lon_to, ref_x, ref_y, ref_lat, ref_lon, ref_color, \
ref_sym, ref_date, font_size, plot_dpi, row_num, col_num, axis_show, cbar_show, title_show, tick_show, title_in, title, \
fig_size_width, fig_size_height, fig_ext, fig_num, fig_w_space, fig_h_space, coords, coastline, resolution, lalo_label, \
lalo_step, scalebar_distance, scalebar_lat, scalebar_lon, show_scalebar, save, output_file, epoch_option_menu, epoch, \
excludes_list_box, ref_date_option_menu, show_ref \
    = None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, \
      None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, \
      None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, \
      None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, \
      None, None

file_base = "/"

use_default = None
current_slider_scale = None
file_info_window = None
file_info = None
inps = None
settings_file = None
settings_file_button = None;

epoch_list = ["All"]
ref_dates_list = [""]

colormaps = ['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap',
             'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r',
             'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2',
             'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r',
             'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn',
             'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral',
             'Spectral_r', 'Vega10', 'Vega10_r', 'Vega20', 'Vega20_r', 'Vega20b', 'Vega20b_r', 'Vega20c', 'Vega20c_r', 'Wistia',
             'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot',
             'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r',
             'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r',
             'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar',
             'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot',
             'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r',
             'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink',
             'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spectral',
             'spectral_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b',
             'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'viridis', 'viridis_r', 'winter', 'winter_r']

projections = ["cea", "mbtfpq", "aeqd", "sinu", "poly", "moerc", "gnom", "moll", "lcc", "tmerc", "nplaea", "gall",
               "npaeqd", "mill", "merc", "stere", "eqdc", "rotpole", "cyl", "npstere", "spstere", "hammer", "geos",
               "nsper", "eck4", "aea", "kav7", "spaeqd", "ortho", "cass", "vandg", "laea", "splaea", "robin",
               " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " "]

unit_options = ["cm", "dm", "m",  "km", "", "cm/yr", "dm/yr", "m/yr", "km/yr"]

attributes = []
update_in_progress = False

settings_to_variable = None


def parse_settings():

    global settings_file

    with open(settings_file.get(), 'r') as the_settings_file:
        for line in the_settings_file:
            if line is not "" and "=" in line:
                parts = line.split("=")
                key = parts[0].rstrip().strip("\t")
                value = parts[1].strip(" \t\t").strip("\n").rstrip()

                if key in settings_to_variable.keys():

                    if value.lower() == 'false':
                        value = 0
                    elif value.lower() == 'true':
                        value = 1

                    if key in ['input.file']:
                        inps.file = value
                        on_file_selection(inps.file)

                    settings_to_variable[key].set(value)


def write_settings_file():

    output = ""

    for setting in sorted(settings_to_variable.keys()):

        value = str(settings_to_variable[setting].get())

        if setting in [lr_flip, ud_flip, wrap, opposite, shading, countours, axis_show, cbar_show, title_show,
                       tick_show, title_in, coastline, lalo_label, show_scalebar] and value is 0:
            value = "False"
        elif setting in [lr_flip, ud_flip, wrap, opposite, shading, countours, axis_show, cbar_show, title_show,
                         tick_show, title_in, coastline, lalo_label, show_scalebar] and value is 1:
            value = "True"

        if value != "":
            output += "{0:35} = {1:10}\n".format(str(setting), str(value))
            # "{0:35} = {1:10}\n".format(name, value)

    f = filedialog.asksaveasfile(mode='w', defaultextension=".txt")
    if f is None:  # asksaveasfile return `None` if dialog closed with "cancel".
        return

    f.write(output)
    f.close()


def cmdLineParse(argv):
    global inps
    parser = argparse.ArgumentParser(description='Display InSAR Product',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=None)

    ##### Input
    infile = parser.add_argument_group('Input File', 'File/Dataset to display')
    infile.add_argument('--file', dest='file', metavar='FILE', help='file for display')
    infile.add_argument('--settings', dest='settings', metavar='FILE', help='settings file to use in setup')

    inps = parser.parse_args(argv)
    return inps


def pick_file():

    global file_info, file_base

    if len(epoch_list) > 0:
        epoch_option_menu['menu'].delete(0, "end")
    if len(ref_dates_list) > 0:
        ref_date_option_menu['menu'].delete(0, "end")

    if h5_file.get() == "":
        filename = filedialog.askopenfilename(initialdir=file_base, title="Select file",
                                              filetypes=(("jpeg files", "*.h5"), ("all files", "*.*"), ("more files", "*.he5")))
        frame.filename = filename

        on_file_selection(filename)

        #return frame.filename
    else:
        h5_file.set("")
        h5_file_short.set("No File Selected")
        pick_h5_file_button.config(text="Select .h5 File")
        epoch_option_menu['menu'].delete(0, 'end')
        ref_date_option_menu['menu'].delete(0, 'end')
        excludes_list_box.delete(0, 'end')
        y_lim_upper.set(0)
        y_lim_lower.set(0)
        file_info = None


def pick_mask():
    if mask_file.get() == "":
        filename = filedialog.askopenfilename(initialdir=file_base, title="Select file",
                                              filetypes=(("jpeg files", "*.h5"), ("all files", "*.*")))
        frame.filename = filename
        if filename != "":
            mask_file.set(frame.filename)
        #return frame.filename
    else:
        mask_file.set("")

    set_mask_short()


def pick_dem():
    if dem_file.get() == "":
        filename = filedialog.askopenfilename(initialdir=file_base, title="Select file",
                                              filetypes=(("jpeg files", "*.h5"), ("all files", "*.*")))
        frame.filename = filename
        if filename != "":
            dem_file.set(frame.filename)
    else:
        dem_file.set("")
        dem_short.set("No File Selected")

    set_dem_short()


def pick_settings():

    global settings_file_button

    if settings_file.get() == "":
        filename = filedialog.askopenfilename(initialdir=file_base, title="Select file",
                                              filetypes=(("TEXT files", "*.txt"), ("all files", "*.*")))
        frame.filename = filename
        if filename != "":
            settings_file.set(frame.filename)
            settings_file_button.config(text=str(frame.filename).split("/")[-1])
            parse_settings()
        #return frame.filename
    else:
        print("BAD")
        settings_file.set("")
        settings_file_button.config(text="Select a Settings File")



def on_file_selection(file):
    global file_base, file_info
    h5_file.set(file)
    inps.file = file
    parts = h5_file.get().split("/")
    file_base = parts[0]
    h5_file_short.set(parts[-1])

    if h5_file.get() != "":
        pick_h5_file_button.config(text="Cancel")
    else:
        h5_file_short.set("No File Selected")
        pick_h5_file_button.config(text="Select .h5 File")
        return

    file_info = info.hdf5_structure_string(h5_file.get())

    set_variables_from_attributes()


def on_configure(event):
    canvas.configure(scrollregion=canvas.bbox('all'))


def update_sliders(unit, setValues=True):
    global current_slider_scale

    scale = 1.0
    new_max = starting_upper_lim
    new_min = starting_lower_lim

    if unit == "m":
        scale = 1.0
    elif unit == "cm":
        scale = 100.0
    elif unit == "mm":
        scale = 1000.0
    elif unit == "dm":
        scale = 10
    elif unit == "km":
        scale = 0.001

    current_slider_scale = scale

    y_lim_upper_slider.configure(to_= new_max*scale)
    y_lim_upper_slider.configure(from_= -1*float(new_max*scale)/4)
    y_lim_lower_slider.configure(to_= float(new_max*scale)/4)
    y_lim_lower_slider.configure(from_= new_min*scale)

    if setValues:
        current_y_lim_upper = y_lim_upper.get()
        current_y_lim_lower = y_lim_lower.get()

        y_lim_upper.set(current_y_lim_upper*scale)
        y_lim_lower.set(current_y_lim_lower*scale)


def scale_sliders(unit):
    global current_slider_scale

    new_max_slider_val = y_lim_upper.get() / current_slider_scale
    new_min_slider_val = y_lim_lower.get() / current_slider_scale

    print(new_max_slider_val)

    if unit in ["m", "m/yr"]:
        current_slider_scale = 1.0
    elif unit in ["cm", "cm/yr"]:
        current_slider_scale = 100.0
    elif unit in ["mm", "mm/yr"]:
        current_slider_scale = 1000.0
    elif unit in ["dm", "dm/yr"]:
        current_slider_scale = 10
    elif unit in ["km", "km/yr"]:
        current_slider_scale = 0.001

    new_max_slider_val *= current_slider_scale
    new_min_slider_val *= current_slider_scale

    update_sliders(unit, setValues=False)

    y_lim_upper_slider.set(new_max_slider_val)
    y_lim_lower_slider.set(new_min_slider_val)



def show_file_info():
    global file_info_window, file_info

    file_info_window = Tk()
    file_info_window.minsize(width=350, height=550)
    file_info_window.maxsize(height=550)
    file_info_window.resizable(width=True, height=False)

    text_box = Text(file_info_window, wrap=NONE)
    text_box.insert(END, file_info)
    text_box.config(height=550)
    text_box.config(state=DISABLED)

    text_box.pack(fill=X)

    def close():
        global file_info_window

        file_info_window.destroy()
        file_info_window = None


    file_info_window.protocol("WM_DELETE_WINDOW", close)


def in_range(parameter, value):

    min_val, max_val = 0, 0

    if parameter is 'X':
        min_val = 0
        max_val = attributes['XMAX']
    elif parameter is 'Y':
        min_val = 0
        max_val = attributes['YMAX']
    elif parameter is 'LAT':
        min_val = compute_lalo(0, 0)[0]
        max_val = compute_lalo(0, attributes['YMAX'])[0]
    elif parameter is 'LON':
        min_val = compute_lalo(0, 0)[1]
        max_val = compute_lalo(attributes['XMAX'], 0)[1]
    else:
        raise ValueError('Illegal Parameter. Accepted parameters are: X, Y, LAT, LON')

    return int(min_val) <= int(value.strip(" ")) <= int(max_val)


def show_plot():

    global file_info_window, attributes

    options = [inps.file]

    print(inps.file)
    print(attributes)

    if epoch.get() != "All" and epoch.get() != 'None':
        options.append(epoch.get())

    options += ["--alpha", str(transparency.get())]

    if mask_file.get() != "":
        options.append("--mask")
        options.append(mask_file.get())

    excludes = [str(excludes_list_box.get(idx)) for idx in excludes_list_box.curselection()]
    if len(excludes) > 0:
        options.append("--exclude")
        for ex in excludes:
            options.append(str(ex))

    if use_default.get() != 1:
        options.append("-m")
        options.append(str(y_lim_lower.get()))
        options.append("-M")
        options.append(str(y_lim_upper.get()))

    if unit.get() != "":
        options.append("-u")
        options.append(unit.get())
    if colormap.get() != "":
        options.append("-c")
        options.append(colormap.get())
    if projection.get() != "":
        options.append("--projection")
        options.append(projection.get())
    if lr_flip.get() == 1:
        options.append("--flip-lr")
    if ud_flip.get() == 1:
        options.append("--flip-ud")
    if wrap.get() == 1:
        options.append("--wrap")
    if opposite.get() == 1:
        options.append("--opposite")

    if dem_file.get() != "":
        options.append("--dem")
        options.append(dem_file.get())
    if shading.get() == 0:
        options.append("--dem-noshade")
    if countours.get() == 0:
        options.append("--dem-nocontour")
    if countour_smoothing.get() != "":
        options.append("--contour-smooth")
        options.append(countour_smoothing.get())
    if countour_step.get() != "":
        options.append("--contour-step")
        options.append(countour_step.get())

    if subset_x_from.get() != "" and subset_x_to.get() != "" and in_range('X', subset_x_from.get()) and in_range('X', subset_x_to.get()):
        options.append("-x")
        options.append(subset_x_from.get())
        options.append(subset_x_to.get())
    elif len(attributes) > 0 and attributes['XMAX']:
        options.append("-x")
        options.append(0)
        options.append(attributes['XMAX'])

    if subset_y_from.get() != "" and subset_y_to.get() != "" and in_range('Y', subset_y_from.get()) and in_range('Y', subset_y_to.get()):
        options.append("-y")
        options.append(subset_y_from.get())
        options.append(subset_y_to.get())
    elif len(attributes) > 0 and attributes['YMAX']:
        options.append("-y")
        options.append(0)
        options.append(attributes['YMAX'])

    if ref_x.get() != "" and ref_y.get() != "":
        options.append("--ref-yx")
        options.append(ref_y.get())
        options.append(ref_x.get())
    if ref_lat.get() != "" and ref_lon.get() != "":
        options.append("--ref-lalo")
        options.append(ref_lat.get())
        options.append(ref_lon.get())

    if show_ref == 0:
        options.append("--noreference")
    if ref_color.get() != "":
        options.append("--ref-color")
        options.append(ref_color.get())
    if ref_sym.get() != "":
        options.append("--ref-symbol")
        options.append(ref_sym.get())
    if ref_date.get() != "":
        options.append("--ref-date")
        options.append(ref_date.get())

    if font_size.get() != "":
        options.append("-s")
        options.append(font_size.get())
    if plot_dpi.get() != "":
        options.append("--dpi")
        options.append(plot_dpi.get())
    if row_num.get() != "":
        options.append("-r")
        options.append(row_num.get())
    if col_num.get() != "":
        options.append("-p")
        options.append(col_num.get())
    if axis_show.get() == 0:
        options.append("--noaxis")
    if tick_show.get() == 0:
        options.append("--notick")
    if title_show.get() == 0:
        options.append("--notitle")
    if cbar_show.get() == 0:
        options.append("--nocbar")
    if title_in.get() == 1:
        options.append("--title-in")
    if title.get() != "":
        options.append("--figtitle")
        options.append(title.get())
    if fig_size_width.get() != "" and fig_size_height.get() != "":
        options.append("--figsize")
        options.append(fig_size_height.get())
        options.append(fig_size_width.get())
    if fig_w_space.get() != "":
        options.append("--wspace")
        options.append(fig_w_space.get())
    if fig_h_space.get() != "":
        options.append("--hspace")
        options.append(fig_h_space.get())
    if fig_ext.get() != "":
        options.append("--figext")
        options.append(fig_ext.get())
    if fig_num.get() != "":
        options.append("--fignum")
        options.append(fig_num.get())
    if coords.get() != "":
        options.append("--coord")
        options.append(coords.get())

    if coastline.get() != 0:
        options.append("--coastline")
    if resolution.get() != "" and resolution.get() != "None":
        options.append("--resolution")
        options.append(resolution.get())
    if lalo_label.get() != 0:
        options.append("--lalo-label")
    if lalo_step.get() != "":
        options.append("--lalo-step")
        options.append(lalo_step.get())
    if scalebar_distance.get() != "" and scalebar_lat.get() != "" and scalebar_lon.get() != "":
        options.append("--scalebar")
        options.append(scalebar_distance.get())
        options.append(scalebar_lat.get())
        options.append(scalebar_lon.get())
    if show_scalebar.get() == 0:
        options.append("--noscalebar")

    if save.get() != 0:
        options.append("--save")
    if output_file.get() != "":
        options.append("-o")

        location_parts = h5_file.get().split("/")
        location = "/".join(location_parts[1:-1])

        options.append("/" + str(location) + "/" + output_file.get())

    if show_info.get() == 1:

        if file_info_window is not None:
            file_info_window.destroy()
            file_info_window = None

        show_file_info()

    print(options)
    if h5_file.get() != "":
        if view.fig is not None:
            view.fig.clear()
        view.main(options)
    else:
        print("No file selected")


def reset_plot():
    set_variables_from_attributes()


def set_dem_short(x, y, z):
    if dem_file.get() != "":
        dem_short.set(dem_file.get().split("/")[-1])
        pick_dem_file_button.config(text="Cancel")
    else:
        dem_short.set("No File Selected")
        pick_dem_file_button.config(text="Select DEM File")


def set_mask_short(x, y, z):
    if mask_file.get() != "":
        mask_short.set(mask_file.get().split("/")[-1])
        pick_mask_file_button.config(text="Cancel")
    else:
        mask_short.set("No File Selected")
        pick_mask_file_button.config(text="Select Mask File")


def set_variables_from_attributes():

    global attributes

    set_sliders()

    dem_file.set("")

    mask_file.set("")

    unit.set(attributes['UNIT'])
    colormap.set('jet')
    projection.set("cyl")
    transparency.set(1.0)


    shading.set(1)
    countours.set(1)
    countour_smoothing.set("3.0")
    countour_step.set("200")

    try:

        subset_x_from.set(0)
        subset_y_from.set(0)
        subset_x_to.set(attributes['XMAX'])
        subset_y_to.set(attributes['YMAX'])

        '''ul_lon, ul_lat, lr_lon, lr_lat = compute_lalo(attributes['WIDTH'], attributes['LENGTH'], all_data=True)

        subset_lat_from.set(ul_lat)
        subset_lon_from.set(ul_lon)
        subset_lat_to.set(lr_lat)
        subset_lon_to.set(lr_lon)'''

        try:
            attributes['REF_X']

            ref_x.set(attributes['REF_X'])
            ref_y.set(attributes["REF_Y"])
            ref_lon_data, ref_lat_data = compute_lalo(ref_x.get(), ref_y.get())
            ref_lat.set(ref_lat_data)
            ref_lon.set(ref_lon_data)
            ref_color.set("b")
            ref_sym.set(".")
            ref_date.set(attributes['REF_DATE'])
        except:
            ref_x.set("")
            ref_y.set("")
            ref_lat.set("")
            ref_lon.set("")


    except KeyError:
        print()


    row_num.set("1")
    col_num.set("1")
    font_size.set("16")
    plot_dpi.set("200")
    fig_size_width.set("6.0")
    fig_size_height.set("8.0")
    title.set(h5_file_short.get())
    fig_w_space.set("1")
    fig_h_space.set("1")
    axis_show.set(1)
    cbar_show.set(1)
    title_show.set(1)
    tick_show.set(1)
    title_in.set(1)

    fig_ext.set(".pdf")
    fig_num.set("1")
    coords.set("geo")

    resolution.set("c")
    show_scalebar.set(1)

    '''     Set ListBox and Epoch Information   '''
    set_epoch_info()


def set_sliders(value=None):
    global starting_lower_lim, starting_upper_lim, y_lim_upper_slider, y_lim_lower_slider

    print(value)
    if value is not None:
        epoch.set(value)

    data = read_file_data(value)

    print(data)

    max_val = numpy.nanmax(data)
    min_val = numpy.nanmin(data)

    print("SLL "+str(min_val))

    starting_upper_lim = max_val * 5
    starting_lower_lim = min_val * 5

    update_sliders("m", setValues=False)

    y_lim_upper.set(max_val)
    y_lim_lower.set(min_val)


def set_epoch_info():
    global epoch, epoch_list, epoch_option_menu, excludes_list_box, attributes, ref_date_option_menu

    excludes_list_box.delete(0, END)
    epoch_option_menu['menu'].delete(0, END)

    for the_epoch in epoch_list:
        epoch_option_menu.children['menu'].add_command(label=the_epoch,
                                                       command=lambda val=the_epoch: set_sliders(val))
        if the_epoch is not "All":
            excludes_list_box.insert(END, the_epoch)

            if attributes['FILE_TYPE'] in ['timeseries']:
                ref_date_option_menu.children['menu'].add_command(label=the_epoch,
                                                                  command=lambda val=the_epoch: ref_date.set(val))
    if attributes['FILE_TYPE'] in ['timeseries']:
        ref_date.set(epoch_list[1])
    else:
        ref_date.set("")

    epoch.set("All")


def read_file_data(epoch=None):
    global epoch_list, attributes, ref_dates_list, use_default

    atr = readfile.read_attribute(h5_file.get())
    file_type = atr['FILE_TYPE']
    epoch_list = ["All"]
    ref_dates_list = [""]
    print(file_type)

    epoch_list = readfile.get_slice_list(h5_file.get())

    if file_type in readfile.multi_dataset_hdf5_file:
        use_default.set(1)

    if epoch and epoch is not "All":
        data, attributes = readfile.read(h5_file.get(), datasetName=epoch)
    else:
        data, attributes = readfile.read(h5_file.get(), datasetName=epoch_list[len(epoch_list)-1])

    return data


def compute_lalo(x, y, all_data=False):
    try:
        x_data = int(float(x))
    except:
        x_data = 0

    try:
        y_data = int(float(y))
    except:
        y_data = 0

    data_box = (0, 0, x_data, y_data)
    print("Data Box:"+str(data_box))
    lalo = ut.coordinate(attributes).box_pixel2geo(data_box)

    formatted_lalo = [str(round(num, 2)) for num in lalo]

    if all_data:
        return formatted_lalo
    else:
        return formatted_lalo[2], formatted_lalo[3]


def compute_xy(lat, lon):
    lat_data = round(float(lat), 4)
    lon_data = round(float(lon), 4)

    data_box = (float(attributes['LON_REF1']), float(attributes['LAT_REF1']), lon_data, lat_data)

    xy = ut.coordinate(attributes).box_geo2pixel(data_box)

    return str(xy[2]), str(xy[3])


def update_subset_lalo(x, y, z):
    global update_in_progress

    if update_in_progress:
        return

    update_in_progress = True
    x_from, x_to, y_from, y_to = subset_x_from.get(), subset_x_to.get(), subset_y_from.get(), subset_y_to.get()

    lon_from, lat_from = compute_lalo(x_from, y_from)
    lon_to, lat_to = compute_lalo(x_to, y_to)

    subset_lat_from.set(lat_from)
    subset_lat_to.set(lat_to)
    subset_lon_from.set(lon_from)
    subset_lon_to.set(lon_to)

    update_in_progress = False


def update_subset_xy(x, y, z):
    global update_in_progress

    if update_in_progress:
        return

    update_in_progress = True

    lat_from, lat_to, lon_from, lon_to = subset_lat_from.get(), subset_lat_to.get(), subset_lon_from.get(), subset_lon_to.get()

    x_from, y_from = compute_xy(lat_from, lon_from)
    x_to, y_to = compute_xy(lat_to, lon_to)

    subset_x_from.set(x_from)
    subset_x_to.set(x_to)
    subset_y_from.set(y_from)
    subset_y_to.set(y_to)

    update_in_progress = False


def update_reference_lalo(x, y, z):
    global update_in_progress

    if update_in_progress:
        return

    update_in_progress = True
    x, y = ref_x.get(), ref_y.get()

    lon, lat = compute_lalo(x, y)

    ref_lat.set(lat)
    ref_lon.set(lon)

    update_in_progress = False


def update_reference_xy(x, y, z):
    global update_in_progress

    if update_in_progress:
        return

    update_in_progress = True

    lat, lon = ref_lat.get(), ref_lon.get()

    x, y = compute_xy(lat, lon)

    ref_x.set(x)
    ref_y.set(y)

    update_in_progress = False


def validate_numbers(action, index, value_if_allowed, prior_value, text, validation_type, trigger_type, widget_name):

    if value_if_allowed == "" or value_if_allowed == "-":
        return True

    if text in '0123456789.-+ ':
        try:
            float(value_if_allowed)
            return True
        except ValueError:
            return False
    else:
        return False


def main():
    global canvas, frame, attributes, update_in_progress, h5_file, h5_file_short, pick_h5_file_button, mask_file, mask_short, \
        pick_mask_file_button, starting_upper_lim, y_lim_upper, y_lim_upper_slider, y_lim_lower, y_lim_lower_slider, unit, \
        colormap, projection, lr_flip, ud_flip, wrap, opposite, transparency, show_info, dem_file, dem_short, \
        pick_dem_file_button, shading, countours, countour_smoothing, countour_step, subset_x_from, subset_x_to, subset_y_from, \
        subset_y_to, subset_lat_from, subset_lat_to, subset_lon_from, subset_lon_to, ref_x, ref_y, ref_lat, ref_lon, font_size, \
        plot_dpi, row_num, col_num, axis_show, cbar_show, title_show, tick_show, title_in, title, fig_size_width, \
        fig_size_height, fig_ext, fig_num, fig_w_space, fig_h_space, coords, coastline, resolution, lalo_label, lalo_step, \
        scalebar_distance, scalebar_lat, scalebar_lon, show_scalebar, save, output_file, ref_color, ref_sym, ref_date, \
        epoch_option_menu, epoch, epoch_list, excludes_list_box, use_default, ref_date_option_menu, ref_dates_list, show_ref,\
        settings_file_button, settings_file, settings_to_variable

    '''     Setup window, widget canvas, and scrollbar. Add Submit Button to top of window      '''
    root = Tk()
    root.minsize(width=365, height=750)
    root.maxsize(width=365, height=750)
    root.resizable(width=False, height=False)

    vcmd_num = (root.register(validate_numbers), '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')

    reset_settings_file_frame = Frame(root)
    reset_settings_file_frame.pack(side=TOP)

    reset_button = Button(reset_settings_file_frame, text="Reset Settings", padx=15, command=lambda: reset_plot())
    reset_button.pack(side=LEFT, pady=(10, 5))

    settings_file = StringVar()
    settings_file_button = Button(reset_settings_file_frame, text="Select Settings File", padx=15, command=lambda: pick_settings())
    settings_file_button.pack(side=LEFT, pady=(10, 5))

    submit_button = Button(root, text="Show Plot", command=lambda: show_plot(), background="green", padx=15)
    submit_button.pack(side=TOP, pady=(10, 20))

    canvas = Canvas(root, width=345, height=680)
    canvas.pack(side=LEFT, anchor='nw')

    scrollbar = Scrollbar(root, command=canvas.yview)
    scrollbar.pack(side=LEFT, fill='y')

    canvas.configure(yscrollcommand=scrollbar.set)
    canvas.bind('<Configure>', on_configure)

    frame = Frame(canvas)
    canvas.create_window((0, 0), window=frame, anchor='nw')

    '''     Frames, Text Variables, and Widgets for selection of the timeseries.h5 file to plot data from.      '''
    pick_h5_file_frame = Frame(frame)

    h5_file = StringVar()
    h5_file_short = StringVar()
    h5_file_short.set("No File Selected")

    pick_h5_file_button = Button(pick_h5_file_frame, text='Select .h5 File', width=15,
                                 command=lambda: pick_file())
    selected_ts_file_label = Label(pick_h5_file_frame, textvariable=h5_file_short)

    '''     Frames, Text Variables, and Widgets for selection of the maskConnComp.h5 file to add a mask to the ata.     '''
    pick_mask_file_frame = Frame(frame)

    mask_file = StringVar()
    mask_file.trace('w', callback=set_mask_short)
    mask_short = StringVar()
    mask_short.set("No File Selected")

    pick_mask_file_button = Button(pick_mask_file_frame, text='Select Mask File', width=15,
                                   command=lambda: pick_mask())
    selected_mask_file_label = Label(pick_mask_file_frame, textvariable=mask_short)

    '''     WIDGETS FOR SHOWING EPOCHS AND EXLUDE DATES     '''

    epoch_labels_frame = Frame(frame)

    epoch_option_menu_label = Label(epoch_labels_frame, text="Epoch", width=10, anchor='w')
    exclude_date_label = Label(epoch_labels_frame, text="Exclude Dates", width=15, anchor='w')

    epoch_frame = Frame(frame)

    epoch = StringVar()
    epoch.set(" ")
    epoch_option_menu = OptionMenu(epoch_frame, epoch, *epoch_list, command=set_sliders)
    epoch_option_menu.config(width=10)

    excludes_list_box = Listbox(epoch_frame, selectmode=MULTIPLE, height=5)

    '''
    |-----------------------------------------------------------------------------------------------------|
    |                                                                                                     |
    |                                WIDGETS TO CONTROL DISPLAY OPTIONS                                   |
    |                                                                                                     |
    |-----------------------------------------------------------------------------------------------------|

    '''

    '''     DISPLAY OPTIONS WIDGETS'''
    display_options_label = Label(frame, text="DISPLAY OPTIONS:", anchor=W)


    '''    WIDGETS FOR UPPER AND LOWER Y-LIM      '''
    starting_upper_lim = 5000

    y_lim_frame = Frame(frame)

    y_lim_upper = DoubleVar()
    y_lim_upper.set(20)

    y_lim_upper_frame = Frame(y_lim_frame)
    y_lim_upper_label = Label(y_lim_upper_frame, text="Maximum", width=8)
    y_lim_upper_slider = Scale(y_lim_upper_frame, from_=0, to=starting_upper_lim, orient=HORIZONTAL, length=150,
                               variable=y_lim_upper, showvalue=0, resolution=0.001)
    y_lim_upper_entry = Entry(y_lim_upper_frame, textvariable=y_lim_upper, width=6, validate='key', validatecommand=vcmd_num)

    y_lim_lower = DoubleVar()
    y_lim_lower.set(-20)

    y_lim_lower_frame = Frame(y_lim_frame)
    y_lim_lower_label = Label(y_lim_lower_frame, text="Minimum", width=8)
    y_lim_lower_slider = Scale(y_lim_lower_frame, from_=0, to=starting_upper_lim, orient=HORIZONTAL, length=150,
                               variable=y_lim_lower, showvalue=0, resolution=0.001)
    y_lim_lower_entry = Entry(y_lim_lower_frame, textvariable=y_lim_lower, width=6, validate='key', validatecommand=vcmd_num)

    use_default = BooleanVar()
    use_default_y_lims = Checkbutton(y_lim_frame, text="Use Default Y Lim", variable=use_default)


    '''     WIDGETS FOR UNIT, COLORMAP, AND PROJECTION    '''
    unit_cmap_projection_labels_frame = Frame(frame)
    unit_label = Label(unit_cmap_projection_labels_frame, text="Unit", width=7, anchor='w')
    colormap_label = Label(unit_cmap_projection_labels_frame, text="Colormap", width=13, anchor='w')
    projection_label = Label(unit_cmap_projection_labels_frame, text="Projection", width=8, anchor='w')

    unit_cmap_projection_frame = Frame(frame)

    unit = StringVar()
    unit.set(" ")

    unit_option_menu = OptionMenu(unit_cmap_projection_frame, unit, *unit_options, command=scale_sliders)
    unit_option_menu.config(width=5)

    colormap = StringVar()
    colormap.set(" ")

    colormap_option_menu = OptionMenu(unit_cmap_projection_frame, colormap, *colormaps)
    colormap_option_menu.config(width=11)

    projection = StringVar()
    projection.set(" ")

    projection_option_menu = OptionMenu(unit_cmap_projection_frame, projection, *projections)
    projection_option_menu.config(width=7)


    '''     WIDGETS FOR FLIPPING, WRAP, AND OPPOSITE    '''
    flip_frame = Frame(frame)

    lr_flip = IntVar()
    ud_flip = IntVar()
    wrap = IntVar()
    opposite = IntVar()

    lr_flip_checkbutton = Checkbutton(flip_frame, text="Flip LR", variable=lr_flip)
    ud_flip_checkbutton = Checkbutton(flip_frame, text="Flip UD", variable=ud_flip)
    wrap_checkbutton = Checkbutton(flip_frame, text="Wrap", variable=wrap)
    opposite_checkbutton = Checkbutton(flip_frame, text="Opposite", variable=opposite)


    '''     WIDGETS FOR TRANSPARENCY'''
    transparency = IntVar()

    transparency_frame = Frame(frame)
    transparency_label = Label(transparency_frame, text="Alpha", width=8)
    transparency_slider = Scale(transparency_frame, from_=0, to=1, resolution=0.1, orient=HORIZONTAL, length=150,
                                variable=transparency, showvalue=0)
    transparency_entry = Entry(transparency_frame, textvariable=transparency, width=6, validate='key', validatecommand=vcmd_num)

    '''     WIDGETS FOR SHOWING INFO SCREEN'''
    show_info = IntVar()
    show_info_checkbutton = Checkbutton(frame, text="Show File Info", variable=show_info)





    '''
        |-----------------------------------------------------------------------------------------------------|
        |                                                                                                     |
        |                         WIDGETS TO CONTROL DEM TOPOGRAPHY OPTIONS                                   |
        |                                                                                                     |
        |-----------------------------------------------------------------------------------------------------|

    '''

    '''     DEM OPTIONS WIDGETS'''
    dem_options_label = Label(frame, text="DEM OPTIONS:", anchor=W)


    '''     WIDGETS FOR DEM TOPOGRAPHY FILE    '''
    pick_dem_file_frame = Frame(frame)

    dem_file = StringVar()
    dem_file.trace('w', callback=set_dem_short)
    dem_short = StringVar()
    dem_short.set("No File Selected")

    pick_dem_file_button = Button(pick_dem_file_frame, text='Select DEM File', width=15, command=lambda: pick_dem())
    selected_dem_file_label = Label(pick_dem_file_frame, textvariable=dem_short)


    '''     WIDGETS FOR DEM SHADING AND CONTOURS     '''
    dem_options_frame = Frame(frame)

    shading = IntVar()

    countours = IntVar()

    dem_shading_checkbutton = Checkbutton(dem_options_frame, text="Show Shaded Relief", variable=shading)
    dem_countours_checkbutton = Checkbutton(dem_options_frame, text="Show Countour Lines", variable=countours)


    '''     WIDGERS FOR DEN CONTOUR SMOOTHING AND STEP'''
    dem_countour_options = Frame(frame)

    countour_smoothing = StringVar()

    countour_step = StringVar()

    dem_countour_smoothing_frame = Frame(dem_countour_options, width=15)
    dem_countour_smoothing_label = Label(dem_countour_smoothing_frame, text="Contour Smoothing: ", anchor='c', width=15)
    dem_countour_smoothing_entry = Entry(dem_countour_smoothing_frame, textvariable=countour_smoothing, width=6, validate='key', validatecommand=vcmd_num)

    dem_countour_step_frame = Frame(dem_countour_options, width=15)
    dem_countour_step_label = Label(dem_countour_step_frame, text="Countour Step: ", anchor='c', width=15)
    dem_countour_step_entry = Entry(dem_countour_step_frame, textvariable=countour_step, width=6, validate='key', validatecommand=vcmd_num)





    '''
            |-----------------------------------------------------------------------------------------------------|
            |                                                                                                     |
            |                                 WIDGETS TO CONTROL SUBSET OPTIONS                                   |
            |                                                                                                     |
            |-----------------------------------------------------------------------------------------------------|

    '''

    '''     SUBSET OPTIONS WIDGETS      '''
    subset_label = Label(frame, text="SUBSET DATA", anchor=W)


    '''     WIDGETS FOR SUBSET X-VALUES'''
    subset_x_frame = Frame(frame)

    subset_x_from = StringVar()
    subset_x_from.trace('w', callback=update_subset_lalo)
    subset_x_from_label = Label(subset_x_frame, text="X         From: ")
    subset_x_from_entry = Entry(subset_x_frame, textvariable=subset_x_from, width=6, validate='key', validatecommand=vcmd_num)

    subset_x_to = StringVar()
    subset_x_to.trace('w', callback=update_subset_lalo)
    subset_x_to_label = Label(subset_x_frame, text="To: ")
    subset_x_to_entry = Entry(subset_x_frame, textvariable=subset_x_to, width=6, validate='key', validatecommand=vcmd_num)

    '''     WIDGETS FOR SUBSET Y-VALUES     '''
    subset_y_frame = Frame(frame)

    subset_y_from = StringVar()
    subset_y_from.trace('w', callback=update_subset_lalo)
    subset_y_from_label = Label(subset_y_frame, text="Y         From: ")
    subset_y_from_entry = Entry(subset_y_frame, textvariable=subset_y_from, width=6, validate='key', validatecommand=vcmd_num)

    subset_y_to = StringVar()
    subset_y_to.trace('w', callback=update_subset_lalo)
    subset_y_to_label = Label(subset_y_frame, text="To: ")
    subset_y_to_entry = Entry(subset_y_frame, textvariable=subset_y_to, width=6, validate='key', validatecommand=vcmd_num)


    '''     WIDGETS FOR SUBSET LAT-VALUES       '''
    subset_lat_frame = Frame(frame)

    subset_lat_from = StringVar()
    subset_lat_from.trace('w', callback=update_subset_xy)
    subset_lat_from_label = Label(subset_lat_frame, text="Lat      From: ")
    subset_lat_from_entry = Entry(subset_lat_frame, textvariable=subset_lat_from, width=6, validate='key', validatecommand=vcmd_num)

    subset_lat_to = StringVar()
    subset_lat_to.trace('w', callback=update_subset_xy)
    subset_lat_to_label = Label(subset_lat_frame, text="To: ")
    subset_lat_to_entry = Entry(subset_lat_frame, textvariable=subset_lat_to, width=6, validate='key', validatecommand=vcmd_num)


    '''     WIDGETS FOR SUBSET LON-VALUES       '''
    subset_lon_frame = Frame(frame)

    subset_lon_from = StringVar()
    subset_lon_from.trace('w', callback=update_subset_xy)
    subset_lon_from_label = Label(subset_lon_frame, text="Lon      From: ")
    subset_lon_from_entry = Entry(subset_lon_frame, textvariable=subset_lon_from, width=6, validate='key', validatecommand=vcmd_num)

    subset_lon_to = StringVar()
    subset_lon_to.trace('w', callback=update_subset_xy)
    subset_lon_to_label = Label(subset_lon_frame, text="To: ")
    subset_lon_to_entry = Entry(subset_lon_frame, textvariable=subset_lon_to, width=6, validate='key', validatecommand=vcmd_num)





    '''
            |-----------------------------------------------------------------------------------------------------|
            |                                                                                                     |
            |                              WIDGETS TO CONTROL REFERENCE OPTIONS                                   |
            |                                                                                                     |
            |-----------------------------------------------------------------------------------------------------|

    '''

    '''     REFERENCE OPTIONS WIDGETS     '''
    reference_label = Label(frame, text="REFERENCE:", anchor=W)


    '''     WIDGETS FOR REFERENCE XY'''
    ref_xy_frame = Frame(frame)

    ref_x = StringVar()
    ref_x.trace('w', callback=update_reference_lalo)
    ref_x_label = Label(ref_xy_frame, text="X:    ")
    ref_x_entry = Entry(ref_xy_frame, textvariable=ref_x, width=6, validate='key', validatecommand=vcmd_num)

    ref_y = StringVar()
    ref_y.trace('w', callback=update_reference_lalo)
    ref_y_label = Label(ref_xy_frame, text="Y:    ")
    ref_y_entry = Entry(ref_xy_frame, textvariable=ref_y, width=6, validate='key', validatecommand=vcmd_num)


    '''     WIDGETS FOR REFERENCE LALO'''
    ref_latlon_frame = Frame(frame)

    ref_lat = StringVar()
    ref_lat.trace('w', callback=update_reference_xy)
    ref_lat_label = Label(ref_latlon_frame, text="Lat: ")
    ref_lat_entry = Entry(ref_latlon_frame, textvariable=ref_lat, width=6, validate='key', validatecommand=vcmd_num)

    ref_lon = StringVar()
    ref_lon.trace('w', callback=update_reference_xy)
    ref_lon_label = Label(ref_latlon_frame, text="Lon: ")
    ref_lon_entry = Entry(ref_latlon_frame, textvariable=ref_lon, width=6, validate='key', validatecommand=vcmd_num)


    '''     WIDGETS FOR SHOWING REFERENCE MARKER        '''
    show_ref_frame = Frame(frame)

    show_ref = IntVar()
    show_ref_checkbutton = Checkbutton(show_ref_frame, text="Show Reference", variable=show_ref)


    '''     WIDGETS FOR REFERENCE OPTIONS       '''
    reference_options_labels_frame = Frame(frame)

    ref_color_label = Label(reference_options_labels_frame, text="Ref Color", width=10, anchor='w')
    ref_symbol_label = Label(reference_options_labels_frame, text="Ref Symbol", width=10, anchor='w')
    ref_date_label = Label(reference_options_labels_frame, text="Ref Date", width=10, anchor='w')

    reference_options_frame = Frame(frame)

    reference_colors = ["b", "g", "r", "m", "c", "y", "k", "w"]
    ref_color = StringVar()
    ref_color.set(" ")
    ref_color_option_menu = OptionMenu(reference_options_frame, ref_color, *reference_colors)
    ref_color_option_menu.config(width=8)

    reference_symbols = [".", ",", "o", "v", "^", "<", ">", "1", "2", "3", "4", "8", "s", "p", "P", "*", "h", "H", "+",
                         "x", "X", "d", "D", "|", "_", " ", " "]
    ref_sym = StringVar()
    ref_sym.set(" ")
    ref_symbol_option_menu = OptionMenu(reference_options_frame, ref_sym,*reference_symbols)
    ref_symbol_option_menu.config(width=8)

    ref_date = StringVar()
    ref_date.set(" ")
    ref_date_option_menu = OptionMenu(reference_options_frame, ref_date, *ref_dates_list)
    ref_date_option_menu.config(width=8)





    '''
            |-----------------------------------------------------------------------------------------------------|
            |                                                                                                     |
            |                                 WIDGETS TO CONTROL FIGURE OPTIONS                                   |
            |                                                                                                     |
            |-----------------------------------------------------------------------------------------------------|

    '''
    '''     FIGURE OPTIONS WIDGETS'''
    figure_label = Label(frame, text="FIGURE:", anchor=W)


    '''     WIDGETS FOR FONT SIZE AND FIGURE DPI'''
    font_dpi_frame = Frame(frame)

    font_size = StringVar()
    font_size_label = Label(font_dpi_frame, text="Font Size:    ")
    font_size_entry = Entry(font_dpi_frame, textvariable=font_size, width=6, validate='key', validatecommand=vcmd_num)

    plot_dpi = StringVar()
    dpi_label = Label(font_dpi_frame, text="DPI:    ")
    dpi_entry = Entry(font_dpi_frame, textvariable=plot_dpi, width=6, validate='key', validatecommand=vcmd_num)


    '''     WIDGETS FOR NUMBER OF ROWS AND NUMBER OF COLUMNS      '''
    row_col_num_frame = Frame(frame)

    row_num = StringVar()
    row_num_label = Label(row_col_num_frame, text="Row Num:   ")
    row_num_entry = Entry(row_col_num_frame, textvariable=row_num, width=6, validate='key', validatecommand=vcmd_num)

    col_num = StringVar()
    col_num_label = Label(row_col_num_frame, text="Col Num:   ")
    col_num_entry = Entry(row_col_num_frame, textvariable=col_num, width=6, validate='key', validatecommand=vcmd_num)


    '''     WIDGETS FOR SHOWING AXIS, COLORBAR, TITLE, AND AXIS TICKS       '''
    axis_cbar_frame = Frame(frame)

    axis_show = IntVar()
    axis_show_checkbutton = Checkbutton(axis_cbar_frame, text="Show Axis", variable=axis_show)

    cbar_show = IntVar()
    cbar_show_checkbutton = Checkbutton(axis_cbar_frame, text="Show Colorbar", variable=cbar_show)

    title_show = IntVar()
    title_show_checkbutton = Checkbutton(axis_cbar_frame, text="Show Title", variable=title_show)

    title_tick_frame = Frame(frame)

    tick_show = IntVar()
    tick_show_checkbutton = Checkbutton(title_tick_frame, text="Show Ticks", variable=tick_show)

    title_in = IntVar()
    title_in_checkbutton = Checkbutton(title_tick_frame, text="Title in Axes", variable=title_in)


    '''     WIDGETS FOR TITLE     '''
    title_input_frame = Frame(frame)

    title = StringVar()
    title_input_label = Label(title_input_frame, text="Figure Title: ")
    title_input_entry = Entry(title_input_frame, textvariable=title, width=25)


    '''     WIDGETS FO DIGURE SIZE      '''
    fig_size_frame = Frame(frame)

    fig_size_label = Label(fig_size_frame, text="Fig Size")

    fig_size_width = StringVar()
    fig_size_width_label = Label(fig_size_frame, text="Width: ")
    fig_size_width_entry = Entry(fig_size_frame, textvariable=fig_size_width, width=6, validate='key', validatecommand=vcmd_num)

    fig_size_height = StringVar()
    fig_size_height_label = Label(fig_size_frame, text="Length: ")
    fig_size_height_entry = Entry(fig_size_frame, textvariable=fig_size_height, width=6, validate='key', validatecommand=vcmd_num)


    '''     WIDGETS FOR FIGURE EXTENIONS AND FIGURE NUMBERS'''
    fig_ext_num_label_frame = Frame(frame)
    fig_ext_label = Label(fig_ext_num_label_frame, text="Fig Ext", width=14, anchor='w')
    fig_num_label = Label(fig_ext_num_label_frame, text="Fig Num", width=14, anchor='w')

    fig_ext_num_frame = Frame(frame)

    figure_extensions = [".emf", ".eps", ".pdf", ".png", ".ps", ".raw", ".rgba", ".svg", ".svgz"]
    fig_ext = StringVar()
    fig_ext.set(" ")
    fig_ext_option_menu = OptionMenu(fig_ext_num_frame, fig_ext, *figure_extensions)
    fig_ext_option_menu.config(width=14)

    figure_numbers = ["1", "2", "3", "4", "5"]
    fig_num = StringVar()
    fig_num.set(" ")
    fig_num_option_menu = OptionMenu(fig_ext_num_frame, fig_num, *figure_numbers)
    fig_num_option_menu.config(width=14)


    '''     WIDGETS FOR FIGURE WIDTH AND HIEGHT SPACE'''
    fig_w_space_frame = Frame(frame)

    fig_w_space = StringVar()
    fig_w_space_label = Label(fig_w_space_frame, text="Fig Width Space: ")
    fig_w_space_entry = Entry(fig_w_space_frame, textvariable=fig_w_space, width=6, validate='key', validatecommand=vcmd_num)

    fig_h_space_frame = Frame(frame)

    fig_h_space = StringVar()
    fig_h_space_label = Label(fig_h_space_frame, text="Fig Height Space:")
    fig_h_space_entry = Entry(fig_h_space_frame, textvariable=fig_h_space, width=6, validate='key', validatecommand=vcmd_num)


    '''     WIDGETS FOR COORDINATE TYPE     '''
    coords_frame = Frame(frame)

    coords = StringVar()
    coords.set(" ")
    coords_label = Label(coords_frame, text="Coords: ", width=5)
    coords_option_menu = OptionMenu(*(coords_frame, coords) + tuple(["radar", "geo"]))
    coords_option_menu.config(width=15)





    '''
            |-----------------------------------------------------------------------------------------------------|
            |                                                                                                     |
            |                                    WIDGETS TO CONTROL MAP OPTIONS                                   |
            |                                                                                                     |
            |-----------------------------------------------------------------------------------------------------|

    '''

    '''     MAP OPTIONS WIDGETS     '''
    map_options_label = Label(frame, text="MAP: ", anchor=W)


    '''     WIDGETS FOR COASTLINE'''
    coastline_res_frame = Frame(frame)

    coastline = IntVar()
    coastline_checkbutton = Checkbutton(coastline_res_frame, text="Show Coastline", variable=coastline)


    '''     WODGETS FOR RESOLUTION      '''
    resolution_label = Label(coastline_res_frame, text="Res: ", width=3, anchor='w')
    resolution = StringVar()
    resolution.set(" ")
    resolution_option_menu = OptionMenu(*(coastline_res_frame, resolution) + tuple(["c", "l", "i", "h", "f", "None"]))
    resolution_option_menu.config(width=8)

    '''     WIDGETS FOR LALO LABEL AND LALO STEP'''
    lalo_settings_frame = Frame(frame)

    lalo_label = IntVar()
    lalo_label_checkbutton = Checkbutton(lalo_settings_frame, text="Show LALO Label", variable=lalo_label)

    lalo_step = StringVar()
    lalo_step_label = Label(lalo_settings_frame, text="LALO Step: ")
    lalo_step_entry = Entry(lalo_settings_frame, textvariable=lalo_step, width=6, validate='key', validatecommand=vcmd_num)


    '''     WIDGETS FOR SCALEBAR DISTANCE, LATITUDE, LONGITUDE      '''
    scalebar_settings = Frame(frame)

    scalebar_distance = StringVar()
    scalebar_lat = StringVar()
    scalebar_lon = StringVar()

    scalebar_label = Label(scalebar_settings, text="Scalebar")
    scalebar_dist_label = Label(scalebar_settings, text="Dist: ")
    scalebar_dist_entry = Entry(scalebar_settings, textvariable=scalebar_distance, width=4, validate='key', validatecommand=vcmd_num)
    scalebar_lat_label = Label(scalebar_settings, text="Lat/Lon: ")
    scalebar_lat_entry = Entry(scalebar_settings, textvariable=scalebar_lat, width=4, validate='key', validatecommand=vcmd_num)
    scalebar_lon_entry = Entry(scalebar_settings, textvariable=scalebar_lon, width=4, validate='key', validatecommand=vcmd_num)


    '''     WIDGETA FOR SHOWING SCALEBAR     '''
    show_scalebar_frame = Frame(frame)

    show_scalebar = IntVar()
    show_scalebar_checkbutton = Checkbutton(show_scalebar_frame, text="Show Scalebar", variable=show_scalebar)





    '''
                |-----------------------------------------------------------------------------------------------------|
                |                                                                                                     |
                |                                 WIDGETS TO CONTROL OUTPUT OPTIONS                                   |
                |                                                                                                     |
                |-----------------------------------------------------------------------------------------------------|

    '''

    '''     OUTPUT OPTIONS WIDGETS      '''
    output_label = Label(frame, text="OUTPUT", anchor=W)

    output_frame = Frame(frame)

    '''     WIDGETS FOR SAVE    '''
    save = IntVar()
    save_checkbutton = Checkbutton(output_frame, text="Save Output", variable=save)

    '''     WIDGETS FOR OUTPU FILE      '''
    output_file = StringVar()
    output_file_label = Label(output_frame, text="Output File: ")
    output_file_entry = Entry(output_frame, textvariable=output_file, width=12)

    save_settings_button = Button(frame, text="Save Settings", command=write_settings_file, padx=15)










    '''
                |-----------------------------------------------------------------------------------------------------|
                |                                                                                                     |
                |                         PACKING AND PLACEMENT COMMANDS FOR ALL WIDGETS                              |
                |                                                                                                     |
                |-----------------------------------------------------------------------------------------------------|

    '''

    pick_h5_file_frame.pack(anchor='w', fill=X, pady=(10, 5), padx=10)
    pick_h5_file_button.pack(side=LEFT, anchor='w', padx=(0, 20))
    selected_ts_file_label.pack(side=LEFT, fill=X)

    pick_mask_file_frame.pack(anchor='w', fill=X)
    pick_mask_file_button.pack(side=LEFT, anchor='w', pady=5, padx=(10, 20))
    selected_mask_file_label.pack(side=LEFT, fill=X)

    epoch_labels_frame.pack(anchor='w', fill=X, pady=5)
    epoch_option_menu_label.pack(side=LEFT, padx=(12, 20))
    exclude_date_label.pack(side=LEFT)

    epoch_frame.pack(anchor='w', fill=X)
    epoch_option_menu.pack(side=LEFT, anchor='n', pady=5, padx=(10, 20))
    excludes_list_box.pack(side=LEFT, pady=5)


    display_options_label.pack(anchor='w', fill=X, pady=(35, 0), padx=10)

    y_lim_frame.pack(fill=X, pady=10, padx=10)

    y_lim_upper_frame.pack(side=TOP, fill=X, pady=(0, 10))
    y_lim_upper_label.pack(side=LEFT)
    y_lim_upper_slider.pack(side=LEFT, padx=10)
    y_lim_upper_entry.pack(side=LEFT)

    y_lim_lower_frame.pack(side=TOP, fill=X)
    y_lim_lower_label.pack(side=LEFT)
    y_lim_lower_slider.pack(side=LEFT, padx=10)
    y_lim_lower_entry.pack(side=LEFT)

    use_default_y_lims.pack(side=LEFT, pady=(10, 5))

    unit_cmap_projection_labels_frame.pack(anchor='w', fill=X, pady=(10, 0), padx=10)
    unit_label.pack(side=LEFT, padx=(0, 20))
    colormap_label.pack(side=LEFT, padx=(0, 20))
    projection_label.pack(side=LEFT)

    unit_cmap_projection_frame.pack(anchor='w', fill=X, pady=(10, 5), padx=10)
    unit_option_menu.pack(side=LEFT, padx=(0, 10), fill=X)
    colormap_option_menu.pack(side=LEFT, padx=(0, 10))
    projection_option_menu.pack(side=LEFT)

    flip_frame.pack(anchor='w', fill=X, pady=10, padx=10)
    lr_flip_checkbutton.pack(side=LEFT, padx=(0, 12))
    ud_flip_checkbutton.pack(side=LEFT, padx=(0, 12))
    wrap_checkbutton.pack(side=LEFT, padx=(0, 12))
    opposite_checkbutton.pack(side=LEFT)

    transparency_frame.pack(side=TOP, fill=X)
    transparency_label.pack(side=LEFT)
    transparency_slider.pack(side=LEFT, padx=10)
    transparency_entry.pack(side=LEFT)

    show_info_checkbutton.pack(anchor='center', pady=10)

    dem_options_label.pack(anchor='w', fill=X, pady=(35, 0), padx=10)

    pick_dem_file_frame.pack(fill=X, pady=10)
    pick_dem_file_button.pack(side=LEFT, anchor='w', pady=5, padx=(10, 20))
    selected_dem_file_label.pack(side=LEFT, fill=X)

    dem_options_frame.pack(anchor='w', fill=X, pady=(0, 10), padx=10)
    dem_shading_checkbutton.pack(side=LEFT, padx=(0, 12))
    dem_countours_checkbutton.pack(side=LEFT)

    dem_countour_options.pack(anchor='w', fill=X, pady=10)

    dem_countour_smoothing_frame.pack(anchor='w', side=LEFT, pady=(0, 10), padx=(20, 10))
    dem_countour_smoothing_label.pack(side=TOP, pady=(0, 10), fill=X)
    dem_countour_smoothing_entry.pack(side=TOP)

    dem_countour_step_frame.pack(anchor='w', side=LEFT, pady=(5, 10), padx=10)
    dem_countour_step_label.pack(side=TOP, pady=(0, 10), fill=X)
    dem_countour_step_entry.pack(side=TOP)

    subset_label.pack(anchor='w', fill=X, pady=(15, 0), padx=10)

    subset_x_frame.pack(anchor='w', fill=X, pady=10, padx=10)

    subset_x_from_label.pack(side=LEFT, padx=(0, 5))
    subset_x_from_entry.pack(side=LEFT, padx=(0, 10))
    subset_x_to_label.pack(side=LEFT, padx=(0, 5))
    subset_x_to_entry.pack(side=LEFT, padx=(0, 10))

    subset_y_frame.pack(anchor='w', fill=X, pady=10, padx=10)

    subset_y_from_label.pack(side=LEFT, padx=(0, 5))
    subset_y_from_entry.pack(side=LEFT, padx=(0, 10))
    subset_y_to_label.pack(side=LEFT, padx=(0, 5))
    subset_y_to_entry.pack(side=LEFT, padx=(0, 10))

    subset_lat_frame.pack(anchor='w', fill=X, pady=10, padx=10)

    subset_lat_from_label.pack(side=LEFT, padx=(0, 5))
    subset_lat_from_entry.pack(side=LEFT, padx=(0, 10))
    subset_lat_to_label.pack(side=LEFT, padx=(0, 5))
    subset_lat_to_entry.pack(side=LEFT, padx=(0, 10))

    subset_lon_frame.pack(anchor='w', fill=X, pady=10, padx=10)

    subset_lon_from_label.pack(side=LEFT, padx=(0, 5))
    subset_lon_from_entry.pack(side=LEFT, padx=(0, 10))
    subset_lon_to_label.pack(side=LEFT, padx=(0, 5))
    subset_lon_to_entry.pack(side=LEFT, padx=(0, 10))

    reference_label.pack(fill=X, padx=10, pady=(35, 10))

    ref_xy_frame.pack(anchor='w', fill=X, pady=10, padx=10)

    ref_x_label.pack(side=LEFT, padx=(0, 5))
    ref_x_entry.pack(side=LEFT, padx=(0, 10))
    ref_y_label.pack(side=LEFT, padx=(0, 5))
    ref_y_entry.pack(side=LEFT, padx=(0, 10))

    ref_latlon_frame.pack(anchor='w', fill=X, pady=10, padx=10)

    ref_lat_label.pack(side=LEFT, padx=(0, 5))
    ref_lat_entry.pack(side=LEFT, padx=(0, 10))
    ref_lon_label.pack(side=LEFT, padx=(0, 5))
    ref_lon_entry.pack(side=LEFT, padx=(0, 10))

    show_ref_frame.pack(anchor='w', fill=X, padx=10, pady=(0, 10))
    show_ref_checkbutton.pack(side=LEFT, pady=10)

    reference_options_labels_frame.pack(anchor='w', fill=X, pady=(0, 0), padx=10)
    ref_color_label.pack(side=LEFT, padx=(0, 20))
    ref_symbol_label.pack(side=LEFT, padx=(0, 20))
    ref_date_label.pack(side=LEFT)

    reference_options_frame.pack(anchor='w', fill=X, pady=(5, 10), padx=10)
    ref_color_option_menu.pack(side=LEFT, padx=(0, 10))
    ref_symbol_option_menu.pack(side=LEFT, padx=(0, 10))
    ref_date_option_menu.pack(side=LEFT)

    figure_label.pack(fill=X, padx=10, pady=(35, 10))

    font_dpi_frame.pack(anchor='w', fill=X, padx=10, pady=10)
    font_size_label.pack(side=LEFT, padx=(0, 5))
    font_size_entry.pack(side=LEFT, padx=(0, 10))
    dpi_label.pack(side=LEFT, padx=(0, 5))
    dpi_entry.pack(side=LEFT)

    row_col_num_frame.pack(anchor='w', fill=X, padx=10, pady=10)
    row_num_label.pack(side=LEFT, padx=(0, 5))
    row_num_entry.pack(side=LEFT, padx=(0, 10))
    col_num_label.pack(side=LEFT, padx=(0, 5))
    col_num_entry.pack(side=LEFT)

    axis_cbar_frame.pack(anchor='w', fill=X, pady=10, padx=10)
    axis_show_checkbutton.pack(side=LEFT, padx=(0, 12))
    cbar_show_checkbutton.pack(side=LEFT, padx=(0, 12))
    title_show_checkbutton.pack(side=LEFT, padx=(0, 12))

    title_tick_frame.pack(anchor='w', fill=X, pady=10, padx=10)
    tick_show_checkbutton.pack(side=LEFT, padx=(0, 12))
    title_in_checkbutton.pack(side=LEFT)

    title_input_frame.pack(anchor='w', fill=X, pady=10, padx=10)
    title_input_label.pack(side=LEFT, padx=(0, 10))
    title_input_entry.pack(side=LEFT)

    fig_size_frame.pack(anchor='w', fill=X, pady=10, padx=10)
    fig_size_label.pack(side=LEFT, padx=(0, 10))
    fig_size_width_label.pack(side=LEFT, padx=(0, 5))
    fig_size_width_entry.pack(side=LEFT, padx=(0, 10))
    fig_size_height_label.pack(side=LEFT, padx=(0, 5))
    fig_size_height_entry.pack(side=LEFT, padx=(0, 10))

    fig_ext_num_label_frame.pack(anchor='w', fill=X, padx=10, pady=(10, 0))
    fig_ext_label.pack(side=LEFT, padx=(0, 20))
    fig_num_label.pack(side=LEFT)

    fig_ext_num_frame.pack(anchor='w', fill=X, padx=10, pady=(5, 10))
    fig_ext_option_menu.pack(side=LEFT, padx=(0, 10))
    fig_num_option_menu.pack(side=LEFT, padx=(0, 10))

    fig_w_space_frame.pack(anchor='w', fill=X, padx=10, pady=10)
    fig_w_space_label.pack(side=LEFT, padx=(0, 5))
    fig_w_space_entry.pack(side=LEFT, padx=(0, 10))

    fig_h_space_frame.pack(anchor='w', fill=X, padx=10, pady=10)
    fig_h_space_label.pack(side=LEFT, padx=(0, 5))
    fig_h_space_entry.pack(side=LEFT)

    coords_frame.pack(anchor='w', fill=X, padx=10, pady=10)
    coords_label.pack(side=LEFT, padx=(0, 10))
    coords_option_menu.pack(side=LEFT)

    map_options_label.pack(fill=X, padx=10, pady=(35, 10))

    coastline_res_frame.pack(anchor='w', fill=X, padx=10, pady=10)
    coastline_checkbutton.pack(side=LEFT, padx=(0, 25))
    resolution_label.pack(side=LEFT, padx=(0, 10))
    resolution_option_menu.pack(side=LEFT)

    lalo_settings_frame.pack(fill=X, padx=10, pady=10)
    lalo_label_checkbutton.pack(side=LEFT, padx=(0, 15))
    lalo_step_label.pack(side=LEFT, padx=(0, 5))
    lalo_step_entry.pack(side=LEFT)

    scalebar_settings.pack(fill=X, padx=10, pady=10)
    scalebar_label.pack(side=LEFT, padx=(0, 10))
    scalebar_dist_label.pack(side=LEFT, padx=(0, 2))
    scalebar_dist_entry.pack(side=LEFT, padx=(0, 5))
    scalebar_lat_label.pack(side=LEFT, padx=(0, 2))
    scalebar_lat_entry.pack(side=LEFT)
    scalebar_lon_entry.pack(side=LEFT)

    show_scalebar_frame.pack(fill=X, padx=10, pady=5)
    show_scalebar_checkbutton.pack(side=LEFT, pady=10)

    output_label.pack(anchor='w', fill=X, padx=10, pady=(35, 10))

    output_frame.pack(anchor='w', fill=X, padx=10, pady=10)
    save_checkbutton.pack(side=LEFT, padx=(0, 10))
    output_file_label.pack(side=LEFT, padx=(0, 5))
    output_file_entry.pack(side=LEFT)

    save_settings_button.pack(side=TOP)

    space = Frame(frame)
    space.config(height=50)
    space.pack(side=LEFT)

    if inps.file:
        on_file_selection(inps.file)

    if inps.settings:
        parse_settings()
       #sys.exit(0)

    settings_to_variable = {

        'input.file': h5_file,
        'input.epoch': epoch,
        'input.mask': mask_file,
        'output.save': save,
        'output.outfile': output_file,

        'display.min': y_lim_lower,
        'display.max': y_lim_upper,
        'display.unit': unit,
        'display.colormap': colormap,
        'display.projection': projection,
        'display.wrap': wrap,
        'display.opposite': opposite,
        'display.flip_lr': lr_flip,
        'display.flip_ud': ud_flip,
        'display.alpha': transparency,

        'dem.dem_file': dem_file,
        'dem.no_shade': shading,
        'dem.no_contour': countours,
        'dem.contour_smooth': countour_smoothing,
        'dem.contour_step': countour_step,

        'subset.x_start': subset_x_from,
        'subset.x_end': subset_x_to,
        'subset.y_start': subset_y_from,
        'subset.y_end': subset_y_to,
        'subset.lat_start': subset_lat_from,
        'subset.lat_end': subset_lat_to,
        'subset_lon_start': subset_lon_from,
        'subset_lon_end': subset_lon_to,

        'ref.ref_date': ref_date,
        'ref.ref_lat': ref_lat,
        'ref.ref_lon': ref_lon,
        'ref.ref_x': ref_x,
        'ref.ref_y': ref_y,
        'ref.no_reference': show_ref,
        'ref.ref_color': ref_color,
        'ref.ref_symbol': ref_sym,

        'fig.font_size': font_size,
        'fig.dpi': plot_dpi,
        'fig.row': row_num,
        'fig.col': col_num,
        'fig.axis': axis_show,
        'fig.ticks': tick_show,
        'fig.title_in': title_in,
        'fig.fig_title': title,
        'fig.fig_size_w': fig_size_width,
        'fig.fig_size_h': fig_size_height,
        'fig.fig_ext': fig_ext,
        'fig.fig_num': fig_num,
        'fig.width_space': fig_w_space,
        'fig.height_space': fig_h_space,
        'fig.coord': coords,

        'map.coastline': coastline,
        'map.resolution': resolution,
        'map.lalo_label': lalo_label,
        'map.lalo_step': lalo_step,
        'map.scalebar_dist': scalebar_distance,
        'map.scalebar_lat': scalebar_lat,
        'map.scalebar_lon': scalebar_lon,
        'map.no_scalebar': show_scalebar

    }

    mainloop()


if __name__ == '__main__':
    inps = cmdLineParse(sys.argv[1:])
    main()
