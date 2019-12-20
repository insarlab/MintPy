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

from mintpy.objects import timeseries
from mintpy.utils import readfile, utils as ut
from mintpy import tsview as ts_view, view, info


canvas, frame, h5_file, h5_file_short, pick_h5_file_button, mask_file, mask_short, \
pick_mask_file_button, starting_upper_lim, y_lim_upper, y_lim_upper_slider, y_lim_lower, y_lim_lower_slider, unit, \
colormap, projection, lr_flip, ud_flip, wrap, opposite, transparency, show_info, dem_file, dem_short, \
pick_dem_file_button, shading, countours, countour_smoothing, countour_step, pix_input_xy_x, pix_input_xy_y, pix_input_lalo_la, \
pix_input_lalo_lo, ref_pix_input_xy_x, ref_pix_input_xy_y, ref_pix_input_lalo_la, ref_pix_input_lalo_lo, ref_x, ref_y, \
ref_lat, ref_lon, ref_color, ref_sym, font_size, title_show, marker_size, edge_width, no_flip, zfirst, title_show, tick_show, \
title_in, title, fig_size_width, fig_size_height, fig_ext, fig_num, fig_w_space, fig_h_space, coords, coastline, resolution, \
lalo_label, lalo_step, scalebar_distance, scalebar_lat, scalebar_lon, show_scalebar, save, output_file, ref_date_option_menu, \
ref_date, excludes_list_box, num, num_option_menu\
    = None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, \
      None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, \
      None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, \
      None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, \
      None

ref_dates_list = ["All"]
num_list = [""]
settings_file = None
settings_file_button = None

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

unit_options = ["cm", "m", "dm", "km", "", "cm/yr", "m/yr", "dm/yr", "km/yr"]

attributes = []
update_in_progress = False
starting_lower_lim = 0
current_slider_scale = 1.0
inps = None
file_base = "/"
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

                    print(key)
                    print(settings_to_variable[key])
                    settings_to_variable[key].set(value)


def write_settings_file():

    output = ""

    for setting in sorted(settings_to_variable.keys()):

        value = str(settings_to_variable[setting].get())

        if setting in [title_show, no_flip, zfirst] and value is 0:
            value = "False"
        elif setting in [title_show, no_flip, zfirst] and value is 1:
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
    global h5_file, h5_file_short, attributes, starting_upper_lim, ref_date_option_menu, ref_dates_list, ref_date, y_lim_upper_slider, \
        y_lim_lower_slider, num_option_menu, num_list, file_base

    if h5_file.get() == "":
        filename = filedialog.askopenfilename(initialdir=file_base, title="Select file",
                                              filetypes=(("jpeg files", "*.h5"), ("all files", "*.*")))
        frame.filename = filename

        on_file_selection(filename)

        #return frame.filename

    else:
        h5_file.set("")
        h5_file_short.set("No File Selected")
        pick_h5_file_button.config(text="Select .h5 File")
        ref_date_option_menu['menu'].delete(1, 'end')
        num_option_menu['menu'].delete(0, 'end')
        excludes_list_box.delete(0, 'end')



def pick_mask():
    if mask_file.get() == "":
        filename = filedialog.askopenfilename(initialdir=file_base, title="Select file",
                                              filetypes=(("jpeg files", "*.h5"), ("all files", "*.*")))
        frame.filename = filename
        mask_file.set(frame.filename)
        mask_short.set(filename.split("/")[-1])
        pick_mask_file_button.config(text="Cancel")
        return frame.filename
    else:
        mask_file.set("")
        mask_short.set("No File Selected")
        pick_mask_file_button.config(text="Select Mask File")

    set_mask_short()


def pick_dem():
    if dem_file.get() == "":
        filename = filedialog.askopenfilename(initialdir=file_base, title="Select file",
                                              filetypes=(("jpeg files", "*.h5"), ("all files", "*.*")))
        frame.filename = filename
        dem_file.set(frame.filename)
        dem_short.set(filename.split("/")[-1])
        pick_dem_file_button.config(text="Cancel")
        return frame.filename
    else:
        dem_file.set("")
        dem_short.set("No File Selected")
        pick_dem_file_button.config(text="Select Topography File")

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


def set_dem_short(x, y, z):
    if dem_file.get() != "":
        dem_short.set(dem_file.get().split("/")[-1])
        pick_dem_file_button.config(text="Cancel")
    else:
        dem_short.set("No File Selected")
        pick_dem_file_button.config(text="Select Topography File")


def set_mask_short(x, y, z):
    if mask_file.get() != "":
        mask_short.set(mask_file.get().split("/")[-1])
        pick_mask_file_button.config(text="Cancel")
    else:
        mask_short.set("No File Selected")
        pick_mask_file_button.config(text="Select Mask File")

def on_configure(event):
    canvas.configure(scrollregion=canvas.bbox('all'))


def on_file_selection(file):
    global file_base
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

    set_variables_from_attributes()


def update_sliders(unit, setValues=True):
    global current_slider_scale

    scale = 1.0
    new_max = starting_upper_lim
    new_min = starting_lower_lim

    if unit in ["m", "m/yr"]:
        scale = 1.0
    elif unit in ["cm", "cm/yr"]:
        scale = 100.0
    elif unit in ["mm", "mm/yr"]:
        scale = 1000.0
    elif unit in ["dm", "dm/yr"]:
        scale = 10.0
    elif unit in ["km", "km/yr"]:
        scale = 0.001

    current_slider_scale = scale

    y_lim_upper_slider.configure(to_= new_max*scale)
    y_lim_upper_slider.configure(from_= -1*(float(new_max*scale)/4))
    y_lim_lower_slider.configure(to_= new_max*scale)
    y_lim_lower_slider.configure(from_= (float(new_min*scale)))

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
        current_slider_scale = 10.0
    elif unit in ["km", "km/yr"]:
        current_slider_scale = 0.001

    new_max_slider_val *= current_slider_scale
    new_min_slider_val *= current_slider_scale

    update_sliders(unit, setValues=False)

    y_lim_upper_slider.set(new_max_slider_val)
    y_lim_lower_slider.set(new_min_slider_val)



def show_file_info(file_info):
    window = Tk()
    window.minsize(width=350, height=550)
    window.maxsize(height=550)
    window.resizable(width=True, height=False)

    text_box = Text(window, wrap=NONE)
    text_box.insert(END, file_info)
    text_box.config(height=550)
    text_box.config(state=DISABLED)

    text_box.pack(fill=X)


def show_plot():

    options = [inps.file]

    if ref_date.get() != "All":
        options.append("--ref-date")
        options.append(ref_date.get())

    options += ["--ylim-mat", str(y_lim_lower.get()), str(y_lim_upper.get())]

    if mask_file.get() != "":
        options.append("--mask")
        options.append(mask_file.get())

    excludes = [str(excludes_list_box.get(idx)) for idx in excludes_list_box.curselection()]
    if len(excludes) > 0:
        options.append("--exclude")
        for ex in excludes:
            options.append(str(ex))

    if num.get() != "":
        options.append("-n")
        options.append(num.get())

    if unit.get() != "" and unit.get() != " ":
        options.append("-u")
        options.append(unit.get())
    if colormap.get() != "" and colormap.get() != " ":
        options.append("-c")
        options.append(colormap.get())

    if dem_file.get() != "":
        options.append("--dem")
        options.append(dem_file.get())

    if font_size.get() != "":
        options.append("-s")
        options.append(font_size.get())
    if title_show.get() != 0:
        options.append("--notitle")
    if marker_size.get() != "":
        options.append("--ms")
        options.append(marker_size.get())
    if edge_width.get() != "":
        options.append("--ew")
        options.append(edge_width.get())
    if no_flip.get() == 1:
        options.append("--no-flip")
    if fig_size_width.get() != "" and fig_size_height.get() != "":
        options.append("--figsize")
        options.append(fig_size_height.get())
        options.append(fig_size_width.get())
    if zfirst.get() == 1:
        options.append("--zf")

    if pix_input_xy_x.get() != "" and pix_input_xy_y.get() != "":
        options.append("--yx")
        options.append(pix_input_xy_y.get())
        options.append(pix_input_xy_x.get())
    if pix_input_lalo_la.get() != "" and pix_input_lalo_lo.get() != "":
        options.append("--lalo")
        options.append(pix_input_lalo_la.get())
        options.append(pix_input_lalo_lo.get())

    if ref_pix_input_xy_x.get() != "" and ref_pix_input_xy_y.get() != "":
        options.append("--ref-yx")
        options.append(ref_pix_input_xy_y.get())
        options.append(ref_pix_input_xy_x.get())
    if pix_input_lalo_la.get() != "" and pix_input_lalo_lo.get() != "":
        options.append("--ref-lalo")
        options.append(ref_pix_input_lalo_la.get())
        options.append(ref_pix_input_lalo_lo.get())

    '''if save.get() != 0:
        options.append("--save")
    if output_file.get() != "":
        options.append("-o")

        location_parts = h5_file.get().split("/")
        location = "/".join(location_parts[1:-1])

        options.append("/" + str(location) + "/" + output_file.get())'''

    if show_info.get() == 1:
        file_info = info.hdf5_structure_string(h5_file.get())
        show_file_info(file_info)

    print(options)
    if h5_file.get() != "":
        if ts_view.fig_v is not None:
            ts_view.p1_x = None
            ts_view.p1_y = None
            ts_view.p2_x = None
            ts_view.p2_y = None
            ts_view.p1_scatter_point = None
            ts_view.p2_scatter_point = None
            ts_view.second_plot_axis_visible = False

            ts_view.fig_v.clear()

        ts_view.main(options)
    else:
        print("No file selected")


def reset_plot():
    set_variables_from_attributes()


def set_variables_from_attributes():
    global atr, attributes, data, num_list, starting_upper_lim, ref_dates_list, ref_date_option_menu, excludes_list_box, num_option_menu

    set_sliders()

    unit.set(attributes['UNIT'])
    colormap.set('jet')

    fig_size_width.set("5")
    fig_size_height.set("10")

    marker_size.set("12.0")
    edge_width.set("1.0")
    font_size.set("12")
    no_flip.set(1)
    zfirst.set(1)
    title_show.set(1)

    pix_input_xy_x.set("300")
    pix_input_xy_y.set("300")

    set_epoch_info()


def set_sliders(value=None):
    global starting_lower_lim, starting_upper_lim, y_lim_upper_slider, y_lim_lower_slider, num, attributes, unit

    if value is not None:
        num.set(value)
    else:
        num.set(1)

    data = read_file_data(value)

    max_val = numpy.nanmax(data)
    min_val = numpy.nanmin(data)

    starting_lower_lim = min_val * 5
    starting_upper_lim = max_val * 5

    unit.set(attributes['UNIT'])
    update_sliders(attributes['UNIT'])
    y_lim_upper.set(max_val)
    y_lim_lower.set(min_val)


def set_epoch_info():
    global excludes_list_box, ref_date_option_menu, ref_date, num_list, num_option_menu, the_epoch, num, ref_dates_list

    num_list.pop(0)
    num_list = list(range(len(ref_dates_list)))

    excludes_list_box.delete(0, END)
    ref_date_option_menu['menu'].delete(0, END)
    num_option_menu['menu'].delete(1, END)

    for the_epoch in ref_dates_list:
        ref_date_option_menu.children['menu'].add_command(label=the_epoch,
                                                          command=lambda val=the_epoch: ref_date.set(val))
        if the_epoch is not "All":
            excludes_list_box.insert(END, the_epoch)

    ref_date.set("All")

    for number in num_list:
        num_option_menu.children['menu'].add_command(label=number,
                                                     command=lambda val=number: set_sliders(val))

    num.set(str(len(num_list) - 3))


def read_file_data(epoch=None):
    global atr, attributes, ref_dates_list

    atr = readfile.read_attribute(h5_file.get())

    file_type = atr['FILE_TYPE']

    ref_dates_list = ["All"]

    h5file = h5py.File(h5_file.get(), 'r')
    if file_type in ['HDFEOS']:
        ref_dates_list += h5file.attrs['DATE_TIMESERIES'].split()
    else:
        ref_dates_list = timeseries(h5_file.get()).get_date_list()

    if epoch and epoch is not "All":
        data, attributes = readfile.read(h5_file.get(), datasetName=ref_dates_list[epoch])
    else:
        data, attributes = readfile.read(h5_file.get(), datasetName=ref_dates_list[len(ref_dates_list) - 1])

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

    lalo = ut.coordinate(attributes).box_pixel2geo(data_box)

    formatted_lalo = [str(round(num, 2)) for num in lalo]

    if all_data:
        return formatted_lalo
    else:
        return formatted_lalo[2], formatted_lalo[3]


def compute_xy(lat, lon):
    lat_data = round(float(lat), 4)
    lon_data = round(float(lon), 4)

    data_box = (float(attributes['X_FIRST']), float(attributes['Y_FIRST']), lon_data, lat_data)

    xy = ut.coordinate(attributes).box_geo2pixel(data_box)

    return str(xy[2]), str(xy[3])


def update_pix_input_lalo(x, y, z):
    global update_in_progress

    if update_in_progress:
        return

    update_in_progress = True
    x_in, y_in = pix_input_xy_x.get(), pix_input_xy_y.get()

    lon, lat = compute_lalo(x_in, y_in)

    pix_input_lalo_la.set(lat)
    pix_input_lalo_lo.set(lon)

    update_in_progress = False


def update_pix_input_xy(x, y, z):
    global update_in_progress

    if update_in_progress:
        return

    update_in_progress = True

    lat, lon = pix_input_lalo_la.get(), pix_input_lalo_lo.get()

    x_comp, y_comp = compute_xy(lat, lon)

    pix_input_xy_x.set(x_comp)
    pix_input_xy_y.set(y_comp)

    update_in_progress = False


def update_ref_pix_input_lalo(x, y, z):
    global update_in_progress

    if update_in_progress:
        return

    update_in_progress = True
    x_in, y_in = ref_pix_input_xy_x.get(), ref_pix_input_xy_y.get()

    lon, lat = compute_lalo(x_in, y_in)

    ref_pix_input_lalo_la.set(lat)
    ref_pix_input_lalo_lo.set(lon)

    update_in_progress = False


def update_ref_pix_input_xy(x, y, z):
    global update_in_progress

    if update_in_progress:
        return

    update_in_progress = True

    lat, lon = ref_pix_input_lalo_la.get(), ref_pix_input_lalo_lo.get()

    x_comp, y_comp = compute_xy(lat, lon)

    ref_pix_input_xy_x.set(x_comp)
    ref_pix_input_xy_y.set(y_comp)

    update_in_progress = False


def validate_numbers(action, index, value_if_allowed,
                     prior_value, text, validation_type, trigger_type, widget_name):

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
    global canvas, frame, h5_file, h5_file_short, pick_h5_file_button, mask_file, mask_short, \
        pick_mask_file_button, starting_upper_lim, y_lim_upper, y_lim_upper_slider, y_lim_lower, y_lim_lower_slider, unit, \
        colormap, projection, lr_flip, ud_flip, wrap, opposite, transparency, show_info, dem_file, dem_short, \
        pick_dem_file_button, shading, countours, countour_smoothing, countour_step, pix_input_xy_x, pix_input_xy_y, pix_input_lalo_la, \
        pix_input_lalo_lo, ref_pix_input_xy_x, ref_pix_input_xy_y, ref_pix_input_lalo_la, ref_pix_input_lalo_lo, ref_x, ref_y, \
        ref_lat, ref_lon, ref_color, ref_sym, font_size, title_show, marker_size, edge_width, no_flip, zfirst, title_show, tick_show, \
        title_in, title, fig_size_width, fig_size_height, fig_ext, fig_num, fig_w_space, fig_h_space, coords, coastline, resolution, \
        lalo_label, lalo_step, scalebar_distance, scalebar_lat, scalebar_lon, show_scalebar, save, output_file, ref_date_option_menu, \
        ref_date, excludes_list_box, num, num_option_menu, settings_file, settings_to_variable, settings_file_button

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
                                 command=lambda: pick_file(), padx=15)
    selected_ts_file_label = Label(pick_h5_file_frame, textvariable=h5_file_short)

    '''     Frames, Text Variables, and Widgets for selection of the maskConnComp.h5 file to add a mask to the data.     '''
    pick_mask_file_frame = Frame(frame)

    mask_file = StringVar()
    mask_file.trace('w', callback=set_mask_short)
    mask_short = StringVar()
    mask_short.set("No File Selected")

    pick_mask_file_button = Button(pick_mask_file_frame, text='Select Mask File', width=15,
                                   command=lambda: pick_mask(), padx=15)
    selected_mask_file_label = Label(pick_mask_file_frame, textvariable=mask_short)

    '''     Frames, Text Variables, and Widgets for selection of the dem.h5 file to add a mask to the data.     '''
    pick_dem_file_frame = Frame(frame)

    dem_file = StringVar()
    dem_file.trace('w', callback=set_dem_short)
    dem_short = StringVar()
    dem_short.set("No File Selected")

    pick_dem_file_button = Button(pick_dem_file_frame, text='Select Topography File', width=15,
                                  command=lambda: pick_dem(), padx=15)
    selected_dem_file_label = Label(pick_dem_file_frame, textvariable=dem_short)

    '''     WIDGETS FOR SHOWING EPOCHS AND EXLUDE DATES     '''

    num_labels_frame = Frame(frame)

    num_option_menu_label = Label(num_labels_frame, text="Plot Number", width=10, anchor='w')

    num_frame = Frame(frame)

    num = StringVar()
    num.set(" ")
    num_option_menu = OptionMenu(num_frame, num, *num_list)
    num_option_menu.config(width=10)

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


    '''     WIDGETS FOR UNIT, COLORMAP, AND PROJECTION    '''
    unit_cmap_projection_labels_frame = Frame(frame)
    unit_label = Label(unit_cmap_projection_labels_frame, text="Unit", width=16, anchor='w')
    colormap_label = Label(unit_cmap_projection_labels_frame, text="Colormap", width=14, anchor='w')

    unit_cmap_projection_frame = Frame(frame)

    unit = StringVar()
    unit.set(" ")
    unit_option_menu = OptionMenu(unit_cmap_projection_frame, unit, *unit_options, command=scale_sliders)
    unit_option_menu.config(width=14)

    colormap = StringVar()
    colormap.set(" ")
    colormap_option_menu = OptionMenu(unit_cmap_projection_frame, colormap, *colormaps)
    colormap_option_menu.config(width=14)


    '''     WIDGETS FO DIGURE SIZE      '''
    fig_size_frame = Frame(frame)

    fig_size_label = Label(fig_size_frame, text="Fig Size")

    fig_size_width = StringVar()
    fig_size_width_label = Label(fig_size_frame, text="Width: ")
    fig_size_width_entry = Entry(fig_size_frame, textvariable=fig_size_width, width=6, validate='key', validatecommand=vcmd_num)

    fig_size_height = StringVar()
    fig_size_height_label = Label(fig_size_frame, text="Length: ")
    fig_size_height_entry = Entry(fig_size_frame, textvariable=fig_size_height, width=6, validate='key', validatecommand=vcmd_num)

    '''     WIDGETS FOR FONT SIZE AND FIGURE DPI'''
    font_title_frame = Frame(frame)

    font_size = StringVar()
    font_size_label = Label(font_title_frame, text="Font Size:    ")
    font_size_entry = Entry(font_title_frame, textvariable=font_size, width=6, validate='key', validatecommand=vcmd_num)

    title_show = IntVar()
    show_title_label = Label(font_title_frame, text="Show Title:    ")
    show_title_checkbutton = Checkbutton(font_title_frame, text="Show Title", variable=title_show)


    '''     WIDGETS FOR SHOWING AXIS, COLORBAR, TITLE, AND AXIS TICKS       '''
    flip_zfirst_frae = Frame(frame)

    no_flip = IntVar()
    no_flip_checkbutton = Checkbutton(flip_zfirst_frae, text="No Flip", variable=no_flip)

    zfirst = IntVar()
    zfirst_checkbutton = Checkbutton(flip_zfirst_frae, text="Zero First", variable=zfirst)


    '''     WIDGETS FOR NUMBER OF ROWS AND NUMBER OF COLUMNS      '''
    mkrsize_edgewidth_frame = Frame(frame)

    marker_size = StringVar()
    marker_size_label = Label(mkrsize_edgewidth_frame, text="Marker Size: ")
    marker_size_entry = Entry(mkrsize_edgewidth_frame, textvariable=marker_size, width=6, validate='key', validatecommand=vcmd_num)

    edge_width = StringVar()
    edge_width_label = Label(mkrsize_edgewidth_frame, text="Edge Width: ")
    edge_width_entry = Entry(mkrsize_edgewidth_frame, textvariable=edge_width, width=6, validate='key', validatecommand=vcmd_num)


    '''     WIDGETS FOR SHOWING EPOCHS AND EXLUDE DATES     '''

    epoch_labels_frame = Frame(frame)

    epoch_option_menu_label = Label(epoch_labels_frame, text="Reference Date", width=15, anchor='w')
    exclude_date_label = Label(epoch_labels_frame, text="Exclude Dates", width=15, anchor='w')

    epoch_frame = Frame(frame)

    ref_date = StringVar()
    ref_date.set(" ")
    ref_date_option_menu = OptionMenu(epoch_frame, ref_date, *ref_dates_list)
    ref_date_option_menu.config(width=15)

    excludes_list_box = Listbox(epoch_frame, selectmode=MULTIPLE, height=5, width=15)

    '''     WIDGETS FOR SHOWING INFO SCREEN'''
    show_info = IntVar()
    show_info_checkbutton = Checkbutton(frame, text="Show File Info", variable=show_info)





    '''
            |-----------------------------------------------------------------------------------------------------|
            |                                                                                                     |
            |                                 WIDGETS TO CONTROL SUBSET OPTIONS                                   |
            |                                                                                                     |
            |-----------------------------------------------------------------------------------------------------|

    '''

    '''     SUBSET OPTIONS WIDGETS      '''
    pixel_input_label = Label(frame, text="PIXEL INPUT", anchor=W)


    '''     WIDGETS FOR SUBSET X-VALUES'''
    pix_input_xy_frame = Frame(frame)

    pix_input_xy_x = StringVar()
    pix_input_xy_x.trace('w', callback=update_pix_input_lalo)
    pix_input_xy_x_label = Label(pix_input_xy_frame, text="X:     ")
    pix_input_xy_x_entry = Entry(pix_input_xy_frame, textvariable=pix_input_xy_x, width=6, validate='key', validatecommand=vcmd_num)

    pix_input_xy_y = StringVar()
    pix_input_xy_y.trace('w', callback=update_pix_input_lalo)
    pix_input_xy_y_label = Label(pix_input_xy_frame, text="Y:    ")
    pix_input_xy_y_entry = Entry(pix_input_xy_frame, textvariable=pix_input_xy_y, width=6, validate='key', validatecommand=vcmd_num)


    '''     WIDGETS FOR SUBSET LAT-VALUES       '''
    pix_input_lalo_frame = Frame(frame)

    pix_input_lalo_la = StringVar()
    pix_input_lalo_la.trace('w', callback=update_pix_input_xy)
    pix_input_lalo_la_label = Label(pix_input_lalo_frame, text="Lat:  ")
    pix_input_lalo_la_entry = Entry(pix_input_lalo_frame, textvariable=pix_input_lalo_la, width=6, validate='key', validatecommand=vcmd_num)

    pix_input_lalo_lo = StringVar()
    pix_input_lalo_lo.trace('w', callback=update_pix_input_xy)
    pix_input_lalo_lo_label = Label(pix_input_lalo_frame, text="Lon:   ")
    pix_input_lalo_lo_entry = Entry(pix_input_lalo_frame, textvariable=pix_input_lalo_lo, width=6, validate='key', validatecommand=vcmd_num)


    '''     SUBSET OPTIONS WIDGETS      '''
    ref_pixel_input_label = Label(frame, text="REFERENCE PIXEL INPUT", anchor=W)

    '''     WIDGETS FOR SUBSET X-VALUES'''
    ref_pix_input_xy_frame = Frame(frame)

    ref_pix_input_xy_x = StringVar()
    ref_pix_input_xy_x.trace('w', callback=update_ref_pix_input_lalo)
    ref_pix_input_xy_x_label = Label(ref_pix_input_xy_frame, text="X:     ")
    ref_pix_input_xy_x_entry = Entry(ref_pix_input_xy_frame, textvariable=ref_pix_input_xy_x, width=6, validate='key', validatecommand=vcmd_num)

    ref_pix_input_xy_y = StringVar()
    ref_pix_input_xy_y.trace('w', callback=update_ref_pix_input_lalo)
    ref_pix_input_xy_y_label = Label(ref_pix_input_xy_frame, text="Y:      ")
    ref_pix_input_xy_y_entry = Entry(ref_pix_input_xy_frame, textvariable=ref_pix_input_xy_y, width=6, validate='key', validatecommand=vcmd_num)

    '''     WIDGETS FOR SUBSET LAT-VALUES       '''
    ref_pix_input_lalo_frame = Frame(frame)

    ref_pix_input_lalo_la = StringVar()
    ref_pix_input_lalo_la.trace('w', callback=update_ref_pix_input_xy)
    ref_pix_input_lalo_la_label = Label(ref_pix_input_lalo_frame, text="Lat:   ")
    ref_pix_input_lalo_la_entry = Entry(ref_pix_input_lalo_frame, textvariable=ref_pix_input_lalo_la, width=6, validate='key', validatecommand=vcmd_num)

    ref_pix_input_lalo_lo = StringVar()
    ref_pix_input_lalo_lo.trace('w', callback=update_ref_pix_input_xy)
    ref_pix_input_lalo_lo_label = Label(ref_pix_input_lalo_frame, text="Lon:   ")
    ref_pix_input_lalo_lo_entry = Entry(ref_pix_input_lalo_frame, textvariable=ref_pix_input_lalo_lo, width=6, validate='key', validatecommand=vcmd_num)





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

    pick_dem_file_frame.pack(fill=X)
    pick_dem_file_button.pack(side=LEFT, anchor='w', pady=5, padx=(10, 20))
    selected_dem_file_label.pack(side=LEFT, fill=X)

    num_labels_frame.pack(anchor='w', fill=X, pady=5)
    num_option_menu_label.pack(side=LEFT, padx=(12, 20))

    num_frame.pack(anchor='w', fill=X)
    num_option_menu.pack(side=LEFT, anchor='n', pady=5, padx=(10, 20))


    pixel_input_label.pack(anchor='w', fill=X, pady=(15, 0), padx=10)

    pix_input_xy_frame.pack(anchor='w', fill=X, pady=10, padx=10)

    pix_input_xy_x_label.pack(side=LEFT, padx=(0, 5))
    pix_input_xy_x_entry.pack(side=LEFT, padx=(0, 10))
    pix_input_xy_y_label.pack(side=LEFT, padx=(0, 5))
    pix_input_xy_y_entry.pack(side=LEFT, padx=(0, 10))

    pix_input_lalo_frame.pack(anchor='w', fill=X, pady=10, padx=10)

    pix_input_lalo_la_label.pack(side=LEFT, padx=(0, 5))
    pix_input_lalo_la_entry.pack(side=LEFT, padx=(0, 10))
    pix_input_lalo_lo_label.pack(side=LEFT, padx=(0, 5))
    pix_input_lalo_lo_entry.pack(side=LEFT, padx=(0, 10))

    ref_pixel_input_label.pack(anchor='w', fill=X, pady=(15, 0), padx=10)

    ref_pix_input_xy_frame.pack(anchor='w', fill=X, pady=10, padx=10)

    ref_pix_input_xy_x_label.pack(side=LEFT, padx=(0, 5))
    ref_pix_input_xy_x_entry.pack(side=LEFT, padx=(0, 10))
    ref_pix_input_xy_y_label.pack(side=LEFT, padx=(0, 5))
    ref_pix_input_xy_y_entry.pack(side=LEFT, padx=(0, 10))

    ref_pix_input_lalo_frame.pack(anchor='w', fill=X, pady=10, padx=10)

    ref_pix_input_lalo_la_label.pack(side=LEFT, padx=(0, 5))
    ref_pix_input_lalo_la_entry.pack(side=LEFT, padx=(0, 10))
    ref_pix_input_lalo_lo_label.pack(side=LEFT, padx=(0, 5))
    ref_pix_input_lalo_lo_entry.pack(side=LEFT, padx=(0, 10))



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

    unit_cmap_projection_labels_frame.pack(anchor='w', fill=X, pady=(10, 0), padx=10)
    unit_label.pack(side=LEFT, padx=(0, 20))
    colormap_label.pack(side=LEFT, padx=(0, 20))

    unit_cmap_projection_frame.pack(anchor='w', fill=X, pady=(10, 5), padx=10)
    unit_option_menu.pack(side=LEFT, padx=(0, 10))
    colormap_option_menu.pack(side=LEFT, padx=(0, 10))

    fig_size_frame.pack(anchor='w', fill=X, pady=10, padx=10)
    fig_size_label.pack(side=LEFT, padx=(0, 10))
    fig_size_width_label.pack(side=LEFT, padx=(0, 5))
    fig_size_width_entry.pack(side=LEFT, padx=(0, 10))
    fig_size_height_label.pack(side=LEFT, padx=(0, 5))
    fig_size_height_entry.pack(side=LEFT, padx=(0, 10))

    font_title_frame.pack(anchor='w', fill=X, padx=10, pady=10)
    font_size_label.pack(side=LEFT, padx=(0, 5))
    font_size_entry.pack(side=LEFT, padx=(0, 30))
    show_title_checkbutton.pack(side=LEFT)

    flip_zfirst_frae.pack(anchor='w', fill=X, pady=10, padx=10)
    no_flip_checkbutton.pack(side=LEFT, padx=(0, 25))
    zfirst_checkbutton.pack(side=LEFT, padx=(0, 25))

    mkrsize_edgewidth_frame.pack(anchor='w', fill=X, padx=10, pady=10)
    marker_size_label.pack(side=LEFT, padx=(0, 0))
    marker_size_entry.pack(side=LEFT, padx=(0, 19))
    edge_width_label.pack(side=LEFT, padx=(0, 0))
    edge_width_entry.pack(side=LEFT)

    epoch_labels_frame.pack(anchor='w', fill=X, pady=5)
    epoch_option_menu_label.pack(side=LEFT, padx=(12, 20))
    exclude_date_label.pack(side=LEFT)

    epoch_frame.pack(anchor='w', fill=X)
    ref_date_option_menu.pack(side=LEFT, anchor='n', pady=5, padx=(10, 20))
    excludes_list_box.pack(side=LEFT, pady=5)

    show_info_checkbutton.pack(anchor='center', pady=10)

    save_settings_button.pack(side=TOP)

    space = Frame(frame)
    space.config(height=50)
    space.pack(side=LEFT)

    if inps.file:
        on_file_selection(inps.file)

    settings_to_variable = {

        "input.file": h5_file,
        "input.num": num,
        "input.mask": mask_file,
        "input.dem": dem_file,

        "output.save": save,
        "output.outfile": output_file,

        "display.min": y_lim_lower,
        "display.max": y_lim_upper,
        "display.unit": unit,
        "display.ref_date": ref_date,
        "display.colormap": colormap,
        "display.no_flip": no_flip,
        #"display.alpha": transparency,
        "display.fig_length": fig_size_height,
        "display.fig_width": fig_size_width,
        "display.no_title": title_show,
        "display.font_size": font_size,
        "display.marker_size": marker_size,
        "display.edge_width": edge_width,

        "pixel.x": pix_input_xy_x,
        "pixel.y": pix_input_xy_y,
        "pixel.ref_x": ref_pix_input_xy_x,
        "pixel.ref_y": ref_pix_input_xy_y

    }

    if inps.settings:
        settings_file.set(inps.settings)
        parse_settings()

    mainloop()


if __name__ == '__main__':
    inps = cmdLineParse(sys.argv[1:])
    main()
