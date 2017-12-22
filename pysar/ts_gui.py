from Tkinter import *
import matplotlib
matplotlib.use('TkAgg')
import tkFileDialog as filedialog
import tsviewer as ts_view


def pick_file():
    filename = filedialog.askopenfilename(initialdir="/User/Joshua/", title="Select file",
                                          filetypes=(("jpeg files", "*.h5"), ("all files", "*.*")))
    root.filename = filename
    timeseries_file.set(root.filename)
    file_short.set(filename.split("/")[-1])
    return root.filename


def pick_mask():
    filename = filedialog.askopenfilename(initialdir="/User/Joshua/", title="Select file",
                                          filetypes=(("jpeg files", "*.h5"), ("all files", "*.*")))
    root.filename = filename
    mask.set(root.filename)
    mask_short.set(filename.split("/")[-1])
    return root.filename


def pick_dem():
    filename = filedialog.askopenfilename(initialdir="/User/Joshua/", title="Select file",
                                          filetypes=(("jpeg files", "*.h5"), ("all files", "*.*")))
    root.filename = filename
    dem.set(root.filename)
    dem_short.set(filename.split("/")[-1])
    return root.filename


def show_plot():

    options = [timeseries_file.get(), "-m", mask.get(), "--dem", dem.get(), "--lalo", start_lat_input.get(),
               start_lon_input.get(), "--ref-lalo", ref_lat_input.get(), ref_lon_input.get(),
               "--ylim", str(y_lim_lower.get()), str(y_lim_upper.get()), "-u", unit.get(), "-c", colormap.get()]

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


root = Tk()
root.minsize(width=350, height=550)
root.maxsize(width=350, height=550)
root.resizable(width=False, height=False)

colormaps = ['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2',
             'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r',
             'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r',
             'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r',
             'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Vega10', 'Vega10_r',
             'Vega20', 'Vega20_r', 'Vega20b', 'Vega20b_r', 'Vega20c', 'Vega20c_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr',
             'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr',
             'bwr_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth',
             'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r',
             'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r',
             'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r',
             'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spectral', 'spectral_r', 'spring', 'spring_r', 'summer',
             'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'viridis', 'viridis_r',
             'winter', 'winter_r']


'''     Frames, Text Variables, and Widgets for selection of the timeseries.h5 file to plot data from.     '''
pick_timeseries_file_frame = Frame(root)

timeseries_file = StringVar()
file_short = StringVar()

pick_ts_file_button = Button(pick_timeseries_file_frame, text='Select Timeseries File', anchor='w', width=15, command=lambda: pick_file())
selected_ts_file_label = Label(pick_timeseries_file_frame, textvariable=file_short)


'''     Frames, Text Variables, and Widgets for selection of the mask.h5 file to add a mask to the ata.     '''
pick_mask_file_frame = Frame(root)

mask = StringVar()
mask_short = StringVar()

pick_mask_file_button = Button(pick_mask_file_frame, text='Select Mask File', anchor='w', width=15, command=lambda: pick_mask())
selected_mask_file_label = Label(pick_mask_file_frame, textvariable=mask_short)


'''     Frames, Text Variables, and Widgets for selection of the topography dem.h5 file to add topography to the data.     '''
pick_dem_file_frame = Frame(root)

dem = StringVar()
dem_short = StringVar()

pick_dem_file_button = Button(pick_dem_file_frame, text='Select Topography File', anchor='w', width=15, command=lambda: pick_dem())
selected_dem_file_label = Label(pick_dem_file_frame, textvariable=dem_short)


'''     Frames, Text Variables, and Widgets for setting latitude and longitude data for plotting.     '''
lat_lon_frame = Frame(root)

'''     Frames, Text Varieables, and Widgets for start lat/long coordinates     '''
start_lat_lon_frame = Frame(lat_lon_frame)

'''     Frames, Text Variables, and Widgtes for start lat coordinate entry     '''
start_lat_frame = Frame(start_lat_lon_frame)

start_lat_input = StringVar()
start_lat_input.set("-0.4606")

start_lat_label = Label(start_lat_frame, text="Starting Latitude", width=17)
start_lat_entry = Entry(start_lat_frame, textvariable=start_lat_input, width=17)

'''     Frames, Text Variables, and Widgtes for start long coordinate entry     '''
start_lon_frame = Frame(start_lat_lon_frame)

start_lon_input = StringVar()
start_lon_input.set("-91.0752")

start_lon_label = Label(start_lon_frame, text="Starting Longitude", width=17)
start_lon_entry = Entry(start_lon_frame, textvariable=start_lon_input, width=17)

'''     Frames, Text Variables, and Widgets for reference lat/long coordinates     '''
ref_lat_lon_frame = Frame(lat_lon_frame, width=10)

'''     Frames, Text Variables, and Widgtes for reference lat coordinate entry     '''
ref_lat_frame = Frame(ref_lat_lon_frame)

ref_lat_input = StringVar()
ref_lat_input.set("-0.4606")

ref_lat_label = Label(ref_lat_frame, text="Reference Latitude", width=17)
ref_lat_entry = Entry(ref_lat_frame, textvariable=ref_lat_input, width=17)

'''     Frames, Text Variables, and Widgtes for reference long coordinate entry     '''
ref_lon_frame = Frame(ref_lat_lon_frame)

ref_lon_input = StringVar()
ref_lon_input.set("-91.0752")

ref_lon_label = Label(ref_lon_frame, text="Reference Longitude", width=17)
ref_lon_entry = Entry(ref_lon_frame, textvariable=ref_lon_input, width=17)


'''     Frames, Text Variables, and Widgets for setting y-lim      '''
y_lim_frame = Frame(root)
y_lim_upper_frame = Frame(y_lim_frame)

y_lim_upper = IntVar()
y_lim_upper.set(20)

y_lim_upper_label = Label(y_lim_upper_frame, text="Upper Y-Lim", width=8)
y_lim_upper_slider = Scale(y_lim_upper_frame, from_= -40, to= 40, orient=HORIZONTAL, length=150, variable=y_lim_upper, showvalue=0)
y_lim_upper_entry = Entry(y_lim_upper_frame, textvariable=y_lim_upper, width=6)

y_lim_lower_frame = Frame(y_lim_frame)

y_lim_lower = IntVar()
y_lim_lower.set(-20)

y_lim_lower_label = Label(y_lim_lower_frame, text="Lower Y-Lim", width=8)
y_lim_lower_slider = Scale(y_lim_lower_frame, from_= -40, to= 40, orient=HORIZONTAL, length=150, variable=y_lim_lower, showvalue=0)
y_lim_lower_entry = Entry(y_lim_lower_frame, textvariable=y_lim_lower, width=6)

'''     Frames, Text Variables, and Widgets for setting extraneous widgets      '''
unit_cmap_frame = Frame(root)

unit = StringVar()
unit.set("cm")
unit_option_menu = apply(OptionMenu, (unit_cmap_frame, unit) + tuple(["cm", "m", "dm", "km"]))
unit_option_menu.config(width=6)

colormap = StringVar()
colormap_option_menu = apply(OptionMenu, (unit_cmap_frame, colormap) + tuple(colormaps))
colormap_option_menu.config(width=10)
colormap.set('hsv')

ref_date = StringVar()
ref_date_option_menu = apply(OptionMenu, (unit_cmap_frame, ref_date) + tuple(["09-11-2009", "09-12-2009", "09-13-2009"]))
ref_date_option_menu.config(width=12)
ref_date.set('09-11-2009')

submit_button = Button(root, text="Show Plot", command=lambda: show_plot())




pick_timeseries_file_frame.pack(anchor='w', fill=X, pady=(10, 5), padx=10)
pick_ts_file_button.pack(side=LEFT, anchor='w', padx=(0, 20))
selected_ts_file_label.pack(side=LEFT, fill=X)

pick_mask_file_frame.pack(anchor='w', fill=X)
pick_mask_file_button.pack(side=LEFT, anchor='w', pady=5, padx=10)
selected_mask_file_label.pack(side=LEFT, fill=X)

pick_dem_file_frame.pack(anchor='w', fill=X)
pick_dem_file_button.pack(side=LEFT, anchor='w', pady=5, padx=10)
selected_dem_file_label.pack(side=LEFT, fill=X)


lat_lon_frame.pack(anchor='w', fill=X, pady=10, padx=10)

start_lat_lon_frame.pack(side=TOP, pady=(0, 10), fill=X)

start_lat_frame.pack(side=LEFT)
start_lat_label.pack(side=TOP)
start_lat_entry.pack(side=TOP)


start_lon_frame.pack(side=LEFT)
start_lon_label.pack(side=TOP)
start_lon_entry.pack(side=TOP)

ref_lat_lon_frame.pack(side=TOP, fill=X)

ref_lat_frame.pack(side=LEFT)
ref_lat_label.pack(side=TOP)
ref_lat_entry.pack(side=TOP)


ref_lon_frame.pack(side=LEFT)
ref_lon_label.pack(side=TOP)
ref_lon_entry.pack(side=TOP)


y_lim_frame.pack(anchor='w', fill=X, pady=10, padx=10)

y_lim_upper_frame.pack(side=TOP, fill=X, pady=(0, 10))
y_lim_upper_label.pack(side=LEFT)
y_lim_upper_slider.pack(side=LEFT, padx=10)
y_lim_upper_entry.pack(side=LEFT)

y_lim_lower_frame.pack(side=TOP, fill=X)
y_lim_lower_label.pack(side=LEFT)
y_lim_lower_slider.pack(side=LEFT, padx=10)
y_lim_lower_entry.pack(side=LEFT)

unit_cmap_frame.pack(anchor='w', fill=X, pady=10, padx=10)
unit_option_menu.pack(side=LEFT, padx=(0, 10))
colormap_option_menu.pack(side=LEFT, padx=(0, 10))
ref_date_option_menu.pack(side=LEFT)

submit_button.pack(anchor='center', pady=20)

mainloop()