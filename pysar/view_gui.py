from Tkinter import *

import h5py
import matplotlib
matplotlib.use('TkAgg')
import tkFileDialog as filedialog
import view as view
import info

def pick_file():
    if h5_file.get() == "":
        filename = filedialog.askopenfilename(initialdir="/User/Joshua/", title="Select file",
                                              filetypes=(("jpeg files", "*.h5"), ("all files", "*.*")))
        frame.filename = filename
        h5_file.set(frame.filename)
        h5_file_short.set(filename.split("/")[-1])
        pick_h5_file_button.config(text="Cancel")
        return frame.filename
    else:
        h5_file.set("")
        h5_file_short.set("No File Selected")
        pick_h5_file_button.config(text="Select .h5 File")


def pick_mask():
    if mask_file.get() == "":
        filename = filedialog.askopenfilename(initialdir="/User/Joshua/", title="Select file",
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


def pick_dem():
    if dem_file.get() == "":
        filename = filedialog.askopenfilename(initialdir="/User/Joshua/", title="Select file",
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


def on_configure(event):
    # update scrollregion after starting 'mainloop'
    # when all widgets are in canvas
    canvas.configure(scrollregion=canvas.bbox('all'))


root = Tk()
root.minsize(width=375, height=750)
root.maxsize(width=375, height=750)
root.resizable(width=False, height=False)

canvas = Canvas(root, width=355)
canvas.pack(side=LEFT, fill=Y)

scrollbar = Scrollbar(root, command=canvas.yview)
scrollbar.pack(side=LEFT, fill='y')

canvas.configure(yscrollcommand = scrollbar.set)
canvas.bind('<Configure>', on_configure)

frame = Frame(canvas)
canvas.create_window((0,0), window=frame, anchor='nw')

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

projections = ["cea", "mbtfpq", "aeqd", "sinu", "poly", "moerc", "gnom", "moll", "lcc", "tmerc", "nplaea", "gall",
               "npaeqd", "mill", "merc", "stere", "eqdc", "rotpole", "cyl", "npstere", "spstere", "hammer", "geos",
               "nsper", "eck4", "aea", "kav7", "spaeqd", "ortho", "class", "vandg", "laea", "splaea", "robin"]



'''     Frames, Text Variables, and Widgets for selection of the timeseries.h5 file to plot data from.     '''
pick_h5_file_frame = Frame(frame)

h5_file = StringVar()
h5_file_short = StringVar()
h5_file_short.set("No File Selected")

pick_h5_file_button = Button(pick_h5_file_frame, text='Select Timeseries File', anchor='w', width=15, command=lambda: pick_file())
selected_ts_file_label = Label(pick_h5_file_frame, textvariable=h5_file_short)


'''     Frames, Text Variables, and Widgets for selection of the mask.h5 file to add a mask to the ata.     '''
pick_mask_file_frame = Frame(frame)

mask_file = StringVar()
mask_short = StringVar()
mask_short.set("No File Selected")

pick_mask_file_button = Button(pick_mask_file_frame, text='Select Mask File', anchor='w', width=15, command=lambda: pick_mask())
selected_mask_file_label = Label(pick_mask_file_frame, textvariable=mask_short)


'''     Frames, Text Variables, and Widgets for setting y-lim      '''
y_lim_frame = Frame(frame)
y_lim_upper_frame = Frame(y_lim_frame)

y_lim_upper = IntVar()
y_lim_upper.set(20)

y_lim_upper_label = Label(y_lim_upper_frame, text="Maximum", width=8)
y_lim_upper_slider = Scale(y_lim_upper_frame, from_= -40, to= 40, orient=HORIZONTAL, length=150, variable=y_lim_upper, showvalue=0)
y_lim_upper_entry = Entry(y_lim_upper_frame, textvariable=y_lim_upper, width=6)

y_lim_lower_frame = Frame(y_lim_frame)

y_lim_lower = IntVar()
y_lim_lower.set(-20)

y_lim_lower_label = Label(y_lim_lower_frame, text="Minimum", width=8)
y_lim_lower_slider = Scale(y_lim_lower_frame, from_= -40, to= 40, orient=HORIZONTAL, length=150, variable=y_lim_lower, showvalue=0)
y_lim_lower_entry = Entry(y_lim_lower_frame, textvariable=y_lim_lower, width=6)


'''     Frames, Text Variables, and Widgets for setting extraneous properties      '''
unit_cmap_projection_frame = Frame(frame)

unit = StringVar()
unit.set("cm")
unit_option_menu = apply(OptionMenu, (unit_cmap_projection_frame, unit) + tuple(["cm", "m", "dm", "km"]))
unit_option_menu.config(width=6)

colormap = StringVar()
colormap_option_menu = apply(OptionMenu, (unit_cmap_projection_frame, colormap) + tuple(colormaps))
colormap_option_menu.config(width=10)
colormap.set('hsv')

projection = StringVar()
projection_option_menu = apply(OptionMenu, (unit_cmap_projection_frame, projection) + tuple(projections))
projection_option_menu.config(width=12)
projection.set("cea")


flip_frame = Frame(frame)

lr_flip = IntVar()
lr_flip_checkbutton = Checkbutton(flip_frame, text="Flip LR", variable=lr_flip)

ud_flip = IntVar()
ud_flip_checkbutton = Checkbutton(flip_frame, text="Flip UD", variable=ud_flip)

wrap = IntVar()
wrap_checkbutton = Checkbutton(flip_frame, text="Wrap", variable=wrap)

opposite = IntVar()
opposite_checkbutton = Checkbutton(flip_frame, text="Opposite", variable=opposite)


transparency_frame = Frame(frame)
transparency_label = Label(transparency_frame, text="Alpha", width=8)
transparency_slider = Scale(transparency_frame, from_= -40, to= 40, orient=HORIZONTAL, length=150, variable=y_lim_upper, showvalue=0)
transparency_entry = Entry(transparency_frame, textvariable=y_lim_upper, width=6)



'''     Frames, Text Variables, and Widgets for selection of the topography dem.h5 file to add topography to the data.     '''
pick_dem_file_frame = Frame(frame)

dem_file = StringVar()
dem_short = StringVar()
dem_short.set("No File Selected")

pick_dem_file_button = Button(pick_dem_file_frame, text='Select Topography File', anchor='w', width=15, command=lambda: pick_dem())
selected_dem_file_label = Label(pick_dem_file_frame, textvariable=dem_short)

dem_options_frame = Frame(frame)

shading = IntVar()
dem_shading_checkbutton = Checkbutton(dem_options_frame, text="Show Shaded Relief", variable=shading)

countours = IntVar()
dem_countours_checkbutton = Checkbutton(dem_options_frame, text="Show Countour Lines", variable=countours)

dem_countour_options = Frame(frame)

dem_countour_smoothing_frame = Frame(dem_countour_options, width=15)

countour_smoothing = StringVar()
dem_countour_smoothing_label = Label(dem_countour_smoothing_frame, text="Contour Smoothing: ", anchor='c', width=15)
dem_countour_smoothing_entry = Entry(dem_countour_smoothing_frame, textvariable=countour_smoothing, width=6)

dem_countour_step_frame = Frame(dem_countour_options, width=15)

countour_step = StringVar()
dem_countour_step_label = Label(dem_countour_step_frame, text="Countour Step: ", anchor='c', width=15)
dem_countour_step_entry = Entry(dem_countour_step_frame, textvariable=countour_step, width=6)


subset_label = Label(frame, text="SUBSET DATA")

subset_x_frame = Frame(frame)

subset_x_from = StringVar()
subset_x_from_label = Label(subset_x_frame, text="X      From: ")
subset_x_from_entry = Entry(subset_x_frame, textvariable=subset_x_from, width=6)

subset_x_to = StringVar()
subset_x_to_label = Label(subset_x_frame, text="To: ")
subset_x_to_entry = Entry(subset_x_frame, textvariable=subset_x_to, width=6)

subset_y_frame = Frame(frame)

subset_y_from = StringVar()
subset_y_from_label = Label(subset_y_frame, text="Y      From: ")
subset_y_from_entry = Entry(subset_y_frame, textvariable=subset_y_from, width=6)

subset_y_to = StringVar()
subset_y_to_label = Label(subset_y_frame, text="To: ")
subset_y_to_entry = Entry(subset_y_frame, textvariable=subset_y_to, width=6)

subset_lat_frame = Frame(frame)

subset_lat_from = StringVar()
subset_lat_from_label = Label(subset_lat_frame, text="Lat      From: ")
subset_lat_from_entry = Entry(subset_lat_frame, textvariable=subset_lat_from, width=6)

subset_lat_to = StringVar()
subset_lat_to_label = Label(subset_lat_frame, text="To: ")
subset_lat_to_entry = Entry(subset_lat_frame, textvariable=subset_lat_to, width=6)

subset_lon_frame = Frame(frame)

subset_lon_from = StringVar()
subset_lon_from_label = Label(subset_lon_frame, text="Lon      From: ")
subset_lon_from_entry = Entry(subset_lon_frame, textvariable=subset_lon_from, width=6)

subset_lon_to = StringVar()
subset_lon_to_label = Label(subset_lon_frame, text="To: ")
subset_lon_to_entry = Entry(subset_lon_frame, textvariable=subset_lon_to, width=6)





pick_h5_file_frame.pack(anchor='w', fill=X, pady=(10, 5), padx=10)
pick_h5_file_button.pack(side=LEFT, anchor='w', padx=(0, 20))
selected_ts_file_label.pack(side=LEFT, fill=X)

pick_mask_file_frame.pack(anchor='w', fill=X)
pick_mask_file_button.pack(side=LEFT, anchor='w', pady=5, padx=(10, 20))
selected_mask_file_label.pack(side=LEFT, fill=X)

y_lim_frame.pack(anchor='w', fill=X, pady=(35, 10), padx=10)

y_lim_upper_frame.pack(side=TOP, fill=X, pady=(0, 10))
y_lim_upper_label.pack(side=LEFT)
y_lim_upper_slider.pack(side=LEFT, padx=10)
y_lim_upper_entry.pack(side=LEFT)

y_lim_lower_frame.pack(side=TOP, fill=X)
y_lim_lower_label.pack(side=LEFT)
y_lim_lower_slider.pack(side=LEFT, padx=10)
y_lim_lower_entry.pack(side=LEFT)

unit_cmap_projection_frame.pack(anchor='w', fill=X, pady=10, padx=10)
unit_option_menu.pack(side=LEFT, padx=(0, 10))
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

pick_dem_file_frame.pack(anchor='w', fill=X, pady=(35, 10))
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

subset_label.pack(anchor='w', fill=X, pady=(25, 10), padx=10)

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


mainloop()