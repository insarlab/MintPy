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
        root.filename = filename
        h5_file.set(root.filename)
        h5_file_short.set(filename.split("/")[-1])
        pick_h5_file_button.config(text="Cancel")
        return root.filename
    else:
        h5_file.set("")
        h5_file_short.set("No File Selected")
        pick_h5_file_button.config(text="Select .h5 File")


def pick_mask():
    if mask_file.get() == "":
        filename = filedialog.askopenfilename(initialdir="/User/Joshua/", title="Select file",
                                              filetypes=(("jpeg files", "*.h5"), ("all files", "*.*")))
        root.filename = filename
        mask_file.set(root.filename)
        mask_short.set(filename.split("/")[-1])
        pick_mask_file_button.config(text="Cancel")
        return root.filename
    else:
        mask_file.set("")
        mask_short.set("No File Selected")
        pick_mask_file_button.config(text="Select Mask File")





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

projections = ["cea", "mbtfpq", "aeqd", "sinu", "poly", "moerc", "gnom", "moll", "lcc", "tmerc", "nplaea", "gall",
               "npaeqd", "mill", "merc", "stere", "eqdc", "rotpole", "cyl", "npstere", "spstere", "hammer", "geos",
               "nsper", "eck4", "aea", "kav7", "spaeqd", "ortho", "class", "vandg", "laea", "splaea", "robin"]



'''     Frames, Text Variables, and Widgets for selection of the timeseries.h5 file to plot data from.     '''
pick_h5_file_frame = Frame(root)

h5_file = StringVar()
h5_file_short = StringVar()
h5_file_short.set("No File Selected")

pick_h5_file_button = Button(pick_h5_file_frame, text='Select Timeseries File', anchor='w', width=15, command=lambda: pick_file())
selected_ts_file_label = Label(pick_h5_file_frame, textvariable=h5_file_short)


'''     Frames, Text Variables, and Widgets for selection of the mask.h5 file to add a mask to the ata.     '''
pick_mask_file_frame = Frame(root)

mask_file = StringVar()
mask_short = StringVar()
mask_short.set("No File Selected")

pick_mask_file_button = Button(pick_mask_file_frame, text='Select Mask File', anchor='w', width=15, command=lambda: pick_mask())
selected_mask_file_label = Label(pick_mask_file_frame, textvariable=mask_short)


'''     Frames, Text Variables, and Widgets for setting y-lim      '''
y_lim_frame = Frame(root)
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
unit_cmap_projection_frame = Frame(root)

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


flip_frame = Frame(root)

lr_flip = IntVar()
lr_flip_checkbutton = Checkbutton(flip_frame, text="Flip LR", variable=lr_flip)

ud_flip = IntVar()
ud_flip_checkbutton = Checkbutton(flip_frame, text="Flip UD", variable=ud_flip)

wrap = IntVar()
wrap_checkbutton = Checkbutton(flip_frame, text="Wrap", variable=wrap)

opposite = IntVar()
opposite_checkbutton = Checkbutton(flip_frame, text="Opposite", variable=opposite)




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

mainloop()