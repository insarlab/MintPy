from Tkinter import *
import matplotlib
matplotlib.use('TkAgg')
import tkFileDialog as filedialog

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

root = Tk()
root.minsize(width=350, height=550)
root.maxsize(width=350, height=550)
root.resizable(width=False, height=False)

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


'''     Frames, Text Variables, and Widgets for start long coordinate entry     '''
start_lon_frame = Frame(start_lat_lon_frame)


'''     Frames, Text Variables, and Widgets for reference lat/long coordinates     '''
ref_lat_lon_frame = Frame(lat_lon_frame, width=10)

'''     Frames, Text Variables, and Widgets for reference lat coordinate entry     '''
ref_lat_frame = Frame(ref_lat_lon_frame)

'''     Frames, Text Variables, and Widgets for reference long coordinate entry     '''
ref_lon_frame = Frame(ref_lat_lon_frame)


'''     Frames, Text Variables, and Widgets for setting y-lim      '''
y_lim_frame = Frame(root)

'''     Frames, Text Variables, and Widgets for setting upper bound on the y-lim      '''
y_lim_upper_frame = Frame(y_lim_frame)

'''     Frames, Text Variables, and Widgets for setting lower bound on the y-lim      '''
y_lim_lower_frame = Frame(y_lim_frame)

'''     Frames, Text Variables, and Widgets for setting extrenous widgets      '''
unit_cmap_frame = Frame(root)






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

start_lon_frame.pack(side=LEFT)

ref_lat_lon_frame.pack(side=TOP, fill=X)

ref_lat_frame.pack(side=LEFT)


ref_lon_frame.pack(side=LEFT)

y_lim_frame.pack(anchor='w', fill=X, pady=10, padx=10)

y_lim_upper_frame.pack(side=TOP, fill=X, pady=(0, 10))

y_lim_lower_frame.pack(side=TOP, fill=X)

unit_cmap_frame.pack(anchor='w', fill=X, pady=10, padx=10)

mainloop()