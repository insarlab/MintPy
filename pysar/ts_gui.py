from Tkinter import *
import matplotlib
matplotlib.use('TkAgg')

root = Tk()
root.minsize(width=350, height=550)
root.maxsize(width=350, height=550)
root.resizable(width=False, height=False)

'''     Frames, Text Variables, and Widgets for selection of the timeseries.h5 file to plot data from.     '''
pick_timeseries_file_frame = Frame(root)


'''     Frames, Text Variables, and Widgets for selection of the mask.h5 file to add a mask to the ata.     '''
pick_mask_file_frame = Frame(root)


'''     Frames, Text Variables, and Widgets for selection of the topography dem.h5 file to add topography to the data.     '''
pick_dem_file_frame = Frame(root)


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


pick_mask_file_frame.pack(anchor='w', fill=X)

pick_dem_file_frame.pack(anchor='w', fill=X)




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