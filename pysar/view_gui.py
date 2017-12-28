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

















pick_h5_file_frame.pack(anchor='w', fill=X, pady=(10, 5), padx=10)
pick_h5_file_button.pack(side=LEFT, anchor='w', padx=(0, 20))
selected_ts_file_label.pack(side=LEFT, fill=X)

pick_mask_file_frame.pack(anchor='w', fill=X)
pick_mask_file_button.pack(side=LEFT, anchor='w', pady=5, padx=(10, 20))
selected_mask_file_label.pack(side=LEFT, fill=X)

mainloop()