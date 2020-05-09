#! /usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Joshua Zahner, April 2018                        #
############################################################


from tkinter import *
from tkinter import ttk
from tkinter import filedialog
import matplotlib
matplotlib.use('TkAgg')
import h5py
from mintpy.utils import readfile
from mintpy import view

index = 2
root = None
file = None
pick_h5_file_button = None
tree = None


def add_to_tree(name, obj):
    global index
    group = tree.insert("", index, text=name)
    subindex = 0
    for key in sorted(obj.attrs.keys()):
        tree.insert(group, subindex, text=key, values=obj.attrs.get(key))
        subindex += 1
    index += 1


def add_attrs_to_tree(attrs):
    group = tree.insert("", 1, text="Attributes", values=())
    subindex = 0
    for key in sorted(attrs.keys()):
        tree.insert(group, subindex, text=key, values=attrs.get(key))
        subindex += 1


def plot_data():

    atr = readfile.read_attribute(file.get())
    file_type = atr['FILE_TYPE']

    datasets = readfile.get_dataset_list(file.get(), file_type)

    item = tree.focus()
    the_item = tree.item(item)
    epoch_num = the_item['text']
    if epoch_num in datasets:
        view.main([file.get(), epoch_num])


def pick_file():
    global tree

    if file.get() == "":
        filename = filedialog.askopenfilename(initialdir="/", title="Select file",
                                              filetypes=(("HDF5 files", "*.h5"),
                                                         ("all files", "*.*"),
                                                         ("more files", "*.he5")))
        if tree is not None:
            tree.destroy()

        file.set(filename)
        setup_tree(file.get())

        pick_h5_file_button.config(text="Clear")

    else:
        file.set("")
        pick_h5_file_button.config(text="Select .h5 File")


def main():
    global tree, file, pick_h5_file_button, root

    root = Tk()
    root.minsize(width=365, height=750)

    reset_settings_file_frame = Frame(root)
    reset_settings_file_frame.pack(side=TOP, fill="none", expand=True)

    show_button = Button(reset_settings_file_frame, text="Show", width=10, command=plot_data)
    show_button.pack(side=LEFT, pady=(10, 5))

    file = StringVar()
    h5_file_short = StringVar()
    h5_file_short.set("No File Selected")

    pick_h5_file_button = Button(reset_settings_file_frame, text='Select .h5 File', width=15, command=pick_file)

    pick_h5_file_button.pack(side=LEFT, pady=(10, 5))

    root.mainloop()


def setup_tree(file):
    global tree, root

    tree = ttk.Treeview(root)
    tree.configure(height=40)
    tree.column("#0", width=150)
    tree["columns"] = "Value"
    tree.column("Value", width=350)
    tree.heading("Value", text="Value")
    tree.insert("", 0, text="File", values=file)
    h5file = h5py.File(file, 'r')
    add_attrs_to_tree(h5file.attrs)
    h5file.visititems(add_to_tree)
    h5file.close()
    tree.pack(side=TOP)


if __name__ == '__main__':
    main()
