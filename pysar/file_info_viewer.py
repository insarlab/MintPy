#! /usr/bin/env python2

from Tkinter import *
import ttk
import h5py

index = 1

file = "/Users/joshua/Desktop/geology_gui/test_data/timeseries.h5"


def add_to_tree(name, obj):
    global index
    group = tree.insert("", index, text=name)
    subindex = 0
    for key in sorted(obj.attrs.keys()):
        tree.insert(group, subindex, text=key, values=obj.attrs.get(key))
        subindex += 1
    index += 1

root = Tk()

tree = ttk.Treeview(root)
tree.configure(height=850)

'''tree["columns"]=("one","two")
tree.column("one", width=100 )
tree.column("two", width=100)
tree.heading("one", text="coulmn A")
tree.heading("two", text="column B")

tree.insert("" , 0,    text="Line 1", values=("1A","1b"))

id2 = tree.insert("", 1, "dir2", text="Dir 2")
tree.insert(id2, "end", "dir 2", text="sub dir 2", values=("2A","2B"))

##alternatively:
tree.insert("", 3, "dir3", text="Dir 3")
tree.insert("dir3", 3, text=" sub dir 3",values=("3A"," 3B"))

tree.pack()'''

tree.column("#0", width=285)
tree["columns"] = ("Value")
tree.column("Value", width=500)
tree.heading("Value", text="Value")
tree.insert("", 0, text="File", values=(file))

h5file = h5py.File(file, 'r')
h5file.visititems(add_to_tree)
h5file.close()

#tree.insert("", 0, text="FILE NAME", values=("Some File Name"))
#tree.insert("", 1, text="FILE NAME 2", values=("Some File Name 2"))

tree.pack()

root.mainloop()
