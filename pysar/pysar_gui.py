#! /usr/bin/env python3

from tkinter import *
import view_gui as vg
import ts_gui as tg


def show_view_gui():
    print("VIEW GUI")
    root.destroy()
    vg.main()



def show_ts_gui():
    print("TS GUI")
    root.destroy()
    tg.main()


root = Tk()
root.minsize(width=200, height=350)
root.maxsize(width=200, height=350)
root.resizable(width=False, height=False)

view_gui_button = Button(root, text="view.py GUI", width=12, anchor='c', command=show_view_gui)
ts_gui_button = Button(root, text="tsviewer.py GUI", width=12, anchor='c', command=show_ts_gui)

view_gui_button.pack(pady=(130, 20))
ts_gui_button.pack()

mainloop()