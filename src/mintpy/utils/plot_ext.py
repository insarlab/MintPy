"""Class wrapped around matplotlib/mintpy for polygon selection."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Sep 2019                           #
############################################################
# Recommend import:
#   from mintpy.utils import plot_ext


import numpy as np
from matplotlib import pyplot as plt
from matplotlib.path import Path
from matplotlib.widgets import PolygonSelector

from mintpy import view


################################# SelectFromCollection class begin ########################################
class SelectFromCollection:
    """Select indices from a matplotlib collection using `PolygonSelector`.

    Selected pixels within the polygon is marked as True and saved in the
    member variable self.mask, in the same size as input AxesImage object
    with all the other pixels marked as False.

    Parameters
    ----------
    ax : :class:`~matplotlib.axes.Axes`
        Axes to interact with.

    collection : :class:`matplotlib.collections.Collection` subclass
        Collection you want to select from.

    Examples
    --------
    import matplotlib.pyplot as plt
    from mintpy.utils import readfile, plot as pp

    fig, ax = plt.subplots()
    data = readfile.read('velocity.h5', datasetName='velocity')[0]
    im = ax.imshow(data)

    selector = pp.SelectFromCollection(ax, im)
    plt.show()
    selector.disconnect()

    plt.figure()
    plt.imshow(selector.mask)
    plt.show()
    """

    def __init__(self, ax, collection, alpha_other=0.3):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.prepare_coordinates()

        self.poly = PolygonSelector(ax, self.onselect)

        msg = "\nSelect points in the figure by enclosing them within a polygon.\n"
        msg += "Press the 'esc' key to start a new polygon.\n"
        msg += "Try hold to left key to move a single vertex.\n"
        msg += "After complete the selection, close the figure/window to continue.\n"
        print(msg)

    def prepare_coordinates(self):
        imgExt = self.collection.get_extent()
        self.length = int(imgExt[2] - imgExt[3])
        self.width  = int(imgExt[1] - imgExt[0])
        yy, xx = np.mgrid[:self.length, :self.width]
        self.coords = np.hstack((xx.reshape(-1, 1),
                                 yy.reshape(-1, 1)))

    def onselect(self, verts):
        self.poly_path = Path(verts)
        self.mask = self.poly_path.contains_points(self.coords).reshape(self.length,
                                                                        self.width)
        self.canvas.draw_idle()

    def disconnect(self):
        self.poly.disconnect_events()
        self.canvas.draw_idle()

## Utility script
def get_poly_mask(fname, datasetName, print_msg=True, view_cmd=''):
    """Get mask of pixels within polygon from interactive selection
    Parameters: data : 2D np.array in size of (length, width)
    Returns:    mask : 2D np.arrat in size of (length, width) in np.bool_ format
    """
    # option 1 - Advanced plot using view.prep/plot_slice()
    cmd = f'view.py {fname} '
    if datasetName:
        cmd += f' {datasetName} '
    if view_cmd:
        cmd += view_cmd

    d_v, atr ,inps = view.prep_slice(cmd)
    ax = plt.subplots(figsize=inps.fig_size)[1]
    inps.fig_coord = 'yx'   # selector works for y/x coord plot only
    ax, inps, im = view.plot_slice(ax, d_v, atr, inps)[0:3]

    ## Option 2 - Simple plot with matplotlib
    #from mintpy.utils import readfile
    #data = readfile.read(fname, datasetName)[0]
    #vlim = np.nanmax(np.abs(data))
    #vmin, vmax = -vlim, vlim
    ## for dataset with non-negative values such as elevation
    #if np.nanmin(data) > 0:
    #    vmin = np.nanmin(data)
    ## plot
    #fig, ax = plt.subplots()
    #im = ax.imshow(data, cmap='jet', vmin=vmin, vmax=vmax)
    #fig.colorbar(im)

    selector = SelectFromCollection(ax, im)
    plt.show()
    selector.disconnect()

    if hasattr(selector, 'mask'):
        mask = selector.mask
        if print_msg:
            print(f'selected polygon: {selector.poly_path}')
    else:
        mask = None
        if print_msg:
            print('no polygon selected.\n')
    return mask
################################## SelectFromCollection class end #########################################
