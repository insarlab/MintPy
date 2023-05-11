"""Forward deformation models."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2018                               #
############################################################
# Recommend usage:
#   from mintpy.simulation import defo_model as defo
#   from mintpy.simulation import simulation as sim


import numpy as np
from matplotlib import pyplot as plt

from mintpy.utils import utils0 as ut0


def mogi(geometry, xloc, nu=0.25):
    """Computes surface displacements, strains, and tilts due to a Mogi source.
    Parameters: geometry : tuple of 4 float, Mogi source geometry: East, North, Depth, Volomn change in SI unit.
                xloc : 2D np.array in size of (2, num_pixel), surface coordinates in east/x and north/y
                nu : float, Poisson's ratio
    Returns:    u : 2D np.array of displacement in size of (3, num_pixel) for Ux, Uy, Uz
                e : 2D np.array of strains      in size of (3, num_pixel) for Exx, Exy, Eyy
                t : 2D np.array of tilts        in size of (2, num_pixel) for dUz/dx, dUz/dy
    Notes: depth denotes an unsigned length and should therefore always be given positive. Keep your units consistent!
    This is a python translation from mogi.m originally written by Peter Cervelli, May 1998.
    """
    xloc = np.array(xloc, np.float32).reshape(2, -1)

    # compute displacements
    num_pixel = xloc.shape[1]
    E = geometry[0] - xloc[0, :]
    N = geometry[1] - xloc[1, :]
    E2 = np.square(E)
    N2 = np.square(N)
    d2 = np.square(geometry[2])
    R = np.sqrt(d2 + E2 + N2)
    C = (nu - 1.) * geometry[3] / np.pi

    R3 = C * np.power(R, -3)
    displacement = np.zeros((3, num_pixel), np.float32)
    displacement[0, :] = np.multiply(E, R3)
    displacement[1, :] = np.multiply(N, R3)
    displacement[2, :] = np.multiply(-1 * geometry[2], R3)

    # compute strains (if necessary)
    R5 = C * np.power(R, -5)
    strain = np.zeros((3, num_pixel), np.float32)
    strain[0, :] = np.multiply(R5, (2. * E2 - N2 - d2))
    strain[1, :] = 3. * np.multiply(R5, np.multiply(E, N))
    strain[2, :] = np.multiply(R5, (2. * N2 - E2 - d2))

    # compute tilts
    tilt = np.zeros((2, num_pixel), np.float32)
    tilt[0, :] = -3. * np.multiply(R5, E) * geometry[2]
    tilt[1, :] = -3. * np.multiply(R5, N) * geometry[2]

    return displacement, strain, tilt


def mogi_los(shape, source_geom, resolution=60., scale=1., display=True):
    """Simulate 2D deformation caused by the overpress of a Mogi source underneath

    Parameters: shape: 2-tuple of int in (length, width) or 2D np.ndarray in size of (length, width) in np.bool_
                source_geom : 4-tuple of float, Mogi source geometry: East, North, Depth, Volomn change in SI unit.
    Returns:    dis_los: 2D np.ndarray in size of (length, width), deformation in LOS direction in meter
    """
    if isinstance(shape, np.ndarray):
        mask = np.multiply(np.array(shape != 0), ~np.isnan(shape))
        shape = mask.shape
    else:
        mask = np.ones(shape, np.bool_)

    length, width = shape
    yy, xx = np.mgrid[0:length:length*1j, 0:width:width*1j]
    yy *= resolution
    xx *= resolution
    xloc = np.vstack((xx.reshape(1, -1), yy.reshape(1, -1)))

    dis_map = mogi(source_geom, xloc)[0]
    dis_e = dis_map[0, :].reshape(length, width)
    dis_n = dis_map[1, :].reshape(length, width)
    dis_u = dis_map[2, :].reshape(length, width)

    dis_los = ut0.enu2los(dis_e, dis_n, dis_u,
                          inc_angle=34.,
                          head_angle=-168.)
    dis_los[mask == 0.] = np.nan
    dis_los *= scale

    if display:
        fig, ax = plt.subplots(1, 4, figsize=[10, 3], sharey=True)
        dmin = np.nanmin(dis_los)
        dmax = np.nanmax(dis_los)
        for i, fig_title in enumerate(['east','north','vertical']):
            ax[i].imshow(dis_map[i, :].reshape(length, width), vmin=dmin, vmax=dmax)
            ax[i].set_title(fig_title)
        im = ax[3].imshow(dis_los, vmin=dmin, vmax=dmax)
        ax[3].set_title('los - SenD')
        fig.subplots_adjust(right=0.90)
        cax = fig.add_axes([0.92, 0.25, 0.01, 0.5])
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label('Displacement [m]')
        plt.show()

    return dis_los
