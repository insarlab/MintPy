############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2018, Zhang Yunjun                          #
# Author:  Zhang Yunjun, 2018                              #
############################################################


import numpy as np


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

