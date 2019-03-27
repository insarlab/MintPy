#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2018, Zhang Yunjun                          #
# Author:  Zhang Yunjun                                    #
############################################################
# This is a python translation of the matlab scripts originally
# written by Ramon Hanssen, May 2000, available through the 
# following website:
#     http://doris.tudelft.nl/software/insarfractal.tar.gz
# Reference:
#   Hanssen, R. F. (2001), Radar interferometry: data interpretation 
# and error analysis, Kluwer Academic Pub, Dordrecht, Netherlands. Chap. 4.7.


import numpy as np
import matplotlib.pyplot as plt


def fractal_surface_atmos(shape=(128, 128), resolution=60., p0=1., regime=(95., 99., 100.),
                          beta=(5./3., 8./3., 2./3.), display=False):
    """Simulate an isotropic 2D fractal surface with a power law behavior.

    This is a python translation of fracsurfatmo.m (Ramon Hanssen, 2000).
    The difference from Hanssen (2000) is:
        1) support non-square shape output
        2) take spatial resolution into consideration
        3) divide p0 with a ratio of 0.073483 to make p0 == C0 from check_power_spectrum_1d

    Parameters: shape   : tuple of 2 int, number of rows and columns
                resolution : float, spatial resolution in meter
                p0      : float, multiplier of spectrum amplitude
                beta    : tuple of 3 float, power law exponents
                    for a 1D profile of the data
                regime  : tuple of 3 float, cumulative percentage
                    of spectrum covered by a specific beta
                display : bool, display simulation result or not
    Returns:    fsurf   : 2D np.array in size of (LENGTH, WIDTH)
    Example:    data, atr = readfile.read_attribute('timeseries.h5',
                                                    datasetName='20171115')
                length, width = int(atr['LENGTH']), int(atr['WIDTH'])
                step = abs(ut.range_ground_resolution(atr))
                p0 = check_power_spectrum_1d(data, resolution=step)
                sim_trop_turb = frac_surf_atmos(shape=(length, width),
                                                resolution=step,
                                                p0=p0)
    """

    beta = np.array(beta, np.float32)
    regime = np.array(regime, np.float32)
    length, width = shape

    # simulate a uniform random signal
    h = np.random.rand(length, width)
    H = np.fft.fftshift(np.fft.fft2(h))

    # scale the spectrum with the power law
    yy, xx = np.mgrid[0:length-1:length*1j,
                      0:width-1:width*1j]
    yy -= int(length / 2. + 0.5)
    xx -= int(width / 2. + 0.5)
    xx *= resolution / 1000.
    yy *= resolution / 1000.
    k = np.sqrt(np.square(xx) + np.square(yy))

    """
    beta+1 is used as beta, since, the power exponent
    is defined for a 1D slice of the 2D spectrum:
    austin94: "Adler, 1981, shows that the surface profile 
      created by the intersection of a plane and a
      2-D fractal surface is itself fractal with 
      a fractal dimension  equal to that of the 2D 
      surface decreased by one.
    """
    beta += 1.

    """
    The power beta/2 is used because the power spectral
    density is proportional to the amplitude squared 
    Here we work with the amplitude, instead of the power
    so we should take sqrt( k.^beta) = k.^(beta/2)  RH
    """
    beta /= 2.
    
    mk = np.max(k)
    k0 = 0
    k1 = (regime[0] / 100.) * mk
    k2 = (regime[1] / 100.) * mk
    k3 = np.max(k)
    
    regime1 = np.multiply(k >  k0, k <= k1)
    regime2 = np.multiply(k >= k1, k <= k2)
    regime3 = np.multiply(k >= k2, k <= k3)
    
    fraction1 = np.power(k[regime1], beta[0])
    fraction2 = np.power(k[regime2], beta[1])
    fraction3 = np.power(k[regime3], beta[2])
    
    fraction = np.zeros(k.shape, np.float32)
    fraction[regime1] = fraction1
    fraction[regime2] = fraction2 / np.min(fraction2) * np.max(fraction[regime1])
    fraction[regime3] = fraction3 / np.min(fraction3) * np.max(fraction[regime2])

    # prevent dividing by zero
    fraction[fraction == 0.] = 1.

    # test
    debug = False
    if debug:
        fig = plt.figure('100')
        plt.loglog(k[regime1], 1./fraction[regime1], '.')
        plt.loglog(k[regime2], 1./fraction[regime2], '.')
        plt.loglog(k[regime3], 1./fraction[regime3], '.')
        #plt.xlabel('Wavenumber')
        #plt.ylabel('Power')
        plt.show()

    Hnew = p0 / 0.073483 * np.divide(H, fraction)
    # create spectral surface by ifft
    fsurf = np.abs(np.fft.ifft2(Hnew))
    # remove mean to get zero-mean data
    fsurf -= np.mean(fsurf)
    fsurf = np.array(fsurf, np.float32)

    if display:
        plt.figure()
        plt.imshow(fsurf)
        plt.colorbar()
        plt.show()
    return fsurf


def power_slope(freq, power):
    """ Derive the slope beta and C0 of an exponential function in loglog scale
    S(k) = np.power(C0, 2) * np.power(k, -beta)
    power = np.power(C0, 2) * np.power(freq, -beta)

    Python translation of pslope.m (Ramon Hanssen, 2000)
    
    Parameters: freq  : 1D / 2D np.array
                power : 1D / 2D np.array
    Returns:    C0    : float, spectral power density at freq == 0.
                beta  : float, slope of power profile
    """
    freq = freq.flatten()
    power = power.flatten()

    # check if there is zero frequency. If yes, remove it.
    if not np.all(freq != 0.):
        idx = freq != 0.
        freq = freq[idx]
        power = power[idx]

    logf = np.log10(freq)
    logp = np.log10(power)

    beta = -1 * np.polyfit(logf, logp, deg=1)[0]

    position = np.interp(0, logf, range(len(logf)))
    logC0 = np.interp(position, range(len(logp)), logp)
    C0 = np.sqrt(np.power(10, logC0))
    return C0, beta


def check_power_spectrum_1d(data, resolution=60., display=False):
    """Check the rotationally averaged 1D spectrum of input 2D matrix
    Check Table 4.5 in Hanssen, 2001 (Page 143) for explaination of outputs.

    Python translation of checkfr.m (Ramon Hanssen, 2000)

    Parameters: data       : 2D np.array (free from NaN value)
                resolution : float, spatial resolution of input data in meters
                display    : bool, display input data and its calculated 1D power spectrum
    Returns:    C0   : float, spectral power density at freq == 0.
                beta : float, slope of power profile
                D2   : fractal dimension
    """

    if display:
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[10, 3])
        im = ax[0].imshow(data)
        plt.colorbar(im, ax=ax[0])

    # use the square part of the matrix for spectrum calculation
    N = min(data.shape)
    data = data[:N, :N]

    # The frequency coordinate
    dx = resolution / 1000.
    Nyquist = 1. / (2. * dx)
    df = 1. / (N * dx)
    freq  = np.arange(-Nyquist, Nyquist, df)
    FY, FX = np.meshgrid(freq, freq)
    k = np.sqrt(np.square(FX) + np.square(FY))

    # rotationally averaged 1D spectrum from 2D spectrum
    fdata2d = np.fft.fftshift(np.fft.fft2(data))
    pow2d = np.abs((1. / N**2) * np.multiply(fdata2d, np.conj(fdata2d)))
    pow1d = np.zeros(int(N/2-1))
    for i in range(len(pow1d)):
        sel = np.multiply(k >= (i+1)*df, k <= (i+2)*df)
        pow1d[i] = np.mean(pow2d[sel])

    freq = np.arange(1, N/2) * df

    # calculate slopes from spectrum
    C0, beta = power_slope(freq, pow1d)
    D2 = (7. - beta + 1.) / 2.

    if display:
        print('C0 = {:.6f}, beta = {:.3f}, Fractal dim = {:.1f}'.format(C0, beta, D2))
        ax[1].loglog(freq, pow1d)
        plt.show()
    return C0, beta, D2


