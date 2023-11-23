"""Simulate tropospheric turbulence."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2019                               #
############################################################
# Recommend usage:
#   from mintpy.simulation import fractal
# This module is based on the matlab scripts written by
#   Ramon Hanssen, May 2000, available at:
#   http://doris.tudelft.nl/software/insarfractal.tar.gz
# Reference:
#   Hanssen, R. F. (2001), Radar interferometry: data interpretation and
#   error analysis, Kluwer Academic Pub, Dordrecht, Netherlands. Chap. 4.7.


import os

import matplotlib.pyplot as plt
import numpy as np

try:
    import pyfftw
except ImportError:
    raise ImportError('Cannot import pyfftw!')


# speedup pyfftw
NUM_THREADS = min(os.cpu_count(), 4)
print(f'using {NUM_THREADS} threads for pyfftw computation.')
pyfftw.config.NUM_THREADS = NUM_THREADS


def fractal_surface_atmos(shape=(128, 128), resolution=60., p0=1., freq0=1e-3,
                          regime=(0.001, 0.999, 1.00), beta=(5./3., 8./3., 2./3.)):
    """Simulate an isotropic 2D fractal surface with a power law behavior, which corresponds with the
    [-5/3, -8/3, -2/3] power law.

    E.g. equation (4.7.28) from Hanssen (2001):
    P_phi(f) =  P2(f/f0) ^ -5/3    for 1.5  <= f0/f <= 50   km       regime[1]-regime[0]
                P1(f/f0) ^ -8/3    for 0.25 <= f0/f <= 1.5  km       regime[0]
                P3(f/f0) ^ -2/3    for 0.02 <= f0/f <= 0.25 km       regime[2]-regime[1]

    regime=[0.001, 0.999, 1.0] for larger scale turbulence
    regime=[0.980, 0.990, 1.0] for middle scale turbulence
    regime=[0.010, 0.020, 1.0] for small  scale turbulence

    This is based on the fracsurfatmo.m written by Ramon Hanssen, 2000.

    Parameters: shape      : tuple of 2 int, number of rows and columns
                resolution : float, spatial resolution in meter
                p0         : float, multiplier of power spectral density in m^2.
                regime     : tuple of 3 float, cumulative percentage of spectrum covered by a specific beta
                             e.g.: (0.60, 0.90, 1.00), (0.95, 0.99, 1.00)
                beta       : tuple of 3 float, power law exponents for a 1D profile of the data
                display    : bool, display simulation result or not
    Returns:    fsurf      : 2D np.array in size of (length, width) in m.
    Example:    data, atr = readfile.read_attribute('timeseriesResidual_ramp.h5', datasetName='20171115')
                step = abs(ut.range_ground_resolution(atr))
                p0 = get_power_spectral_density(data, resolution=step)[0]
                tropo = fractal_surface_atmos(shape=data.shape, resolution=step, p0=p0)
    """

    beta = np.array(beta, np.float32)
    regime = np.array(regime, np.float32)
    length, width = shape

    # simulate a uniform random signal
    h = np.random.rand(length, width)
    H = pyfftw.interfaces.numpy_fft.fft2(h)
    H = pyfftw.interfaces.numpy_fft.fftshift(H)

    # scale the spectrum with the power law
    yy, xx = np.mgrid[0:length-1:length*1j,
                      0:width-1:width*1j].astype(np.float32)
    yy -= np.rint(length/2)
    xx -= np.rint(width/2)
    xx *= resolution
    yy *= resolution
    k = np.sqrt(np.square(xx) + np.square(yy))    #pixel-wise distance in m

    """
    beta+1 is used as beta, since, the power exponent
    is defined for a 1D slice of the 2D spectrum:
    austin94: "Adler, 1981, shows that the surface profile
      created by the intersection of a plane and a
      2-D fractal surface is itself fractal with
      a fractal dimension equal to that of the 2D
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

    kmax = np.max(k)
    k1 = max(regime[0] * kmax, 4 * resolution)
    k2 = regime[1] * kmax

    regime1 = k <= k1
    regime2 = np.multiply(k >= k1, k <= k2)
    regime3 = k >= k2

    fraction1 = np.power(k[regime1], beta[0])
    fraction2 = np.power(k[regime2], beta[1])
    fraction3 = np.power(k[regime3], beta[2])

    fraction = np.zeros(k.shape, np.float32)
    fraction[regime1] = fraction1
    fraction[regime2] = fraction2 / np.min(fraction2) * np.max(fraction[regime1])
    fraction[regime3] = fraction3 / np.min(fraction3) * np.max(fraction[regime2])

    # prevent dividing by zero
    fraction[fraction == 0.] = 1.

    # get the fractal spectrum and transform to spatial domain
    Hfrac = np.divide(H, fraction)
    fsurf = pyfftw.interfaces.numpy_fft.ifft2(Hfrac)
    fsurf = np.abs(fsurf).astype(np.float32)
    fsurf -= np.mean(fsurf)

    # calculate the power spectral density of 1st realization
    p1 = get_power_spectral_density(fsurf, resolution=resolution, freq0=freq0)[0]

    # scale the spectrum to match the input power spectral density.
    Hfrac *= np.sqrt(p0/p1)
    fsurf = pyfftw.interfaces.numpy_fft.ifft2(Hfrac)
    fsurf = np.abs(fsurf).astype(np.float32)
    fsurf -= np.mean(fsurf)
    return fsurf


def get_power_spectral_density(data, resolution=60., freq0=1e-3, display=False, outfig=None):
    """Get the radially averaged 1D spectrum (power density) of input 2D matrix
    Check Table 4.5 in Hanssen, 2001 (Page 143) for explanation of outputs.

    Python translation of checkfr.m (Ramon Hanssen, 2000)

    Parameters: data       : 2D np.array (free from NaN value), displacement in m.
                resolution : float, spatial resolution of input data in meters
                freq0      : float, reference spatial frequency in cycle / m.
                display    : bool, display input data and its calculated 1D power spectrum
    Returns:    p0   : float, power spectral density at reference frequency in m^2
                beta : float, slope of power profile in loglog scale
                freq : 1D np.array for frequency in cycle/m
                psd1d: 1D np.array for the power spectral density in m^2
    """

    if display:
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=[10, 3])
        im = axs[0].imshow(data*100, cmap='jet')
        cbar = plt.colorbar(im, ax=axs[0])
        cbar.set_label('cm')

    # use the square part of the matrix for spectrum calculation
    data = crop_data_max_square_p2(data)
    N = data.shape[0]

    # calculate the normalized power spectrum (spectral density)
    fdata2d = pyfftw.interfaces.numpy_fft.fft2(data)
    fdata2d = pyfftw.interfaces.numpy_fft.fftshift(fdata2d)
    psd2d = np.abs(np.multiply(fdata2d, np.conj(fdata2d))) / (N**2)

    # The frequency coordinate in cycle / m
    freq = np.fft.fftfreq(N, d=resolution)
    freq = np.roll(freq, shift=np.ceil((N-1)/2).astype(int)) #shift 0 to the center
    ky, kx = np.meshgrid(freq, freq)

    # calculate the radially average spectrum
    freq, psd1d = radial_average_spectrum(kx, ky, psd2d)

    # calculate slopes from spectrum
    p0, beta = power_slope(freq, psd1d, freq0=freq0)
    # D2 = (7. - beta + 1.) / 2.

    if display:
        ax = axs[1]
        # plot
        ax.loglog(freq*1e3, psd1d*1e4)

        # reference frequency
        ax.axvline(freq0*1e3, linestyle='--', color='k')

        # axis format
        ax.set_xlabel('Wavenumber [cycle/km]')
        ax.set_ylabel('Power '+r'$[cm^2]$')
        msg = r'$f_0=$'+f'{freq0*1e3:.3f} '+r'$km^{-1}$'
        msg += '\n'+fr'$p_0={p0*1e4:.1f}\/cm^2$'
        msg += '\n'+fr'$\beta={beta:.2f}$'
        ax.text(0.6, 0.55, msg, transform=ax.transAxes, fontsize=12)
        fig.tight_layout()

        if outfig:
            fig.savefig(outfig, bbox_inches='tight', transparent=True, dpi=300)
            print('save figure to', outfig)
        plt.show()
    return p0, beta, freq, psd1d


def crop_data_max_square_p2(data):
    """Grab the max portion of the input 2D matrix that it's:
        1. square in shape
        2. dimension as a power of 2
    """
    # get max square size in a power of 2
    N = min(data.shape)
    N = np.power(2, int(np.log2(N)))

    # find corner with least number of zero values
    flag = data != 0
    if data.shape[0] > data.shape[1]:
        num_top = np.sum(flag[:N, :N])
        num_bottom = np.sum(flag[-N:, :N])
        if num_top > num_bottom:
            data = data[:N, :N]
        else:
            data = data[-N:, :N]
    else:
        num_left = np.sum(flag[:N, :N])
        num_right = np.sum(flag[:N, -N:])
        if num_left > num_right:
            data = data[:N, :N]
        else:
            data = data[:N, -N:]
    return data


def power_slope(freq, psd, freq0=1e-3):
    """ Derive the slope beta and p0 of an exponential function in loglog scale
    p = p0 * (freq/freq0)^(-beta)

    Python translation of pslope.m (Ramon Hanssen, 2000)

    Parameters: freq  : 1D / 2D np.array in cycle / m.
                psd   : 1D / 2D np.array for the power spectral density
                freq0 : reference frequency in cycle / m.
    Returns:    p0    : float, power spectral density at reference frequency
                        in the same unit as the input psd.
                beta  : float, slope of power profile in loglog scale
    """
    freq = freq.flatten()
    psd = psd.flatten()

    # check if there is zero frequency. If yes, remove it.
    if not np.all(freq != 0.):
        idx = freq != 0.
        freq = freq[idx]
        psd = psd[idx]

    # convert to log-log scale
    logf = np.log10(freq)
    logp = np.log10(psd)

    # fit a linear line
    beta = -1 * np.polyfit(logf, logp, deg=1)[0]

    # interpolate psd at reference frequency
    if freq0 < freq[0] or freq0 > freq[-1]:
        raise ValueError(f'input frequency of interest {freq0} is out of range ({freq[0]}, {freq[-1]})')
    position = np.interp(np.log10(freq0), logf, range(len(logf)))
    logp0 = np.interp(position, range(len(logp)), logp)
    p0 = np.power(10, logp0)
    return p0, beta


def radial_average_spectrum(kx, ky, pds2d):
    """Calculate the radially averaged power spectrum
    Parameters: kx    - 2D np.ndarray in size of (N,N) for frequency in x direction
                kx    - 2D np.ndarray in size of (N,N) for frequency in y direction
                psd2d - 2D np.ndarray in size of (N,N) for 2D power spectral density
    Returns:    freq  - 1D np.ndarray in size of int(N/2 - 1) for frequency in radial direction
                psd1d - 1D np.ndarray in size of int(N/2 - 1) for power spectral density
    """
    N = kx.shape[0]
    df = np.unique(kx)[np.unique(kx) > 0][0]
    k = np.sqrt(np.square(kx) + np.square(ky))

    # rotationally averaged 1D spectrum from 2D spectrum
    psd1d = np.zeros(int(N/2-1))
    for i in range(len(psd1d)):
        ind = np.multiply(k >= (i+0.5)*df, k <= (i+1.5)*df)
        psd1d[i] = np.mean(pds2d[ind])

    # Only consider one half of spectrum (due to symmetry)
    freq = np.arange(1, N/2) * df
    return freq, psd1d


def recon_power_spectral_density(N, step, p0, beta, f0=1e-4):
    """Reconstruct 1D power spectral density from input p0 and beta

    Parameters: N    - int, min size of the 2D matrix
                step - float, spatial resolution of the 2D matrix in meters
                p0   - float / 1D np.ndarray, power spectral density in m^2 at frequency of f0
                beta - float / 1D np.ndarray, power spectra slope in loglog scale
                f0   - float, reference spatial frequency in cycle / m
    Returns:    f    - 1D np.ndarray, spatial frequency sequence in cycle / m in size of (num_freq)
                p    - 1D/2D np.ndarray, power spectral density in size of (num_psd, num_freq)
    Examples:   atr = readfile.read_attribute('ERA5')
                N = min(int(atr['LENGTH']), int(atr['WIDTH']))
                step = abs(ut.range_ground_resolution(atr))
                freq, psd = recon_power_spectral_density(N, step, p0, beta)
    """
    # frequency for x-axis after FFT
    f = np.fft.fftfreq(N, d=step)
    df = np.unique(f)[np.unique(f) > 0][0]
    f = np.arange(1, N/2) * df
    f = np.array(f, dtype=np.float32).reshape(1,-1)

    # p0, beta --> psd
    p0 = np.array(p0, dtype=np.float32).reshape(-1,1)
    beta = np.array(beta, dtype=np.float32).reshape(-1,1)

    logf = np.log10(f)
    logf0 = np.log10(f0)
    logp = -beta * (logf - logf0) + np.log10(p0)
    p = np.power(10, logp)
    f = f.flatten()
    return f, p
