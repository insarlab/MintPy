#! /usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Dec 2019                           #
############################################################
# Recommended usage:
#   from mintpy.simulation import decorrelation as decor
#   OR from mintpy.simulation import simulation as sim


import math
import numpy as np
import numpy.polynomial.polynomial as poly
import scipy.stats as stats
from matplotlib import pyplot as plt
from mintpy.utils import ptime


######################################## Statistics ############################################
def phase_pdf_ds(L, coherence=None, phi_num=1000, epsilon=1e-3):
    """Marginal PDF of interferometric phase for distributed scatterers (DS)
    Eq. 66 (Tough et al., 1995) and Eq. 4.2.23 (Hanssen, 2001)
    Parameters: L         - int, number of independent looks
                coherence - 1D np.array for the range of coherence, with value < 1.0 for valid operation
                phi_num   - int, number of phase sample for the numerical calculation
    Returns:    pdf       - 2D np.array, phase PDF in size of (phi_num, len(coherence))
                coherence - 1D np.array for the range of coherence
    Example:
        epsilon = 1e-4
        coh = np.linspace(0., 1-epsilon, 1000)
        pdf, coh = phase_pdf_ds(1, coherence=coh)
    """

    def gamma(x):
        """
        Gamma function equivalent to scipy.special.gamma(x)
        :param x: float
        :return: float
        """
        # This function replaces scipy.special.gamma(x).
        # It is needed due to a bug where Dask workers throw an exception in which they cannot
        # find `scipy.special.gamma(x)` even when it is imported.

        # When the output of the gamma function is undefined, scipy.special.gamma(x) returns float('inf')
        # whereas math.gamma(x) throws an exception.
        try:
            return math.gamma(x)
        except ValueError:
            return float('inf')


    if coherence is None:
        coherence = np.linspace(0., 1.-epsilon, 1000)
    coherence = np.array(coherence, np.float64).reshape(1, -1)
    phi = np.linspace(-np.pi, np.pi, phi_num, dtype=np.float64).reshape(-1, 1)

    # Phase PDF - Eq. 4.2.32 (Hanssen, 2001)
    A = np.power((1-np.square(coherence)), L) / (2*np.pi)
    A = np.tile(A, (phi_num, 1))
    B = gamma(2*L - 1) / ((gamma(L))**2 * 2**(2*(L-1)))

    beta = np.multiply(np.abs(coherence), np.cos(phi))
    C = np.divide((2*L - 1) * beta, np.power((1 - np.square(beta)), L+0.5))
    C = np.multiply(C, (np.pi/2 + np.arcsin(beta)))
    C += 1 / np.power((1 - np.square(beta)), L)

    sumD = 0
    if L > 1:
        for r in range(int(L)-1):
            D = gamma(L-0.5) / gamma(L-0.5-r)
            D *= gamma(L-1-r) / gamma(L-1)
            D *= (1 + (2*r+1)*np.square(beta)) / np.power((1 - np.square(beta)), r+2)
            sumD += D
        sumD /= (2*(L-1))

    pdf = B*C + sumD
    pdf = np.multiply(A, pdf)
    return pdf, coherence.flatten()


def phase_variance_ds(L,  coherence=None, epsilon=1e-3):
    """Interferometric phase variance for distributed scatterers (DS)
    Eq. 2.1.2 (Box et al., 2015) and Eq. 4.2.27 (Hanssen, 2001)
    Parameters: L         - int, number of independent looks
                coherence - 1D np.array for the range of coherence, with value < 1.0 for valid operation
    Returns:    var       - 1D np.array, phase variance in size of (len(coherence)) in radians^2
                coherence - 1D np.array for the range of coherence
    Example:
        epsilon = 1e-4
        coh = np.linspace(0., 1-epsilon, 1000)
        var, coh = phase_variance_ds(1, coherence=coh)
    """
    if coherence is None:
        coherence = np.linspace(0., 1.-epsilon, 1000, dtype=np.float64)
    phi_num = len(coherence)

    phi = np.linspace(-np.pi, np.pi, phi_num, dtype=np.float64).reshape(-1, 1)
    phi_step = 2*np.pi/phi_num

    pdf, coherence = phase_pdf_ds(L, coherence=coherence, phi_num=phi_num)
    var = np.sum(np.multiply(np.square(np.tile(phi, (1, len(coherence)))), pdf)*phi_step, axis=0)

    # assign negative value (when coherence is very close to 1, i.e. 0.999) to the min positive value
    flag = var <= 0
    if not np.all(flag):
        var[flag] = np.nanmin(var[~flag])
    else:
        var[flag] = np.finfo(np.float64).eps
    return var, coherence


def phase_variance_ps(L, coherence=None, epsilon=1e-3):
    """
    the Cramer-Rao bound (CRB) of phase variance
    Given by Eq. 25 (Rodriguez and Martin, 1992)and Eq 4.2.32 (Hanssen, 2001)
    Valid when coherence is close to 1.

    Parameters: L         - int, number of independent looks applied
                coherence - float / np.ndarray, spatial coherence
                epsilon   - float, a small number to avoid denominator to be zero
    Returns:    var       - float / np.ndarray, phase variance in radians^2
                coherence - float / np.ndarray, spatial coherence (for reference purpose)
    """
    if coherence is None:
        coherence = np.linspace(0.9, 1.-epsilon, 1000, dtype=np.float64)
    var = (1 - coherence**2) / (2 * int(L) * coherence**2)
    return var, coherence


########################################## Simulations #########################################
def calibrate_coherence4phase_pdf_bias(coh, L, coh_step=0.02, num_sample=1e5, print_msg=True):
    """Calculate coherence for the bias caused by the imperfect phase PDF of DS.

    By comparing the estimated coherence from the simulated decorrelation noise based on its PDF
    with the specified true coherence.

    A very large sample size is used to avoid the overestimation bias due to limited sample size.
    No other phase components are added to avoid the underestimation bias due to the underlying fringe.

    Parameters: coh        - 1/2/3D np.ndarray of float32 for spatial coherence
                L          - int, number of independent looks
                coh_step   - float, step of coherence to generate lookup table
                num_sample - int, number of samples used for coherence estimation
    Returns:    coh_cal    - 1/2/3D np.ndarray of float32 for calibrated coherence
    """

    ## 1. calculate the bias due to phase PDF numerically
    # list of coherence
    num_coh = int(1 / coh_step)
    coh_sim = np.linspace(1, coh_step, num=num_coh)
    coh_sim -= coh_step/2   # use value in the bin center

    # loop over different coherence
    coh_est = np.zeros(num_coh, dtype=np.float32)
    for i in range(num_coh):
        # sim decor noise
        sim_pha = sample_decorrelation_phase(coh_sim[i], L=L, size=num_sample)

        # form interferometric phase
        sim_int = np.array(np.exp(1j*sim_pha), dtype=np.complex64)

        # calc coherence from interferometric phase
        coh_est[i] = np.abs(np.mean(sim_int))

    ## 2. fit a polynomial curve to calibrate the bias
    if print_msg:
        print('calibrate the coherence for the phase PDF bias with a polynomial of degree 10')
    ffit = poly.Polynomial(poly.polyfit(coh_est, coh_sim, deg=10, rcond=1e-15))
    coh_cal = ffit(coh)

    return coh_cal


def coherence2decorrelation_phase(coh, L, coh_step=0.01, num_repeat=1, scale=1.0, display=False, print_msg=True):
    """Simulate decorrelation phase based on coherence array/matrix
    based on the phase PDF of DS from Tough et al. (1995).

    Parameters: coh        - 1/2/3D np.ndarray of float32 for spatial coherence
                L          - int, number of independent looks
                coh_step   - float, step of coherence to generate lookup table
                num_repeat - int, number of repeatetion
    Returns:    pha        - np.ndarray of float32 for decorrelation phase in radians
                                for num_repeat == 1, pha.shape = coh.shape
                                for num_repeat >  1, pha.shape = (coh.size, num_repeat)
    """
    shape_orig = coh.shape
    coh = coh.reshape(-1,1)
    num_coh = coh.size

    # check number of looks
    L = int(L)
    msg = 'number of independent looks L={}'.format(L)
    if L > 80:
        L = 80
        msg += ', use L=80 to avoid dividing by 0 in calculation with negligible effect'
    if print_msg:
        print(msg)

    # code for debug
    debug_mode = False
    if debug_mode is True:
        decor = sample_decorrelation_phase(0.4, L, size=int(1e4), scale=scale)
        return decor

    # initiate output matrix
    pha = np.zeros((num_coh, num_repeat), dtype=np.float32)

    # sampling strategy
    num_step = int(1 / coh_step)
    if num_coh <= num_step * 2:
        # for small size --> loop through each input coherence
        for i in range(num_coh):
            pha[i,:] = sample_decorrelation_phase(coh[i], L, size=num_repeat, scale=scale)

    else:
        # for large size --> loop through coherence lookup table and map into input coherence matrix
        prog_bar = ptime.progressBar(maxValue=num_step, print_msg=print_msg)
        for i in range(num_step):
            # find index within the coherence intervals
            coh_i = i * coh_step
            flag = np.multiply(coh >= coh_i, coh < coh_i + coh_step).flatten()
            num_coh_i = np.sum(flag)
            if num_coh_i > 0:
                pha_i = sample_decorrelation_phase(coh_i, L, size=num_coh_i*num_repeat, scale=scale)
                pha[flag,:] = pha_i.reshape(-1, num_repeat)
            prog_bar.update(i+1, suffix='{:.3f}'.format(coh_i))
        prog_bar.close()

    if num_repeat == 1:
        pha = pha.reshape(shape_orig)

    # plot
    if display and len(shape_orig) == 2:
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=[6,3])
        axs[0].imshow(coh.reshape(shape_orig), vmin=0, vmax=1, cmap='gray')
        axs[1].imshow(pha[:,0].reshape(shape_orig), vmin=-np.pi, vmax=np.pi, cmap='jet')
        plt.show()

    return pha


def sample_decorrelation_phase(coherence, L, size=1, phi_num=1000, display=False, scale=1.0, font_size=12):
    '''Sample decorrelation phase based on PDF determined by L and coherence value
    Parameters: coherence - float, spatial coherence
                L         - int, number of independent looks
                phi_num   - int, sample number
    Returns:    phase     - 1D np.array in size of (size,), sampled phase in radians
    Examples:
        decor_noise = sample_decorrelation_phase(0.7, L=1, size=1e4, display=True)
    '''
    size = int(size)
    phiMax = np.pi * float(scale)

    # numerical solution of phase PDF for distributed scatterers
    pdf = phase_pdf_ds(int(L), coherence, phi_num=phi_num)[0].flatten()

    # generate phase distribution
    phi = np.linspace(-phiMax, phiMax, phi_num+1, endpoint=True)
    phi_dist = stats.rv_histogram((pdf, phi))

    # sample from the distribution
    phase = phi_dist.rvs(size=size)

    if display:
        fig, ax = plt.subplots(figsize=[5,3])
        ax.hist(phase, bins=50, density=True, label='Sample\nHistogram\n(norm)')
        ax.plot(phi, phi_dist.pdf(phi), label='PDF')
        ax.plot(phi, phi_dist.cdf(phi), label='CDF')
        ax.set_xlabel('Phase', fontsize=font_size)
        ax.set_ylabel('Probability', fontsize=font_size)
        ax.set_title(r'L = %d, $\gamma$ = %.2f, sample size = %d' % (L, coherence, size), fontsize=font_size)
        ax.set_xlim([-np.pi, np.pi])
        ax.set_xticks([-np.pi, 0, np.pi])
        ax.set_xticklabels([r'-$\pi$', '0', r'$\pi$'], fontsize=font_size)
        ax.tick_params(direction='in', labelsize=font_size)
        ax.legend(fontsize=font_size)
        plt.savefig('DecorNoiseSampling.jpg', bbox_inches='tight', dpi=600)
        plt.show()
    return phase


#################################### Weight Functions #####################################
def coherence2phase_variance(coherence, L=32, scatter='DS', epsilon=1e-3, print_msg=False):
    """Convert coherence to phase variance, depending on the type of scatterers

    For DS, it is equation (66) from Tough et al. (1995).
    For PS, it is equation (25) from Rodriguez and Martin (1992).

    Parameters: coherence - 1/2/3D np.ndarray of float32 for spatial coherence
                L         - int, number of independent looks
                scatter   - str, type of scatterers, PS or DS
    Returns:    variance  - 1/2/3D np.ndarray of float32 for the phase variance in radians^2
    """
    lineStr = '    number of independent looks L={}'.format(L)
    if L > 80:
        L = 80
        lineStr += ', use L=80 to avoid dividing by 0 in calculation with negligible effect'
    if print_msg:
        print(lineStr)

    coh_num = 1000
    coh_min = 0.0 + epsilon
    coh_max = 1.0 - epsilon
    coh_lut = np.linspace(coh_min, coh_max, coh_num)
    coh_min = np.min(coh_lut)
    coh_max = np.max(coh_lut)
    coh_step = (coh_max - coh_min) / (coh_num - 1)

    coherence = np.array(coherence)
    coherence[coherence < coh_min] = coh_min
    coherence[coherence > coh_max] = coh_max
    coherence_idx = np.array((coherence - coh_min) / coh_step, np.int16)

    scatter = scatter.upper()
    if scatter == 'DS':
        var_lut = phase_variance_ds(int(L), coh_lut)[0]
        variance = var_lut[coherence_idx]

    elif scatter == 'PS':
        variance = phase_variance_ps(int(L), coherence)[0]

    else:
        raise ValueError('un-recognized scatterer type: {}'.format(scatter))

    return variance


def coherence2fisher_info_index(data, L=32, epsilon=1e-3):
    """Convert coherence to Fisher information index (Seymour & Cumming, 1994, IGARSS)"""
    print('    number of independent looks L={}'.format(L))
    # prepare input data
    if data.dtype != np.float64:
        data = np.array(data, np.float64)
    data[data > 1-epsilon] = 1-epsilon
    # calculation
    data = 2.0 * int(L) * np.square(data) / (1 - np.square(data))
    return data


def coherence2weight(coh_data, weight_func='var', L=20, epsilon=5e-2, print_msg=True):
    """Convert spaital coherence to weight given the weight function name"""
    coh_data[np.isnan(coh_data)] = epsilon
    coh_data[coh_data < epsilon] = epsilon
    coh_data = np.array(coh_data, np.float64)

    # make sure L >= 1
    L = max(L, 1)

    # Calculate Weight matrix
    weight_func = weight_func.lower()
    if 'var' in weight_func:
        if print_msg:
            print('convert coherence to weight using inverse of phase variance')
            print('    with phase PDF for distributed scatterers from Tough et al. (1995)')
        weight = 1.0 / coherence2phase_variance(coh_data, L, print_msg=print_msg)

    elif any(i in weight_func for i in ['coh', 'lin']):
        if print_msg:
            print('use coherence as weight directly (Perissin & Wang, 2012; Tong et al., 2016)')
        weight = coh_data

    elif any(i in weight_func for i in ['fim', 'fisher']):
        if print_msg:
            print('convert coherence to weight using Fisher Information Index (Seymour & Cumming, 1994)')
        weight = coherence2fisher_info_index(coh_data, L)

    elif weight_func in ['no', 'sbas', 'uniform']:
        weight = None

    else:
        raise Exception('Un-recognized weight function: %s' % weight_func)

    # weight data is calculated in float64 for numerical precision
    # and saved in float32 after calculation for computing speed efficiency.
    if weight is not None:
        weight = np.array(weight, np.float32)
    del coh_data
    return weight

