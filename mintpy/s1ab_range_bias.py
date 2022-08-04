#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Apr 2022                           #
############################################################


import os
import sys
import numpy as np
from matplotlib import pyplot as plt, ticker, colors

from mintpy.objects import timeseries
from mintpy.utils import readfile, writefile, s1_utils, plot as pp
from mintpy.utils.arg_utils import create_argument_parser



####################################################################################
REFERENCE = """reference:
  Yunjun, Z., Fattahi, H., Pi, X., Rosen, P., Simons, M., Agram, P., & Aoki, Y. (2022). Range
    Geolocation Accuracy of C-/L-band SAR and its Implications for Operational Stack Coregistration.
    IEEE Trans. Geosci. Remote Sens., 60, doi:10.1109/TGRS.2022.3168509.
"""

EXAMPLE = """example:
  # Requires a text file named "SAFE_files.txt" containing all Sentinel-1 SAFE filenames.
  # It is generated in ISCE-2/topsStack by default, and could be generated as below if missing:
  # ls ./SLC > SAFE_files.txt

  # 1. compute the S1A/B range bias
  # based on partially corrected TS file, for a more accurate estimation
  s1ab_range_bias.py timeseriesRg_SET_ERA5.h5 -a compute
  s1ab_range_bias.py timeseriesRg_SET_ERA5.h5 -a compute -b data
  s1ab_range_bias.py timeseriesRg_SET_ERA5.h5 -a compute -b data --force
  s1ab_range_bias.py timeseriesRg_SET_ERA5.h5 -a compute -b data --nodisplay

  # 2. correct for the S1A/B range bias [from the 1st/raw TS file]
  s1ab_range_bias.py timeseriesRg.h5 -a correct
"""

def create_parser(subparsers=None):
    synopsis = 'Sentinel-1 A/B range bias correction'
    epilog = REFERENCE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    # input/output files
    parser.add_argument('ts_file', help='Range offset timeseries file to be corrrected, e.g. timeseriesRg_SET_ERA5.h5.')
    parser.add_argument('-g', '--geom', '--geometry', dest='geom_file', help='geometry file including datasets:\nheight')
    parser.add_argument('-m', '--mask', dest='mask_file', help='mask file')

    parser.add_argument('-s', '--safe-list', dest='safe_list_file',
                        help='path to the SAFE_files.txt file, default: in the parent dir of mintpy work dir.')
    parser.add_argument('-o', '--outfile', dest='ts_cor_file',
                        help='Output file name for corrected time-series. Default: add "_S1Bias" suffix.')

    # config
    parser.add_argument('-a', '--action', dest='action', choices={'compute', 'correct'}, default='compute',
                        help='Action to be executed:\n'
                             'compute - estimate the S1A/B range bias and write to HDF5 file.\n'
                             'correct - correct the input TS file using the bias file.')
    parser.add_argument('-b','--method','--bias-method', dest='bias_method', choices={'hardwired', 'data'}, default='hardwired',
                        help='Bias estimation method (default: %(default)s):\n'
                             'hardwired - use hardwired values from section VII-A in Yunjun et al. (2022)\n'
                             'data      - estimate from the input TS file, using the same method as in Yunjun et al. (2022)')
    parser.add_argument('--force', dest='force', action='store_true', help='Force to re-generate the S1Bias.h5 file.')

    # figure
    fig = parser.add_argument_group('Plot the bias estimation result', 'For "--bias-method data" ONLY')
    fig.add_argument('--save', dest='save_fig', action='store_true', help='save the figure')
    fig.add_argument('--nodisplay', dest='disp_fig', action='store_false', help='save and do not display the figure')

    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    inps.mintpy_dir = os.path.dirname(inps.ts_file)

    # --geom
    if not inps.geom_file:
        inps.geom_file = os.path.join(inps.mintpy_dir, 'inputs', 'geometryRadar.h5')
    if not os.path.isfile(inps.geom_file):
        raise FileNotFoundError(f'No geometry file found in: {inps.geom_file}!')

    # --mask
    if not inps.mask_file:
        inps.mask_file = os.path.join(inps.mintpy_dir, 'maskResInv.h5')
    if not os.path.isfile(inps.mask_file):
        inps.mask_file = None

    # --save/nodisplay
    if not inps.disp_fig:
        inps.save_fig = True
    if not inps.disp_fig:
        plt.switch_backend('Agg')

    return inps


####################################################################################
def estimate_s1ab_range_bias(ts_file, mask_file=None, safe_list_file=None):
    """Estimate the S1A/B range bias based on the time series file.

    Parameters: ts_file        - str, path of the time series range offset file
                mask_file      - str, path of the mask file, e.g., maskResInv.h5 file.
                safe_list_file - str, path of the SAFE_files.txt file
    Returns:    bias_list      - list of float32, median bias per subswath in meters
                bias_est       - 2D np.ndarray in size of (length, width) in float32, bias in meters    
                mask_list      - list of 2D np.ndarray in size of (length, width) in bool
    """
    mintpy_dir = os.path.dirname(ts_file)
    print(f'Estimating S1A/B range bias from file: {ts_file}')

    # read time series
    ts_data = readfile.read(ts_file)[0]
    length, width = ts_data.shape[-2:]

    # read mask
    mask_file = mask_file if mask_file else os.path.join(mintpy_dir, 'maskResInv.h5')
    if mask_file and os.path.isfile(mask_file):
        mask = readfile.read(mask_file)[0]
    else:
        mask = np.ones((length, width), dtype=np.bool_)

    # estimate bias - 2D map
    bias_est_poi = s1_utils.estimate_s1ab_bias(
        mintpy_dir,
        ts_data[:, mask],
        safe_list_file=safe_list_file)[0]

    if bias_est_poi is None:
        print('Exit without estimating S1A/B range bias from the input time series.')
        return None, None, None

    bias_est = np.ones((length, width), dtype=np.float32) * np.nan
    bias_est[mask] = bias_est_poi

    # estimate bias - median
    geom_file = os.path.join(mintpy_dir, 'inputs', 'geometryRadar.h5')
    flag = readfile.read(geom_file, datasetName='height')[0] != 0
    mask_list = s1_utils.get_subswath_masks(flag, cut_overlap_in_half=False)[:3]
    bias_list = [np.nanmedian(bias_est[x]) for x in mask_list]
    print('IW1 : {:.3f} m'.format(bias_list[0]))
    print('IW2 : {:.3f} m'.format(bias_list[1]))
    print('IW3 : {:.3f} m'.format(bias_list[2]))

    return bias_list, bias_est, mask_list


def plot_s1ab_range_bias_est(bias_list, bias_est, mask_list, out_dir=None,
                             save_fig=False, disp_fig=True):
    """Plot the S1A/B range bias estimation results.

    Parameters: bias_list - list of float, mean S1A/B range bias in meters
                bias_est  - 2D np.ndarray in float32, pixelwised S1A/B range bias in meters
                mask_list - list of 2D np.ndarray in bool, mask array for IW1/2/3
                out_dir   - str, output directory for the plotted figure.
    """
    vmin, vmax = 7, 14
    font_size = 12
    out_dir = out_dir if out_dir else os.getcwd()

    ## figure 1 - map
    fig_size = pp.auto_figure_size(ds_shape=bias_est.shape, disp_cbar=True, print_msg=True)
    fig, ax = plt.subplots(figsize=fig_size)
    cmap = colors.LinearSegmentedColormap.from_list('magma_t', plt.get_cmap('magma')(np.linspace(0.3, 1.0, 100)))
    im = ax.imshow(bias_est*100., cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest')

    # axis format
    ax.tick_params(which='both', direction='out', bottom=True, top=False, left=True, right=False)
    ax.set_xlabel('Range [pixel]', fontsize=font_size)
    ax.set_ylabel('Azimuth [pixel]', fontsize=font_size)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label('S1A/B range bias [cm]', fontsize=font_size)

    # output
    if save_fig:
        out_fig = os.path.join(out_dir, 's1ab_range_bias_map.pdf')
        print('save figure to file', out_fig)
        plt.savefig(out_fig, bbox_inches='tight', transparent=True, dpi=300)

    ## figure 2 - histogram
    clist = [cmap((x*100. - vmin) / (vmax - vmin)) for x in bias_list]
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=[6, 2])
    for bias, mask, c in zip(bias_list, mask_list, clist):
        ax.hist(bias_est[mask].flatten()*100, bins=70, range=(vmin, vmax), density=False, alpha=0.7, color=c)
        ax.axvline(bias*100, color='k')
    # plot median value
    ax.tick_params(which='both', direction='out', bottom=True, top=True, left=True, right=True)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.set_xlim(vmin, vmax)
    ax.set_xlabel('S1A/B range bias [cm]')
    ax.set_ylabel('# of pixels')
    fig.tight_layout()

    # output
    if save_fig:
        out_fig = os.path.join(out_dir, 's1ab_range_bias_hist.pdf')
        print('save figure to file', out_fig)
        plt.savefig(out_fig, bbox_inches='tight', transparent=True, dpi=300)

    if disp_fig:
        print('showing ....')
        plt.show()
    else:
        plt.close()

    return


def write_s1ab_bias_file(bias_file, bias_list, geom_file, force=False):
    """Write estimated S1A/B range bias to HDF5 file.

    Parameters: bias_file - str, path to the S1A/B range bias file
                bias_list - list of float, constant S1A/B range bias per IW1/2/3
                geom_file - str, path to the geometry file
                force     - bool, overwrite existing bias file.
    Returns:    bias_file - str, path to the S1A/B range bias file
    """
    # run or skip
    if os.path.isfile(bias_file) and not force:
        print(f'S1Bias file exists in: {bias_file}, skip re-writing.')
        return bias_file

    # get the list of masks for IW1/2/3
    flag = readfile.read(geom_file, datasetName='height')[0] != 0
    mask_list = s1_utils.get_subswath_masks(flag, cut_overlap_in_half=False)[:3]

    bias_mat = np.zeros(flag.shape, dtype=np.float32) * np.nan
    for bias, mask in zip(bias_list, mask_list):
        bias_mat[mask] = bias

    # write file
    atr = readfile.read_attribute(geom_file)
    atr['FILE_TYPE'] = 'offset'
    atr['UNIT'] = 'm'
    print('writing S1A/B range bias to file: {}'.format(bias_file))
    writefile.write(bias_mat, out_file=bias_file, metadata=atr)

    return bias_file


def correct_s1ab_range_bias(ts_file, bias_file, ts_cor_file=None, safe_list_file=None):
    """Correct input time series for the S1A/B range bias.

    Parameters: ts_file        - str, path to the range offset time series file
                bias_file      - str, path to the S1A/B range bias file
                ts_cor_file    - str, path to the corrected range offset time series file
                safe_list_file - str, path to the SAFE_files.txt file
    Returns:    ts_cor_file    - str, path to the corrected range offset time series file
    """

    if not os.path.isfile(bias_file):
        msg = f'No bias file found in: {bias_file}!'
        msg += '\nRe-run with "--action compute" to generate it.'
        raise FileNotFoundError(msg)

    date_list = timeseries(ts_file).get_date_list()
    num_date = len(date_list)

    # date info for Sentinel-1B
    mintpy_dir = os.path.dirname(os.path.dirname(bias_file))
    s1b_date_list_file = s1_utils.get_s1ab_date_list_file(mintpy_dir, safe_list_file)[1]
    s1b_date_list = np.loadtxt(s1b_date_list_file, dtype=str).tolist()
    s1b_flag = np.array([x in s1b_date_list for x in date_list], dtype=np.bool_)

    # read data
    ts_data = readfile.read(ts_file)[0].reshape(num_date, -1)
    bias = readfile.read(bias_file)[0].flatten()

    # correct bias
    mask = ts_data == 0
    ts_data[s1b_flag] -= np.tile(bias.reshape(1, -1), (np.sum(s1b_flag), 1))
    ts_data[mask] = 0                    # Do not change zero value in the input TS file
    ts_data[:, np.isnan(bias)] = np.nan  # set to nan for pixels with nan in bias file 

    # write file
    if not ts_cor_file:
        ts_cor_file = '{}_S1Bias.h5'.format(os.path.splitext(ts_file)[0])
    atr = readfile.read_attribute(ts_file)
    length = int(atr['LENGTH'])
    width = int(atr['WIDTH'])
    writefile.write(ts_data.reshape(num_date, length, width),
                    out_file=ts_cor_file,
                    metadata=atr,
                    ref_file=ts_file)

    return ts_cor_file



####################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    # default bias file path
    inps.bias_file = os.path.join(os.path.dirname(inps.geom_file), 'S1Bias.h5')

    # calculate the S1A/B range bias
    if inps.action == 'compute':
        if inps.bias_method == 'hardwired':
            # option 1 - use the hardwired value from section VII-A in Yunjun et al. (2022)
            bias_list = [0.087, 0.106, 0.123]   # m
            print('Used hardwired S1A/B range bias values from Yunjun et al. (2022):')
            print('IW1 : {:.3f} m'.format(bias_list[0]))
            print('IW2 : {:.3f} m'.format(bias_list[1]))
            print('IW3 : {:.3f} m'.format(bias_list[2]))

        else:
            # option 2 - estimate from the time series of its dataset itself
            # estimate optimal (median) value for each subswath from SenDT156
            bias_list, bias_est, mask_list = estimate_s1ab_range_bias(
                ts_file=inps.ts_file,
                mask_file=inps.mask_file,
                safe_list_file=inps.safe_list_file)

            # plot the estimation result
            if bias_list:
                plot_s1ab_range_bias_est(
                    bias_list,
                    bias_est,
                    mask_list,
                    out_dir=os.path.dirname(inps.ts_file),
                    save_fig=inps.save_fig,
                    disp_fig=inps.disp_fig)

        # write S1Bias.h5 file
        if bias_list:
            write_s1ab_bias_file(
                bias_file=inps.bias_file,
                bias_list=bias_list,
                geom_file=inps.geom_file,
                force=inps.force)

    # correct time series range offset file
    elif inps.action == 'correct':
        correct_s1ab_range_bias(
            ts_file=inps.ts_file,
            bias_file=inps.bias_file,
            ts_cor_file=inps.ts_cor_file,
            safe_list_file=inps.safe_list_file)

    return


####################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
