############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import os
import sys
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
        from matplotlib import pyplot as plt
        plt.switch_backend('Agg')

    return inps


####################################################################################
def main(iargs=None):
    from ..s1ab_range_bias import estimate_s1ab_range_bias, write_s1ab_bias_file, correct_s1ab_range_bias, plot_s1ab_range_bias_est
    
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


####################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
