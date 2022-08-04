#! /usr/bin/env python
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Sara Mirzaee, Zhang Yunjun, Bhuvan Varugu, 2018  #
############################################################


import os
import sys
import re
import h5py
import numpy as np
from skimage.transform import resize
from scipy.interpolate import RegularGridInterpolator as RGI

from mintpy.objects import timeseries
from mintpy.utils import ptime, readfile, writefile, utils as ut
from mintpy.utils.arg_utils import create_argument_parser


############################################################################
REFERENCE = """references:
  Yu, C., Li, Z., Penna, N. T., & Crippa, P. (2018). Generic atmospheric correction model for Interferometric
    Synthetic Aperture Radar observations. Journal of Geophysical Research: Solid Earth, 123(10), 9202-9222.
  Yu, C., Li, Z., & Penna, N. T. (2018). Interferometric synthetic aperture radar atmospheric correction
    using a GPS-based iterative tropospheric decomposition model. Remote Sensing of Environment, 204, 109-121.
"""

DIR_DEMO = """--dir ./GACOS
  20060624.ztd
  20060624.ztd.rsc
  20061225.ztd
  20061225.ztd.rsc
  ...
  OR
  20060624.ztd.tif
  20061225.ztd.tif
  ...
"""

EXAMPLE = """example:
  tropo_gacos.py -f timeseries.h5 -g inputs/geometryRadar.h5 --dir ./GACOS
  tropo_gacos.py -f geo/geo_timeseries.h5 -g geo/geo_geometryRadar.h5 --dir ./GACOS
"""


def create_parser(subparsers=None):
    synopsis = 'Tropospheric correction using GACOS (http://www.gacos.net) delays'
    epilog = REFERENCE + '\n' + DIR_DEMO + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('-f', '--file', dest='dis_file', required=True,
                        help='timeseries HDF5 file, i.e. timeseries.h5')
    parser.add_argument('-g', '--geom', dest='geom_file', required=True,
                        help='geometry file.')
    parser.add_argument('--dir','--GACOS-dir', dest='GACOS_dir', default='./GACOS',
                        help='directory to downloaded GACOS delays data (default: %(default)s).')
    parser.add_argument('-o', dest='cor_dis_file',
                        help='Output file name for trospheric corrected timeseries.')

    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    inps.GACOS_dir = os.path.abspath(inps.GACOS_dir)
    print('Use GACOS products at directory:', inps.GACOS_dir)

    # check input files - existance
    for key in ['dis_file', 'geom_file']:
        fname = vars(inps)[key]
        if fname and not os.path.isfile(fname):
            raise FileNotFoundError('input file not exist: {}'.format(fname))

    # check input files - processors & coordinates
    atr1 = readfile.read_attribute(inps.dis_file)
    atr2 = readfile.read_attribute(inps.geom_file)
    coord1 = 'geo' if 'Y_FIRST' in atr1.keys() else 'radar'
    coord2 = 'geo' if 'Y_FIRST' in atr2.keys() else 'radar'
    proc = atr1.get('PROCESSOR', 'isce')

    if coord1 == 'radar' and proc in ['gamma', 'roipac']:
        msg = 'Radar-coded file from {} is NOT supported!'.format(proc)
        msg += '\n    Try to geocode the time-series and geometry files and re-run with them instead.'
        raise ValueError(msg)

    if coord1 != coord2:
        n = max(len(os.path.basename(i)) for i in [inps.dis_file, inps.geom_file])
        msg = 'Input time-series and geometry file are NOT in the same coordinate!'
        msg += '\n    file {f:<{n}} coordinate: {c}'.format(f=os.path.basename(inps.dis_file),  n=n, c=coord1)
        msg += '\n    file {f:<{n}} coordinate: {c}'.format(f=os.path.basename(inps.geom_file), n=n, c=coord2)
        raise ValueError(msg)

    # default output filenames
    inps.tropo_file = os.path.join(os.path.dirname(inps.geom_file), 'GACOS.h5')
    if not inps.cor_dis_file:
        inps.cor_dis_file = inps.dis_file.split('.')[0] + '_GACOS.h5'

    return inps


############################################################################
def get_delay_geo(ztd_file, atr, cos_inc_angle):
    """calc single path tropo delay in line-of-sight direction

    Parameters: ztd_file      - str, path of zenith delay file
                atr           - dict, dictionary of attribute for output file
                cos_inc_angle - 2D np.ndarray in float32, cos(inc_angle)
    Returns:    delay         - 2D np.ndarray in float32, LOS delay
    """

    # get geo_box from ts_file attributes
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    geo_box = ut.coordinate(atr).box_pixel2geo(pixel_box=(0, 0, width, length))

    # geo_box --> pix_box in ztd file
    atr_ztd = readfile.read_attribute(ztd_file)
    pix_box = ut.coordinate(atr_ztd).box_geo2pixel(geo_box)

    # read ztd file
    delay = readfile.read(ztd_file, box=pix_box)[0]

    # interpolate/resample into the same resolution as ts_file
    delay = resize(delay, (length, width),
                   order=1,
                   mode='constant',
                   anti_aliasing=True,
                   preserve_range=True)

    # project from zenith to line-of-sight
    delay /= cos_inc_angle

    # reverse the sign for consistency between different phase correction steps/methods
    delay *= -1

    return delay


def get_delay_radar(ztd_file, cos_inc_angle, pts_new):
    """calc single path tropo delay in line-of-sight direction

    Parameters: ztd_file      - str, path of zenith delay file
                cos_inc_angle - 2D np.ndarray in (len, wid) in float32, cos(inc_angle)
                pts_new       - 2D np.ndarray in (len*wid, 2) in float32
    Returns:    delay         - 2D np.ndarray in float32, LOS delay
    """
    # read ztd file
    delay_ztd, atr_ztd = readfile.read(ztd_file)
    # flip to be consistent with the reversed lats
    delay_ztd = np.flipud(delay_ztd)

    # pixel coordinates in ztd file
    lats, lons = ut.get_lat_lon(atr_ztd, dimension=1)
    # set lats in ascending order as required by RGI
    lats = np.flipud(lats)
    pts_ztd = ((lats.flatten(),
                lons.flatten()))

    # resample in pts_new coordinates
    RGI_func = RGI(pts_ztd, delay_ztd,
                   method='nearest',
                   bounds_error=False,
                   fill_value=0)
    delay = RGI_func(pts_new)
    delay = delay.reshape(cos_inc_angle.shape)

    # project from zenith to line-of-sight
    delay /= cos_inc_angle

    # reverse the sign for consistency between different phase correction steps/methods
    delay *= -1

    return delay


def calculate_delay_timeseries(tropo_file, dis_file, geom_file, GACOS_dir):
    """calculate delay time-series and write to HDF5 file"""

    ## get list of dates
    atr = readfile.read_attribute(dis_file)
    ftype = atr['FILE_TYPE']
    if ftype == 'timeseries':
        date_list = timeseries(dis_file).get_date_list()

    elif ftype == '.unw':
        date12 = readfile.read_attribute(dis_file)['DATE12']
        date_list = ptime.yyyymmdd(date12.split('-'))

    else:
        raise ValueError('un-supported displacement file type: {}'.format(ftype))

    # list of dates --> list of ztd files
    ztd_files = []
    flag = np.ones(len(date_list), dtype=np.bool_)
    for i, date_str in enumerate(date_list):
        fnames = [os.path.join(GACOS_dir, '{}{}'.format(date_str, fext)) for fext in ['.ztd', '.ztd.tif']]
        fnames = [f for f in fnames if os.path.exists(f)]
        if len(fnames) > 0:
            ztd_files.append(fnames[0])
        else:
            print('WARNING: NO ztd file found for {}! Continue without it.'.format(date_str))
            flag[i] = False

    # update date_list to be consistent with ztd_files
    if np.any(flag == 0):
        date_list = np.array(date_list)[flag].tolist()

    ## update_mode
    def get_dataset_size(fname):
        atr = readfile.read_attribute(fname)
        return (atr['LENGTH'], atr['WIDTH'])

    def run_or_skip(ztd_files, tropo_file, geom_file):
        print('update mode: ON')
        print('output file: {}'.format(tropo_file))
        flag = 'skip'

        # check existance and modification time
        if ut.run_or_skip(out_file=tropo_file, in_file=ztd_files, print_msg=False) == 'run':
            flag = 'run'
            print('1) output file either do NOT exist or is NOT newer than all ZTD files.')

        else:
            print('1) output file exists and is newer than all ZTD files.')

            # check dataset size in space / time
            date_list = [str(re.findall('\d{8}', i)[0]) for i in ztd_files]
            if (get_dataset_size(tropo_file) != get_dataset_size(geom_file)
                    or any(i not in timeseries(tropo_file).get_date_list() for i in date_list)):
                flag = 'run'
                print(('2) output file does NOT have the same len/wid as the geometry file {}'
                       ' or does NOT contain all dates').format(geom_file))
            else:
                print('2) output file has the same len/wid as the geometry file and contains all dates')

                # check if output file is fully written
                with h5py.File(tropo_file, 'r') as f:
                    if np.all(f['timeseries'][-1,:,:] == 0):
                        flag = 'run'
                        print('3) output file is NOT fully written.')
                    else:
                        print('3) output file is fully written.')

        # result
        print('run or skip: {}'.format(flag))
        return flag

    if run_or_skip(ztd_files, tropo_file, geom_file) == 'skip':
        return


    ## prepare output file

    # metadata
    atr['FILE_TYPE'] = 'timeseries'
    atr['UNIT'] = 'm'

    # remove metadata related with double reference
    # because absolute delay is calculated and saved
    for key in ['REF_DATE', 'REF_X', 'REF_Y', 'REF_LAT', 'REF_LON']:
        if key in atr.keys():
            atr.pop(key)

    # instantiate time-series
    length, width = int(atr['LENGTH']), int(atr['WIDTH'])
    num_date = len(date_list)
    dates = np.array(date_list, dtype=np.string_)
    ds_name_dict = {
        "date"       : [dates.dtype, (num_date,), dates],
        "timeseries" : [np.float32,  (num_date, length, width), None],
    }
    writefile.layout_hdf5(tropo_file, ds_name_dict, metadata=atr)


    ## calculate phase delay

    # read geometry
    print('read incidenceAngle from file: {}'.format(geom_file))
    inc_angle = readfile.read(geom_file, datasetName='incidenceAngle')[0]
    cos_inc_angle = np.cos(inc_angle * np.pi / 180.0)

    if 'Y_FIRST' in atr.keys():
        # No need for data in geo-coordinates
        pts_new = None

    else:
        # Get pixel lat/lon for data in radar-coordinates
        print('get pixel coordinates in geometry file')
        lats, lons = ut.get_lat_lon(atr, geom_file)
        pts_new = np.hstack((lats.reshape(-1, 1),
                             lons.reshape(-1, 1)))

    # loop for date-by-date IO
    prog_bar = ptime.progressBar(maxValue=num_date)
    for i in range(num_date):
        date_str = date_list[i]
        ztd_file = ztd_files[i]

        # calc delay
        if 'Y_FIRST' in atr.keys():
            delay = get_delay_geo(ztd_file, atr, cos_inc_angle)

        else:
            delay = get_delay_radar(ztd_file, cos_inc_angle, pts_new)

        # write delay to file
        block = [i, i+1, 0, length, 0, width]
        writefile.write_hdf5_block(tropo_file,
                                   data=delay,
                                   datasetName='timeseries',
                                   block=block,
                                   print_msg=False)

        prog_bar.update(i + 1, suffix=os.path.basename(ztd_file))
    prog_bar.close()

    return tropo_file


def correct_timeseries(dis_file, tropo_file, cor_dis_file):
    # diff.py can handle different reference in space and time
    # between the absolute tropospheric delay and the double referenced time-series
    print('\n------------------------------------------------------------------------------')
    print('correcting relative delay for input time-series using diff.py')
    from mintpy import diff

    iargs = [dis_file, tropo_file, '-o', cor_dis_file, '--force']
    print('diff.py', ' '.join(iargs))
    diff.main(iargs)
    return cor_dis_file


def correct_single_ifgram(dis_file, tropo_file, cor_dis_file):
    print('\n------------------------------------------------------------------------------')
    print('correcting relative delay for input interferogram')

    print('read data from {}'.format(dis_file))
    data, atr = readfile.read(dis_file, datasetName='phase')
    date1, date2 = ptime.yyyymmdd(atr['DATE12'].split('-'))

    print('calc tropospheric delay for {}-{} from {}'.format(date1, date2, tropo_file))
    tropo  = readfile.read(tropo_file, datasetName='timeseries-{}'.format(date2))[0]
    tropo -= readfile.read(tropo_file, datasetName='timeseries-{}'.format(date1))[0]
    tropo *= -4. * np.pi / float(atr['WAVELENGTH'])

    print('write corrected data to {}'.format(cor_dis_file))
    writefile.write(data - tropo, cor_dis_file, atr)
    return cor_dis_file


############################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    # calculate tropo delay and savee to h5 file
    calculate_delay_timeseries(
        tropo_file=inps.tropo_file,
        dis_file=inps.dis_file,
        geom_file=inps.geom_file,
        GACOS_dir=inps.GACOS_dir)

    # correct tropo delay from dis time-series
    ftype = readfile.read_attribute(inps.dis_file)['FILE_TYPE']
    if ftype == 'timeseries':
        correct_timeseries(
            dis_file=inps.dis_file,
            tropo_file=inps.tropo_file,
            cor_dis_file=inps.cor_dis_file)

    elif ftype == '.unw':
        correct_single_ifgram(
            dis_file=inps.dis_file,
            tropo_file=inps.tropo_file,
            cor_dis_file=inps.cor_dis_file)
    else:
        print('input file {} is not timeseries nor .unw, correction is not supported yet.'.format(ftype))

    return


############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
