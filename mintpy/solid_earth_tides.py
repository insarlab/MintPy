#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Sep 2020                           #
############################################################
# Recomend import:
#   from mintpy import solid_earth_tides as SET


import os
import sys
import time
import argparse
import warnings
import numpy as np

try:
    import pysolid
except ImportError:
    raise ImportError('Can not import pysolid! Check https://github.com/insarlab/PySolid.')

from mintpy.objects import timeseries
from mintpy.objects.resample import resample
from mintpy.defaults.template import get_template_content
from mintpy.utils import (
    ptime,
    readfile,
    writefile,
    utils as ut,
    attribute as attr,
)


###############################################################
TEMPLATE = get_template_content('correct_SET')

EXAMPLE = """example:
  solid_earth_tides.py timeseries.h5 -g inputs/geometryRadar.h5
  solid_earth_tides.py timeseries.h5 -g inputs/geometryGeo.h5
  solid_earth_tides.py geo/geo_timeseries_ERA5_demErr.h5 -g geo/geo_geometryRadar.h5
"""

REFERENCE = """reference:
  Milbert, D., SOLID EARTH TIDE, http://geodesyworld.github.io/SOFTS/solid.htm, Accessd 2020 September 6.
  Fattahi, H., Z. Yunjun, X. Pi, P. S. Agram, P. Rosen, and Y. Aoki (2020), Absolute geolocation of SAR 
    Big-Data: The first step for operational InSAR time-series analysis, AGU Fall Meeting 2020, 1-17 Dec 2020.
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Solid Earth tides correction',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog='{}\n{}\n{}'.format(REFERENCE, TEMPLATE, EXAMPLE))

    parser.add_argument('dis_file', help='timeseries HDF5 file, i.e. timeseries.h5')
    parser.add_argument('-g','--geomtry', dest='geom_file', type=str, required=True,
                        help='geometry file including incidence/azimuthAngle.')

    parser.add_argument('--verbose', dest='verbose', action='store_true', help='Verbose message.')
    parser.add_argument('--update', dest='update_mode', action='store_true', help='Enable update mode.')

    # output
    parser.add_argument('--set-file', dest='set_file', help='line-of-sight solid earth tide file name')
    parser.add_argument('-o', dest='cor_dis_file', help='Output file name for the corrected timeseries.')

    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check coordinates of time-series / geometry files
    atr1 = readfile.read_attribute(inps.dis_file)
    atr2 = readfile.read_attribute(inps.geom_file)
    geo_or_rdr1 = 'Y_FIRST' in atr1.keys()
    geo_or_rdr2 = 'Y_FIRST' in atr2.keys()
    if geo_or_rdr1 != geo_or_rdr2:
        msg = 'input time-series and geometry file do not have the same coordinates!'
        msg += 'time-series file in geo-coordinates: {}'.format(geo_or_rdr1)
        msg += 'geometry    file in geo-coordinates: {}'.format(geo_or_rdr2)
        raise ValueError(msg)

    # default SET filename
    if not inps.set_file:
        geom_dir = os.path.dirname(inps.geom_file)
        fname = 'SET.h5'
        if os.path.basename(inps.dis_file).startswith('geo_'):
            fname = 'geo_SET.h5'
        inps.set_file = os.path.join(geom_dir, fname)

    # default corrected time-series filename
    if not inps.cor_dis_file:
        dis_dir = os.path.dirname(inps.dis_file)
        fbase, fext = os.path.splitext(os.path.basename(inps.dis_file))
        inps.cor_dis_file = os.path.join(dis_dir, '{}_SET{}'.format(fbase, fext))

    return inps


###############################################################
def prepare_los_geometry(geom_file):
    """Prepare LOS geometry data/info in geo-coordinates
    Parameters: geom_file  - str, path of geometry file
    Returns:    inc_angle  - 2D np.ndarray, incidence angle in radians
                head_angle - 2D np.ndarray, heading   angle in radians
                atr        - dict, metadata in geo-coordinate
    """

    print('read/prepare LOS geometry from file: {}'.format(geom_file))
    atr = readfile.read_attribute(geom_file)

    print('read incidence / azimuth angle')
    inc_angle = readfile.read(geom_file, datasetName='incidenceAngle')[0]
    az_angle  = readfile.read(geom_file, datasetName='azimuthAngle')[0]

    # geocode inc/az angle data if in radar-coord
    if 'Y_FIRST' not in atr.keys():
        print('-'*50)
        print('geocoding the incidence / azimuth angle ...')
        res_obj = resample(lut_file=geom_file, src_file=geom_file)
        res_obj.open()
        res_obj.prepare()

        # resample data
        box = res_obj.src_box_list[0]
        inc_angle = res_obj.run_resample(src_data=inc_angle[box[1]:box[3], box[0]:box[2]])
        az_angle  = res_obj.run_resample(src_data=az_angle[box[1]:box[3], box[0]:box[2]])

        # update attribute
        atr = attr.update_attribute4radar2geo(atr, res_obj=res_obj)

    # azimuth angle --> heading angle
    head_angle = ut.azimuth2heading_angle(az_angle)

    # unit: degree to radian
    inc_angle *= np.pi / 180.
    head_angle *= np.pi / 180.

    return inc_angle, head_angle, atr


def calc_solid_earth_tides_timeseries(date_list, geom_file, out_file, update_mode=True, verbose=False):
    """Calculate the time-series of solid Earth tides in LOS direction.
    Parameters: date_list - list of str, dates in YYYYMMDD
                geom_file - str, path of geometry file in geo coordinates
                out_file  - str, output time-sereis file
    Returns:    out_file  - str, output time-sereis file
    """

    if update_mode and os.path.isfile(out_file):
        print('update mode: ON')
        print('skip re-calculating and use existing file: {}'.format(out_file))
        return out_file

    # prepare LOS geometry: geocoding if in radar-coordinates
    inc_angle, head_angle, atr_geo = prepare_los_geometry(geom_file)

    # get LOS unit vector
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        unit_vec = [
            np.sin(inc_angle) * np.cos(head_angle) * -1,
            np.sin(inc_angle) * np.sin(head_angle),
            np.cos(inc_angle),
        ]

    # initiate data matrix
    num_date = len(date_list)
    length = int(atr_geo['LENGTH'])
    width = int(atr_geo['WIDTH'])
    ts_tide = np.zeros((num_date, length, width), dtype=np.float32)

    # loop for calc
    print('\n'+'-'*50)
    print('calculating solid Earth tides using solid.for (D. Milbert, 2018) ...')
    prog_bar = ptime.progressBar(maxValue=num_date, print_msg=not verbose)
    for i in range(num_date):
        date_str = date_list[i]

        # calculate tide in ENU direction
        (tide_e,
         tide_n,
         tide_u) = pysolid.calc_solid_earth_tides_grid(date_str, atr_geo,
                                                       display=False,
                                                       verbose=verbose)

        # convert ENU to LOS direction
        # sign convention: positive for motion towards satellite
        ts_tide[i,:,:] = (tide_e * unit_vec[0]
                          + tide_n * unit_vec[1]
                          + tide_u * unit_vec[2])

        prog_bar.update(i+1, suffix='{} ({}/{})'.format(date_list[i], i+1, num_date))
    prog_bar.close()

    # radar-coding if input in radar-coordinates
    atr = readfile.read_attribute(geom_file)
    if 'Y_FIRST' not in atr.keys():
        print('radar-coding the LOS tides time-series ...')
        res_obj = resample(lut_file=geom_file)
        res_obj.open()
        res_obj.src_meta = atr_geo
        res_obj.prepare()

        # resample data
        box = res_obj.src_box_list[0]
        ts_tide = res_obj.run_resample(src_data=ts_tide[:,
                                                        box[1]:box[3],
                                                        box[0]:box[2]])

    ## output
    # attribute
    atr['FILE_TYPE'] = 'timeseries'
    atr['UNIT'] = 'm'
    for key in ['REF_Y', 'REF_X', 'REF_DATE']:
        if key in atr.keys():
            atr.pop(key)

    # write
    ds_dict = {}
    ds_dict['date'] = np.array(date_list, dtype=np.string_)
    ds_dict['timeseries'] = ts_tide
    writefile.write(ds_dict, out_file=out_file, metadata=atr)

    return out_file


def correct_timeseries(dis_file, set_file, cor_dis_file):
    """Correct time-series for the solid Earth tides."""
    # diff.py can handle different reference in space and time
    # between the absolute solid Earth tides and the double referenced time-series
    print('\n------------------------------------------------------------------------------')
    print('correcting relative delay for input time-series using diff.py')
    from mintpy import diff

    iargs = [dis_file, set_file, '-o', cor_dis_file]
    print('diff.py', ' '.join(iargs))
    diff.main(iargs)
    return cor_dis_file


###############################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    start_time = time.time()

    # calc SET - prepare
    print('read date list from file: {}'.format(inps.dis_file))
    date_list = timeseries(inps.dis_file).get_date_list()

    # calc SET - run
    calc_solid_earth_tides_timeseries(date_list, inps.geom_file,
                                      out_file=inps.set_file,
                                      update_mode=inps.update_mode,
                                      verbose=inps.verbose)

    # correct SET
    correct_timeseries(dis_file=inps.dis_file,
                       set_file=inps.set_file,
                       cor_dis_file=inps.cor_dis_file)

    m, s = divmod(time.time() - start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.\n'.format(m, s))
    return inps.cor_dis_file

###############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
