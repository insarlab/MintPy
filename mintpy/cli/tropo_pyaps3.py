#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import os
import sys

from mintpy.utils.arg_utils import create_argument_parser

###############################################################
REFERENCE = """reference:
  Jolivet, R., R. Grandin, C. Lasserre, M.-P. Doin and G. Peltzer (2011), Systematic InSAR tropospheric
    phase delay corrections from global meteorological reanalysis data, Geophys. Res. Lett., 38, L17311,
    doi:10.1029/2011GL048757
  Jolivet, R., P. S. Agram, N. Y. Lin, M. Simons, M. P. Doin, G. Peltzer, and Z. Li (2014), Improving
    InSAR geodesy using global atmospheric models, Journal of Geophysical Research: Solid Earth, 119(3),
    2324-2341, doi:10.1002/2013JB010588.
  Hersbach, H., Bell, B., Berrisford, P., Hirahara, S., Horányi, A., Muñoz-Sabater, J., et al. (2020).
    The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society, 146(730),
    1999–2049, doi:10.1002/qj.3803
"""

EXAMPLE = """example:
  # download datasets, calculate tropospheric delays and correct time-series file.
  tropo_pyaps3.py -f timeseries.h5 -g inputs/geometryRadar.h5
  tropo_pyaps3.py -f filt_fine.unw -g ../../../mintpy/inputs/geometryRadar.h5

  # download datasets, calculate tropospheric delays
  tropo_pyaps3.py -d date.list         --hour 12 -m ERA5  -g inputs/geometryGeo.h5
  tropo_pyaps3.py -d 20151002 20151003 --hour 12 -m MERRA -g inputs/geometryRadar.h5

  # download datasets (covering the whole world)
  tropo_pyaps3.py -d date.list --hour 12
  tropo_pyaps3.py -d SAFE_files.txt
  # download datasets (covering the area of interest)
  tropo_pyaps3.py -d SAFE_files.txt -g inputs/geometryRadar.h5
"""

SAFE_FILE = """SAFE_files.txt:
    /data/SanAndreasSenDT42/SLC/S1B_IW_SLC__1SDV_20191117T140737_20191117T140804_018968_023C8C_82DC.zip
    /data/SanAndreasSenDT42/SLC/S1A_IW_SLC__1SDV_20191111T140819_20191111T140846_029864_036803_69CA.zip
    ...
"""

DATA_INFO = """Global Atmospheric Models:
  re-analysis_dataset      coverage  temp_resolution  spatial_resolution       latency       assimilation
  --------------------------------------------------------------------------------------------------------
  ERA5(T)  (ECMWF)          global       hourly        0.25 deg (~31 km)   3 months (5 days)    4D-Var
  ERA-Int  (ECMWF)          global       6-hourly      0.75 deg (~79 km)        2 months        4D-Var
  MERRA(2) (NASA Goddard)   global       6-hourly     0.5*0.625 (~50 km)       2-3 weeks        3D-Var
  NARR     (NOAA, working from Jan 1979 to Oct 2014)

Notes for data access:
  For MERRA2, you need an Earthdata account, and pre-authorize the "NASA GESDISC DATA ARCHIVE" application
      following https://disc.gsfc.nasa.gov/earthdata-login.
  For ERA5 from CDS, you need to agree to the Terms of Use of every datasets that you intend to download.
"""

WEATHER_DIR_DEMO = """--weather-dir ~/data/aux
atmosphere/
    /ERA5
        ERA5_N20_N40_E120_E140_20060624_14.grb
        ERA5_N20_N40_E120_E140_20060924_14.grb
        ...
    /MERRA
        merra-20110126-06.nc4
        merra-20110313-06.nc4
        ...
"""


def create_parser(subparsers=None):
    synopsis = 'Tropospheric correction using weather models via PyAPS'
    epilog = REFERENCE + '\n' + DATA_INFO + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('-f', '--file', dest='dis_file',
                        help='timeseries HDF5 file, i.e. timeseries.h5')
    parser.add_argument('-d', '--date-list', dest='date_list', type=str, nargs='*',
                        help='List of dates in YYYYMMDD or YYMMDD format. It can be:\n'
                             'a) list of strings in YYYYMMDD or YYMMDD format OR\n'
                             'b) a text file with the first column as list of date in YYYYMMDD or YYMMDD format OR\n'
                             'c) a text file with Sentinel-1 SAFE filenames\ne.g.: '+SAFE_FILE)
    parser.add_argument('--hour', type=str, help='time of data in HH, e.g. 12, 06')
    parser.add_argument('-o','--output', dest='cor_dis_file',
                        help='Output file name for trospheric corrected timeseries.')

    # delay calculation
    delay = parser.add_argument_group('delay calculation')
    delay.add_argument('-m', '--model', '-s', dest='tropo_model', default='ERA5',
                       choices={'ERA5'}, # {'ERA5', 'MERRA', 'NARR'},
                       help='source of the atmospheric model (default: %(default)s).')
    delay.add_argument('--delay', dest='delay_type', default='comb', choices={'comb', 'dry', 'wet'},
                       help='Delay type to calculate, comb contains both wet and dry delays (default: %(default)s).')

    delay.add_argument('-w', '--dir', '--weather-dir', dest='weather_dir', default='${WEATHER_DIR}',
                       help='parent directory of downloaded weather data file (default: %(default)s).\n' +
                            'e.g.: '+WEATHER_DIR_DEMO)
    delay.add_argument('-g','--geomtry', dest='geom_file', type=str,
                       help='geometry file including height, incidenceAngle and/or latitude and longitude')
    delay.add_argument('--custom-height', dest='custom_height', type=float,
                       help='[for testing] specify a custom height value for delay calculation.')

    delay.add_argument('--tropo-file', dest='tropo_file', type=str,
                       help='tropospheric delay file name')
    delay.add_argument('--verbose', dest='verbose', action='store_true', help='Verbose message.')
    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check + default: -w / --weather-dir option (expand special symbols)
    inps.weather_dir = os.path.expanduser(inps.weather_dir)
    inps.weather_dir = os.path.expandvars(inps.weather_dir)
    if inps.weather_dir == '${WEATHER_DIR}':
        # fallback to current dir if env var WEATHER_DIR is not defined
        inps.weather_dir = './'
    inps.weather_dir = os.path.abspath(inps.weather_dir)

    # check: existence of input files
    for fname in [inps.dis_file, inps.geom_file]:
        if fname and not os.path.isfile(fname):
            raise FileNotFoundError(f'input file not exist: {fname}')

    # check: required options (for date/time): --file OR --date-list
    if (not inps.dis_file and not inps.date_list):
        raise SystemExit('ERROR: --file OR --date-list is required.\n\n'+EXAMPLE)

    # check: --hour option (ensure 2 chars)
    if inps.hour:
        hour = int(float(inps.hour))
        inps.hour = f'{hour:02d}'

    # default: --tropo-file option
    if inps.geom_file and not inps.tropo_file:
        geom_dir = os.path.dirname(inps.geom_file)
        inps.tropo_file = os.path.join(geom_dir, f'{inps.tropo_model}.h5')

    # default: -o / --output option
    if inps.dis_file and not inps.cor_dis_file:
        fbase, fext = os.path.splitext(inps.dis_file)
        inps.cor_dis_file = f'{fbase}_{inps.tropo_model}{fext}'

    return inps


###############################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.tropo_pyaps3 import run_tropo_pyaps3

    # run
    run_tropo_pyaps3(inps)


###############################################################
if __name__ == '__main__':
    main(sys.argv[1:])
