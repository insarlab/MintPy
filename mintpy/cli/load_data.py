############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################

import os
import sys
import time

from mintpy.defaults import auto_path
from mintpy.defaults.template import get_template_content
from mintpy.utils.arg_utils import create_argument_parser


#################################################################
PROCESSOR_LIST = ['isce', 'aria', 'hyp3', 'gmtsar', 'snap', 'gamma', 'roipac', 'cosicorr']

IFG_DSET_NAME2TEMPLATE_KEY = {
    'unwrapPhase'     : 'mintpy.load.unwFile',
    'coherence'       : 'mintpy.load.corFile',
    'connectComponent': 'mintpy.load.connCompFile',
    'wrapPhase'       : 'mintpy.load.intFile',
    'magnitude'       : 'mintpy.load.magFile',
}

ION_DSET_NAME2TEMPLATE_KEY = {
    'unwrapPhase'     : 'mintpy.load.ionUnwFile',
    'coherence'       : 'mintpy.load.ionCorFile',
    'connectComponent': 'mintpy.load.ionConnCompFile',
}

OFF_DSET_NAME2TEMPLATE_KEY = {
    'azimuthOffset'   : 'mintpy.load.azOffFile',
    'azimuthOffsetStd': 'mintpy.load.azOffStdFile',
    'rangeOffset'     : 'mintpy.load.rgOffFile',
    'rangeOffsetStd'  : 'mintpy.load.rgOffStdFile',
    'offsetSNR'       : 'mintpy.load.offSnrFile',
}

GEOM_DSET_NAME2TEMPLATE_KEY = {
    'height'          : 'mintpy.load.demFile',
    'latitude'        : 'mintpy.load.lookupYFile',
    'longitude'       : 'mintpy.load.lookupXFile',
    'azimuthCoord'    : 'mintpy.load.lookupYFile',
    'rangeCoord'      : 'mintpy.load.lookupXFile',
    'incidenceAngle'  : 'mintpy.load.incAngleFile',
    'azimuthAngle'    : 'mintpy.load.azAngleFile',
    'shadowMask'      : 'mintpy.load.shadowMaskFile',
    'waterMask'       : 'mintpy.load.waterMaskFile',
    'bperp'           : 'mintpy.load.bperpFile',
}

DEFAULT_TEMPLATE = """template:
########## 1. Load Data (--load to exit after this step)
{}\n
{}\n
{}""".format(
    auto_path.AUTO_PATH_GAMMA,
    auto_path.AUTO_PATH_ISCE_STRIPMAP,
    auto_path.AUTO_PATH_ISCE_TOPS,
)

TEMPLATE = get_template_content('load_data')

NOTE = """NOTE:
  For interferogram, unwrapPhase is required, the other dataset are optional, including coherence, connectComponent, wrapPhase, etc.
  The unwrapPhase metadata file requires DATE12 attribute in YYMMDD-YYMMDD format.
  All path of data file must contain the reference and secondary date, either in file name or folder name.
"""

EXAMPLE = """example:
  # MUST run in the mintpy working directory!

  # show example template file for ISCE/ROI_PAC/GAMMA products
  load_data.py -H

  # load & write the following HDF5 files:
  # ./inputs/ifgramStack.h5   for interferogram        stack
  # ./inputs/ionStack.h5      for ionosphere           stack
  # ./inputs/offsetStack.h5   for range/azimuth offset stack
  # ./inputs/geometryRadar.h5 for geometry in radar coordinates
  # ./inputs/geometryGeo.h5   for geometry in geo   coordinates
  load_data.py -t smallbaselineApp.cfg
  load_data.py -t smallbaselineApp.cfg GalapagosSenDT128.txt --project GalapagosSenDT128

  # load geometry ONLY
  smallbaselineApp.py SaltonSeaSenDT173.txt -g
  load_data.py -t smallbaselineApp.cfg --geom
"""


def create_parser(subparsers=None):
    """Create command line parser."""
    synopsis = 'Load stacks of interferograms to HDF5 files'
    epilog = TEMPLATE + '\n' + NOTE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('-H', dest='print_example_template', action='store_true',
                        help='Print/Show the example template file for loading.')
    parser.add_argument('-t', '--template', dest='template_file', type=str, nargs='+',
                        help='template file(s) with path info.')
    parser.add_argument('--geom','--geometry', dest='only_load_geometry', action='store_true',
                        help='Load the geometry file(s) ONLY.')

    # options from template file name & content
    parser.add_argument('--project', type=str, dest='PROJECT_NAME',
                        help='project name of dataset for INSARMAPS Web Viewer')
    parser.add_argument('--processor', type=str, dest='processor', choices=PROCESSOR_LIST,
                        help='InSAR processor/software of the file', default='isce')
    parser.add_argument('--enforce', '-f', dest='updateMode', action='store_false',
                        help='Disable the update mode, or skip checking dataset already loaded.')
    parser.add_argument('--compression', choices={'gzip', 'lzf', None}, default=None,
                        help='compress loaded geometry while writing HDF5 file, default: None.')

    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # check --template option
    if inps.template_file:
        pass

    elif inps.print_example_template:
        print(DEFAULT_TEMPLATE)
        sys.exit(0)

    else:
        parser.print_usage()
        print(('{}: error: one of the following arguments are required:'
               ' -t/--template, -H'.format(os.path.basename(__file__))))
        print('{} -H to show the example template file'.format(os.path.basename(__file__)))
        sys.exit(1)

    return inps



#################################################################
def main(iargs=None):
    from ..load_data import (
        get_extra_metadata,
        read_inps2dict,
        read_inps_dict2geometry_dict_object,
        read_subset_box,
        read_inps_dict2ifgram_stack_dict_object,
        prepare_metadata,
        run_or_skip,
    )

    inps = cmd_line_parse(iargs)
    start_time = time.time()

    # read input options
    iDict = read_inps2dict(inps)

    ## prepare metadata
    prepare_metadata(iDict)
    extraDict = get_extra_metadata(iDict)

    # skip data writing for aria as it is included in prep_aria
    if iDict['processor'] == 'aria':
        return

    ## search & write data files
    print('-'*50)
    print('updateMode : {}'.format(iDict['updateMode']))
    print('compression: {}'.format(iDict['compression']))
    print('multilook x/ystep: {}/{}'.format(iDict['xstep'], iDict['ystep']))
    print('multilook method : {}'.format(iDict['method']))
    kwargs = dict(updateMode=iDict['updateMode'], xstep=iDict['xstep'], ystep=iDict['ystep'])

    # read subset info [need the metadata from above]
    iDict = read_subset_box(iDict)

    # geometry in geo / radar coordinates
    geom_dset_name2template_key = {
        **GEOM_DSET_NAME2TEMPLATE_KEY,
        **IFG_DSET_NAME2TEMPLATE_KEY,
        **OFF_DSET_NAME2TEMPLATE_KEY,
    }
    geom_geo_obj, geom_radar_obj = read_inps_dict2geometry_dict_object(iDict, geom_dset_name2template_key)
    geom_geo_file = os.path.abspath('./inputs/geometryGeo.h5')
    geom_radar_file = os.path.abspath('./inputs/geometryRadar.h5')

    if run_or_skip(geom_geo_file, geom_geo_obj, iDict['box4geo'], **kwargs) == 'run':
        geom_geo_obj.write2hdf5(
            outputFile=geom_geo_file,
            access_mode='w',
            box=iDict['box4geo'],
            xstep=iDict['xstep'],
            ystep=iDict['ystep'],
            compression='lzf')

    if run_or_skip(geom_radar_file, geom_radar_obj, iDict['box'], **kwargs) == 'run':
        geom_radar_obj.write2hdf5(
            outputFile=geom_radar_file,
            access_mode='w',
            box=iDict['box'],
            xstep=iDict['xstep'],
            ystep=iDict['ystep'],
            compression='lzf',
            extra_metadata=extraDict)

    # observations: ifgram, ion or offset
    # loop over obs stacks
    stack_ds_name2tmpl_key_list = [
        IFG_DSET_NAME2TEMPLATE_KEY,
        ION_DSET_NAME2TEMPLATE_KEY,
        OFF_DSET_NAME2TEMPLATE_KEY,
    ]
    stack_files = ['ifgramStack.h5', 'ionStack.h5', 'offsetStack.h5']
    stack_files = [os.path.abspath(os.path.join('./inputs', x)) for x in stack_files]
    for ds_name2tmpl_opt, stack_file in zip(stack_ds_name2tmpl_key_list, stack_files):

        # initiate dict objects
        stack_obj = read_inps_dict2ifgram_stack_dict_object(iDict, ds_name2tmpl_opt)

        # use geom_obj as size reference while loading ionosphere
        geom_obj = None
        if os.path.basename(stack_file).startswith('ion'):
            geom_obj = geom_geo_obj if iDict['geocoded'] else geom_radar_obj

        # write dict objects to HDF5 files
        if run_or_skip(stack_file, stack_obj, iDict['box'], geom_obj=geom_obj, **kwargs) == 'run':
            stack_obj.write2hdf5(
                outputFile=stack_file,
                access_mode='w',
                box=iDict['box'],
                xstep=iDict['xstep'],
                ystep=iDict['ystep'],
                mli_method=iDict['method'],
                compression=iDict['compression'],
                extra_metadata=extraDict,
                geom_obj=geom_obj)

    # time info
    m, s = divmod(time.time()-start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.\n'.format(m, s))

    return


#################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
