############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import os
import sys
import time
import glob
from mintpy.utils.arg_utils import create_argument_parser


####################################################################################
TEMPLATE = """template options:
  ########## 1. load_data
  ## no   - save   0% disk usage, fast [default]
  ## lzf  - save ~57% disk usage, relative slow
  ## gzip - save ~62% disk usage, very slow [not recommend]
  mintpy.load.processor      = aria  #[isce, aria, snap, gamma, roipac], auto for isce
  mintpy.load.updateMode     = auto  #[yes / no], auto for yes, skip re-loading if HDF5 files are complete
  mintpy.load.compression    = auto  #[gzip / lzf / no], auto for no.
  ##---------interferogram datasets:
  mintpy.load.unwFile        = ../stack/unwrapStack.vrt
  mintpy.load.corFile        = ../stack/cohStack.vrt
  mintpy.load.connCompFile   = ../stack/connCompStack.vrt
  mintpy.load.magFile        = ../stack/ampStack.vrt        # optional
  ##---------geometry datasets:
  mintpy.load.demFile        = ../DEM/SRTM_3arcsec.dem
  mintpy.load.incAngleFile   = ../incidenceAngle/*.vrt
  mintpy.load.azAngleFile    = ../azimuthAngle/*.vrt
  mintpy.load.waterMaskFile  = ../mask/watermask.msk
  ##---------subset (optional):
  ## if both yx and lalo are specified, use lalo option
  mintpy.subset.yx           = auto    #[y0:y1,x0:x1 / no], auto for no
  mintpy.subset.lalo         = auto    #[lat0:lat1,lon0:lon1 / no], auto for no
  ##---------multilook (optional):
  ## multilook while loading data with the specified method, to reduce dataset size
  ## nearest, mean and median methods are applicable to interferogram/ionosphere/offset stack(s), except for:
  ## connected components and all geometry datasets, for which nearest is hardwired.
  mintpy.multilook.method    = auto    #[nearest, mean, median], auto for nearest - lines/rows skipping approach
  mintpy.multilook.ystep     = auto    #[int >= 1], auto for 1 - no multilooking
  mintpy.multilook.xstep     = auto    #[int >= 1], auto for 1 - no multilooking
"""

EXAMPLE = """example:
  prep_aria.py -t smallbaselineApp.cfg    # recommended
  prep_aria.py -t SanFranSenDT42.txt
  prep_aria.py -s ../stack/ -d ../DEM/SRTM_3arcsec.dem -i ../incidenceAngle/*.vrt
  prep_aria.py -s ../stack/ -d ../DEM/SRTM_3arcsec.dem -i ../incidenceAngle/*.vrt -a ../azimuthAngle/*.vrt -w ../mask/watermask.msk

  # download / extract / prepare inteferograms stack from ARIA using ARIA-tools:
  # reference: https://github.com/aria-tools/ARIA-tools
  ariaDownload.py -b '37.25 38.1 -122.6 -121.75' --track 42
  ariaTSsetup.py -f 'products/*.nc' -b '37.25 38.1 -122.6 -121.75' --mask Download --num_threads 4 --verbose
"""

def create_parser(subparsers=None):
    """Command line parser."""
    synopsis = 'Prepare ARIA processed products for MintPy.'
    epilog = TEMPLATE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('-t','--template', dest='template_file', type=str,
                        help='template file with the options')
    parser.add_argument('-o', '--output', type=str, nargs=2, dest='outfile',
                        default=['./inputs/ifgramStack.h5',
                                 './inputs/geometryGeo.h5'],
                        help='output HDF5 file')
    parser.add_argument('--update', dest='updateMode', action='store_true',
                        help='Enable the update mode: checking dataset already loaded.')
    parser.add_argument('--compression', choices={'gzip', 'lzf', None}, default=None,
                        help='HDF5 file compression, default: %(default)s')

    # ifgramStack
    stack = parser.add_argument_group('interferogram stack')
    stack.add_argument('-s','--stack-dir', dest='stackDir', type=str,
                       help='The directory which contains stack VRT files.')
    stack.add_argument('-u','--unwrap-stack-name', dest='unwFile', type=str,
                       default="unwrapStack.vrt",
                       help='Name of the stack VRT file of unwrapped data.\n'+
                            'default: %(default)s')
    stack.add_argument('-c','--coherence-stack-name', dest='corFile', type=str,
                       default="cohStack.vrt",
                       help='Name of the stack VRT file of coherence data.\n'+
                            'default: %(default)s')
    stack.add_argument('-l','--conn-comp-stack-name', dest='connCompFile', type=str,
                       default="connCompStack.vrt",
                       help='Name of the stack VRT file of connected component data.\n' +
                            'default: %(default)s')
    stack.add_argument('--amp-stack-name','--amplitude-stack-name', dest='magFile', type=str,
                       default="ampStack.vrt",
                       help='Name of the stack VRT file of interferogram amplitude data (optional).\n' +
                            'default: %(default)s')

    # geometryGeo
    geom = parser.add_argument_group('geometry')
    geom.add_argument('-d','--dem', dest='demFile', type=str,
                      help='Name of the DEM file')
    geom.add_argument('-i','--incidence-angle', dest='incAngleFile', type=str,
                      help='Name of the incidence angle file')
    geom.add_argument('-a','--az-angle','--azimuth-angle', dest='azAngleFile', type=str,
                      help='Name of the azimuth angle file.')
    geom.add_argument('-w','--water-mask', dest='waterMaskFile', type=str,
                      help='Name of the water mask file')
    return parser


def cmd_line_parse(iargs = None):
    from mintpy.prep_aria import read_template2inps

    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # default multilook options
    iDict = vars(inps)
    iDict['xstep'] = int(iDict.get('xstep', 1))
    iDict['ystep'] = int(iDict.get('ystep', 1))
    iDict['method'] = str(iDict.get('method', 'nearest'))

    # --template
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)
    print('multilook x/ystep: {}/{}'.format(iDict['xstep'], iDict['ystep']))
    print('multilook method : {}'.format(iDict['method']))

    # --stack-dir
    if inps.stackDir is not None:
        inps.stackDir     = os.path.abspath(inps.stackDir)
        inps.corFile      = os.path.join(inps.stackDir, os.path.basename(inps.corFile))
        inps.unwFile      = os.path.join(inps.stackDir, os.path.basename(inps.unwFile))
        inps.connCompFile = os.path.join(inps.stackDir, os.path.basename(inps.connCompFile))

    # translate wildcard path input with search result
    # if not exist, raise error for required datasets
    #               set to None for the other datasets
    ds_keys = [key for key in list(iDict.keys()) if key.endswith('File')]
    required_ds_keys = ['unwFile', 'corFile', 'demFile', 'incAngleFile']

    print('search input data file info:')
    max_digit = max([len(i) for i in ds_keys])
    for key in ds_keys:
        fname = iDict[key]

        # search for wildcard pattern
        if fname:
            fnames = glob.glob(fname)
        else:
            fnames = []

        # user the first element if more than one exist
        if len(fnames) > 0:
            iDict[key] = fnames[0]
            print('{k:<{w}} : {f}'.format(k=key, w=max_digit, f=fnames[0]))

        elif key in required_ds_keys:
            # raise exception if any required DS is missing
            raise SystemExit('ERROR: no file found for {} in input path: "{}"!'.format(key, iDict[key]))

        else:
            iDict[key] = None

    return inps


####################################################################################
def main(iargs=None):
    import numpy as np
    from mintpy.utils import writefile, attribute as attr
    from mintpy.prep_aria import extract_metadata, read_subset_box, run_or_skip, write_ifgram_stack, write_geometry

    inps = cmd_line_parse(iargs)
    start_time = time.time()

    if inps.updateMode:
        print('update mode: ON')
    else:
        print('update mode: OFF')

    # extract metadata
    meta = extract_metadata(inps.unwFile)
    box, meta = read_subset_box(inps.template_file, meta)
    if inps.xstep * inps.ystep > 1:
        meta = attr.update_attribute4multilook(meta, lks_y=inps.ystep, lks_x=inps.xstep)

    length = int(meta["LENGTH"])
    width = int(meta["WIDTH"])
    num_pair = int(meta["NUMBER_OF_PAIRS"])

    # prepare output directory
    out_dir = os.path.dirname(inps.outfile[0])
    os.makedirs(out_dir, exist_ok=True)

    ########## output file 1 - ifgramStack
    # define dataset structure for ifgramStack
    dsNameDict = {
        "date"             : (np.dtype('S8'), (num_pair, 2)),
        "dropIfgram"       : (np.bool_,       (num_pair,)),
        "bperp"            : (np.float32,     (num_pair,)),
        "unwrapPhase"      : (np.float32,     (num_pair, length, width)),
        "coherence"        : (np.float32,     (num_pair, length, width)),
        "connectComponent" : (np.int16,       (num_pair, length, width)),
    }
    if inps.magFile is not None:
        dsNameDict['magnitude'] = (np.float32, (num_pair, length, width))

    if run_or_skip(inps, dsNameDict, out_file=inps.outfile[0]) == 'run':
        # initiate h5 file with defined structure
        meta['FILE_TYPE'] = 'ifgramStack'
        writefile.layout_hdf5(inps.outfile[0], dsNameDict,
                              metadata=meta,
                              compression=inps.compression)

        # write data to h5 file in disk
        write_ifgram_stack(inps.outfile[0],
                           unwStack=inps.unwFile,
                           cohStack=inps.corFile,
                           connCompStack=inps.connCompFile,
                           ampStack=inps.magFile,
                           box=box,
                           xstep=inps.xstep,
                           ystep=inps.ystep,
                           mli_method=inps.method)

    ########## output file 2 - geometryGeo
    # define dataset structure for geometry
    dsNameDict = {
        "height"             : (np.float32, (length, width)),
        "incidenceAngle"     : (np.float32, (length, width)),
        "slantRangeDistance" : (np.float32, (length, width)),
    }
    if inps.azAngleFile is not None:
        dsNameDict["azimuthAngle"] = (np.float32, (length, width))
    if inps.waterMaskFile is not None:
        dsNameDict["waterMask"]    = (np.bool_,   (length, width))

    if run_or_skip(inps, dsNameDict, out_file=inps.outfile[1]) == 'run':
        # initiate h5 file with defined structure
        meta['FILE_TYPE'] = 'geometry'
        writefile.layout_hdf5(inps.outfile[1], dsNameDict,
                              metadata=meta,
                              compression=inps.compression)

        # write data to disk
        write_geometry(inps.outfile[1],
                       demFile=inps.demFile,
                       incAngleFile=inps.incAngleFile,
                       azAngleFile=inps.azAngleFile,
                       waterMaskFile=inps.waterMaskFile,
                       box=box,
                       xstep=inps.xstep,
                       ystep=inps.ystep)

    print('-'*50)

    # time info
    m, s = divmod(time.time()-start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.'.format(m, s))


####################################################################################
if __name__=="__main__":
    main(sys.argv[1:])
