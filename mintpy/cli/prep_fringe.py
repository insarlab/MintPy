############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import os
import sys
import glob
from mintpy.utils import arg_utils


####################################################################################
EXAMPLE = """example:
  prep_fringe.py -u './PS_DS/unwrap/*.unw' -c ./PS_DS/tcorr_ds_ps.bin -g ./geometry -m '../reference/IW*.xml' -b ../baselines -o ./mintpy

  cd ~/data/SanAndreasSenDT42/fringe
  prep_fringe.py

  ## example commands after prep_fringe.py
  reference_point.py timeseries.h5 -y 500 -x 1150
  generate_mask.py temporalCoherence.h5 -m 0.7 -o maskTempCoh.h5
  tropo_pyaps3.py -f timeseries.h5 -g inputs/geometryRadar.h5
  remove_ramp.py timeseries_ERA5.h5 -m maskTempCoh.h5 -s linear
  dem_error.py timeseries_ERA5_ramp.h5 -g inputs/geometryRadar.h5
  timeseries2velocity.py timeseries_ERA5_ramp_demErr.h5
  geocode.py velocity.h5 -l inputs/geometryRadar.h5
"""

def create_parser(subparsers=None):
    """Command Line Parser"""
    synopsis = "Prepare FRInGE products for MintPy"
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = arg_utils.create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('-u', '--unw-file', dest='unwFile', type=str, default='./PS_DS/unwrap/*.unw',
                        help='path pattern of unwrapped interferograms (default: %(default)s).')
    parser.add_argument('-c', '--coh-file', dest='cohFile', type=str, default='./PS_DS/tcorr_ds_ps.bin',
                        help='temporal coherence file (default: %(default)s).')
    parser.add_argument('--ps-mask', dest='psMaskFile', type=str, default='./ampDispersion/ps_pixels',
                        help='PS pixels file (default: %(default)s).')
    parser.add_argument('-g', '--geom-dir', dest='geomDir', type=str, default='./geometry',
                        help='FRInGE geometry directory (default: %(default)s).\n'
                             'This is used to grab 1) bounding box\n'
                             '                 AND 2) geometry source directory where the binary files are.')

    parser.add_argument('-m', '--meta-file', dest='metaFile', type=str, default='../reference/IW*.xml',
                        help='metadata file (default: %(default)s).\n'
                             'e.g.: ./reference/IW1.xml        for ISCE/topsStack OR\n'
                             '      ./referenceShelve/data.dat for ISCE/stripmapStack')
    parser.add_argument('-b', '--baseline-dir', dest='baselineDir', type=str, default='../baselines',
                        help='baseline directory (default: %(default)s).')

    parser.add_argument('-o', '--out-dir', dest='outDir', type=str, default='./mintpy',
                        help='output directory (default: %(default)s).')

    parser.add_argument('-r','--range', dest='lks_x', type=int, default=1,
                        help='number of looks in range direction, for multilooking applied after fringe processing.\n'
                             'Only impacts metadata. (default: %(default)s).')
    parser.add_argument('-a','--azimuth', dest='lks_y', type=int, default=1,
                        help='number of looks in azimuth direction, for multilooking applied after fringe processing.\n'
                             'Only impacts metadata. (default: %(default)s).')

    parser.add_argument('--geom-only', action='store_true',
                        help='Only create the geometry file (useful for geocoding a watermask).')

    parser = arg_utils.add_subset_argument(parser, geo=False)

    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # in case meta_file is input as wildcard
    inps.metaFile = sorted(glob.glob(inps.metaFile))[0]

    return inps


####################################################################################
def main(iargs=None):
    from .. import subset
    from ..utils import  isce_utils, utils as ut, attribute as attr
    from ..prep_fringe import (
        read_vrt_info, prepare_metadata,
        prepare_geometry,
        prepare_timeseries,
        prepare_temporal_coherence,
        prepare_ps_mask,
        prepare_stack,
    )
    inps = cmd_line_parse(iargs)

    # translate input options
    processor = isce_utils.get_processor(inps.metaFile)
    src_box, geom_src_dir = read_vrt_info(os.path.join(inps.geomDir, 'lat.vrt'))

    # metadata
    meta = prepare_metadata(inps.metaFile, geom_src_dir, src_box, nlks_x=inps.lks_x, nlks_y=inps.lks_y)


    # subset - read pix_box for fringe file
    pix_box = subset.subset_input_dict2box(vars(inps), meta)[0]
    pix_box = ut.coordinate(meta).check_box_within_data_coverage(pix_box)
    print('input subset in y/x: {}'.format(pix_box))

    # subset - update src_box for isce file and meta
    src_box = (pix_box[0] + src_box[0],
               pix_box[1] + src_box[1],
               pix_box[2] + src_box[0],
               pix_box[3] + src_box[1])
    meta = attr.update_attribute4subset(meta, pix_box)
    print('input subset in y/x with respect to the VRT file: {}'.format(src_box))


    ## output directory
    for dname in [inps.outDir, os.path.join(inps.outDir, 'inputs')]:
        os.makedirs(dname, exist_ok=True)

    ## output filename
    ts_file      = os.path.join(inps.outDir, 'timeseries.h5')
    tcoh_file    = os.path.join(inps.outDir, 'temporalCoherence.h5')
    ps_mask_file = os.path.join(inps.outDir, 'maskPS.h5')
    stack_file   = os.path.join(inps.outDir, 'inputs/ifgramStack.h5')
    if 'Y_FIRST' in meta.keys():
        geom_file = os.path.join(inps.outDir, 'inputs/geometryGeo.h5')
    else:
        geom_file = os.path.join(inps.outDir, 'inputs/geometryRadar.h5')

    ## 1 - geometry (from SLC stacks before fringe, e.g. ISCE2)
    prepare_geometry(
        outfile=geom_file,
        geom_dir=geom_src_dir,
        box=src_box,
        metadata=meta)

    if inps.geom_only:
        return ts_file, tcoh_file, ps_mask_file, geom_file

    ## 2 - time-series (from fringe)
    prepare_timeseries(
        outfile=ts_file,
        unw_file=inps.unwFile,
        metadata=meta,
        processor=processor,
        baseline_dir=inps.baselineDir,
        box=pix_box)

    ## 3 - temporal coherence and mask for PS (from fringe)
    prepare_temporal_coherence(
        outfile=tcoh_file,
        infile=inps.cohFile,
        metadata=meta,
        box=pix_box)

    prepare_ps_mask(
        outfile=ps_mask_file,
        infile=inps.psMaskFile,
        metadata=meta,
        box=pix_box)

    ## 4 - ifgramStack for unwrapped phase and connected components
    prepare_stack(
        outfile=stack_file,
        unw_file=inps.unwFile,
        metadata=meta,
        processor=processor,
        baseline_dir=inps.baselineDir,
        box=pix_box)

    print('Done.')


####################################################################################
if __name__=="__main__":
    main(sys.argv[1:])
