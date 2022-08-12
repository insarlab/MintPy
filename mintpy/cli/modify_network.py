############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import os
import sys

from mintpy.defaults.template import get_template_content
from mintpy.utils.arg_utils import create_argument_parser


###############################  Usage  ################################
TEMPLATE = get_template_content('modify_network')

REFERENCE = """reference:
  Yunjun, Z., Fattahi, H. and Amelung, F. (2019), Small baseline InSAR time series analysis:
  Unwrapping error correction and noise reduction, Computers & Geosciences, 133, 104331,
  doi:10.1016/j.cageo.2019.104331.

  Chaussard, E., Bürgmann, R., Fattahi, H., Nadeau, R. M., Taira, T., Johnson, C. W. and Johanson, I.
  (2015), Potential for larger earthquakes in the East San Francisco Bay Area due to the direct
  connection between the Hayward and Calaveras Faults, Geophysical Research Letters, 42(8),
  2734-2741, doi:10.1002/2015GL063575.

  Kang, Y., Lu, Z., Zhao, C., Xu, Y., Kim, J. W., & Gallegos, A. J. (2021).InSAR monitoring
  of creeping landslides in mountainous regions: A case study in Eldorado National Forest,
  California. Remote Sensing of Environment, 258, 112400. doi:10.1016/j.rse.2021.112400
"""

EXAMPLE = """example:
  modify_network.py inputs/ifgramStack.h5 -t smallbaselineApp.cfg
  modify_network.py inputs/ifgramStack.h5 --reset
  modify_network.py inputs/ifgramStack.h5 --manual
"""


def create_parser(subparsers=None):
    synopsis = 'Modify the network of interferograms'
    epilog = REFERENCE + '\n' + TEMPLATE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', help='Files to modify/drop network, e.g. inputs/ifgramStack.h5.')
    parser.add_argument('-t', '--template', dest='template_file',
                        help='Template file with input options')
    parser.add_argument('--reset', action='store_true',
                        help='restore all interferograms in the file, by marking all dropIfgram=True')
    parser.add_argument('--noaux', dest='update_aux', action='store_false',
                        help='Do not update auxiliary files, e.g.\n' +
                             'maskConnComp.h5 or avgSpatialCoh.h5 from ifgramStack.h5')

    # 1. temp/perp baseline, num of conn., dates, pair index, etc.
    parser.add_argument('--max-tbase', dest='tempBaseMax',
                        type=float, help='max temporal baseline in days')
    parser.add_argument('--max-pbase', dest='perpBaseMax',
                        type=float, help='max perpendicular baseline in meters')
    parser.add_argument('--max-conn-num', dest='connNumMax', type=int,
                        help='max number of connections/neighbors per acquisition')
    parser.add_argument('-r', '--reference', dest='referenceFile',
                        help='Reference hdf5 / list file with network information.\n'
                             'i.e. ifgramStack.h5, date12_list.txt')
    parser.add_argument('--exclude-ifg-index', dest='excludeIfgIndex', nargs='*',
                        help='index of interferograms to remove/drop.\n1 as the first')
    parser.add_argument('--exclude-date', dest='excludeDate', nargs='*',
                        help='date(s) to remove/drop, all interferograms included date(s) will be removed')
    parser.add_argument('--start-date', '--min-date', dest='startDate',
                        help='remove/drop interferograms with date earlier than start-date in YYMMDD or YYYYMMDD format')
    parser.add_argument('--end-date', '--max-date', dest='endDate',
                        help='remove/drop interferograms with date later than end-date in YYMMDD or YYYYMMDD format')

    # 2. coherence-based network
    cohBased = parser.add_argument_group('Data-driven network modification', 'Drop/modify network based on data')
    # 2.1 coherence-based
    cohBased.add_argument('--coherence-based', dest='coherenceBased', action='store_true',
                          help='Enable coherence-based network modification (default: %(default)s).')
    cohBased.add_argument('--min-coherence', dest='minCoherence', type=float, default=0.7,
                          help='Minimum coherence value (default: %(default)s).')
    # 2.2 area-ratio-based
    cohBased.add_argument('--area-ratio-based', dest='areaRatioBased', action='store_true',
                          help='Enable area ratio-based network modification (default: %(default)s).')
    cohBased.add_argument('--min-area-ratio', dest='minAreaRatio', type=float, default=0.75,
                          help='Minimum area ratio value (default: %(default)s).')
    # common parameters
    cohBased.add_argument('--no-mst', dest='keepMinSpanTree', action='store_false',
                          help='Do not keep interferograms in Min Span Tree network based on inversed mean coherene')
    cohBased.add_argument('--mask', dest='maskFile', default='waterMask.h5',
                          help='Mask file used to calculate the spatial coherence '
                               '(default: waterMask.h5 or None)')
    cohBased.add_argument('--aoi-yx', dest='aoi_pix_box', type=int, nargs=4, metavar=('X0', 'Y0', 'X1', 'Y1'), default=None,
                          help='AOI in row/column range for coherence calculation (default: %(default)s).')
    cohBased.add_argument('--aoi-lalo', dest='aoi_geo_box', type=float, nargs=4, metavar=('W', 'S', 'E', 'N'), default=None,
                          help='AOI in lat/lon range for coherence calculation (default: %(default)s).')
    cohBased.add_argument('--lookup', dest='lookupFile',
                          help='Lookup table/mapping transformation file for geo/radar coordinate conversion.\n' +
                               'Needed for mask AOI in lalo')

    # 3. manual selection
    manual = parser.add_argument_group('Manual Network', 'Manually select/drop/modify network')
    manual.add_argument('--manual', action='store_true',
                        help='display network to manually choose line/interferogram to remove')
    return parser


def cmd_line_parse(iargs=None):
    from ..utils import utils as ut
    from ..modify_network import read_template2inps, read_input_index_list

    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    if not inps.lookupFile:
        inps.lookupFile = ut.get_lookup_file()

    # Convert index : input to continous index list
    if inps.excludeIfgIndex:
        inps.excludeIfgIndex = read_input_index_list(inps.excludeIfgIndex, stackFile=inps.file)
    else:
        inps.excludeIfgIndex = []

    # required input arguments
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)
    elif all(not i for i in [inps.referenceFile, inps.tempBaseMax, inps.perpBaseMax, inps.connNumMax,
                             inps.excludeIfgIndex, inps.excludeDate, inps.coherenceBased, inps.areaRatioBased,
                             inps.startDate, inps.endDate, inps.reset, inps.manual]):
        msg = 'No input option found to remove interferogram, exit.\n'
        msg += 'To manually modify network, please use --manual option '
        raise Exception(msg)

    if not os.path.isfile(inps.maskFile):
        inps.maskFile = None
    return inps


#########################  Main Function  ##############################
def main(iargs=None):
    from ..modify_network import reset_network, get_date12_to_drop
    from ..objects import ifgramStack
    from ..utils import utils as ut

    inps = cmd_line_parse(iargs)

    if inps.reset:
        print('--------------------------------------------------')
        reset_network(inps.file)
        return inps.file

    inps.date12_to_drop = get_date12_to_drop(inps)

    if inps.date12_to_drop is not None:
        ifgramStack(inps.file).update_drop_ifgram(inps.date12_to_drop)
        ut.touch('coherenceSpatialAvg.txt')
        print('Done.')
    return


########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
